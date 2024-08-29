//
// Created by Yury Lysogorskiy on 01.12.23.
//
//#ifdef GRACE

#include "pair_grace.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"

#include <cstring>
#include <exception>
#include <numeric>
#include "yaml-cpp/yaml.h"

#include "ace-evaluator/ace_arraynd.h"

#include "utils_pace.h"

//#define PACE_TP_OMP 1

const std::string DEFAULT_INPUT_PREFIX = "serving_default_";


namespace LAMMPS_NS {

    struct ACETPImpl {
        ACETPImpl() : model(nullptr) {}

        ~ACETPImpl() {
            delete model;
        }

        cppflow::model *model;
    };

    bool check_tf_graph_input_presented(cppflow::model* model, const std::string& op_name) {
        try {
            auto op_shape = model->get_operation_shape(op_name);
            return true;
        } catch (std::runtime_error& exc) {
            return false;
        }
    }


}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairGRACE::PairGRACE(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new ACETPImpl;

    scale = nullptr;

    chunksize = 4096;

    data_timer.init();
    tp_timer.init();

    no_virial_fdotr_compute = 1;
}


/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairGRACE::~PairGRACE() {
    if (copymode) return;

    delete aceimpl;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }
    auto data_t = data_timer.as_microseconds();
    auto tp_t = tp_timer.as_microseconds();
    std::cout << "Data prep. timer: " << data_t << " mcs" << std::endl;
    std::cout << "TP call timer: " << tp_t << " mcs" << std::endl;
    std::cout << "Data prep time fraction: " << ((double) data_t / (data_t + tp_t) * 1e2) << " %" << std::endl;
}

/* ---------------------------------------------------------------------- */
void PairGRACE::allocate() {
    allocated = 1;
    int n = atom->ntypes + 1;

    memory->create(setflag, n, n, "pair:setflag");
    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(scale, n, n, "pair:scale");
    map = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairGRACE::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "GRACE potentials require 'metal' units");

    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE]\n");
    }

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else if (strcmp(arg[iarg], "padding") == 0) {
            neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);

            iarg += 2;
        } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
            pad_verbose = true;
            iarg += 1;
        } else
            error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
    }

    do_padding = (neigh_padding_fraction > 0);
    if (do_padding && comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Neighbour padding is ON, padding fraction: {}\n", neigh_padding_fraction);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGRACE::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);

    auto potential_path = utils::get_potential_file_path(arg[2]);

    //load potential file
    delete aceimpl->model;
    //load potential file
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);
    // load cppflow model
    aceimpl->model = new cppflow::model(potential_path);

    // read elements from metadata.yaml
    YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
    auto elements_yaml = metadata_yaml["chemical_symbols"];
    elements_name = elements_yaml.as<std::vector<std::string>>();
    nelements = (int) elements_name.size();
    for (int mu = 0; mu < nelements; mu++) {
        elements_to_index_map[elements_name.at(mu)] = mu;
    }
    cutoff = metadata_yaml["cutoff"].as<double>();
    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE] Model loaded\n");
    }


    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    const int n = atom->ntypes;
    element_type_mapping.resize(n + 1);
    // elements to species-type map
    for (int i = 1; i <= n; i++) {
        char *elemname = arg[2 + i];
        if (strcmp(elemname, "NULL") == 0) {
            // species_type=-1 value will not reach ACE Evaluator::compute_atom,
            // but if it will ,then error will be thrown there
            element_type_mapping[i] = -1;
            map[i] = -1;
            if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Skipping LAMMPS atom type #{}(NULL)\n", i);
        } else {
            int atomic_number = PACE::AtomicNumberByName(elemname);
            if (atomic_number == -1) error->all(FLERR, "[GRACE] '{}' is not a valid element\n", elemname);
            int mu = elements_to_index_map.at(elemname);
            if (mu != -1) {
                if (comm->me == 0)
                    utils::logmesg(lmp, "[GRACE] Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                                   elemname, mu);
                map[i] = mu;
                // set up LAMMPS atom type to ACE species  mapping for ace evaluator
                element_type_mapping[i] = mu;
            } else {
                error->all(FLERR, "[GRACE] Element {} is not supported by ACE-potential from file {}", elemname,
                           potential_path);
            }
        }
    }

    // initialize scale factor
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) scale[i][j] = 1.0;
    }

//    auto operations_vec = aceimpl->model->get_operations();
//    std::cout<<"List of operations "<<std::endl;
//    for(const auto& op_name: operations_vec){
//        std::cout<<"Operation `"<<op_name<<"`"<<endl;
//    }

    //
    has_map_atoms_to_structure_op = check_tf_graph_input_presented(aceimpl->model,
                                                                   "serving_default_map_atoms_to_structure");
//    std::cout<<"[DEBUG] has_map_atoms_to_structure_op="<<has_map_atoms_to_structure_op<<endl;

    has_nstruct_total_op = check_tf_graph_input_presented(aceimpl->model,
                                                                   "serving_default_n_struct_total");
//    std::cout<<"[DEBUG] has_nstruct_total_op="<<has_nstruct_total_op<<endl;


    has_mu_i_op = check_tf_graph_input_presented(aceimpl->model,
                                                          "serving_default_mu_i");
//    std::cout<<"[DEBUG] has_mu_i_op="<<has_mu_i_op<<endl;
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACE::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");

    // request a full neighbor list
    neighbor->add_request(this, NeighConst::REQ_FULL);

    // request atom map (maybe?)
    if (atom->map_style == Atom::MAP_NONE) {
        atom->map_init();
        atom->map_set();
    }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGRACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    return cutoff;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairGRACE::extract(const char *str, int &dim) {
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}


/* ---------------------------------------------------------------------- */
/**
 * signature_def['serving_default']:
  The given SavedModel SignatureDef contains the following input(s):
    +inputs['atomic_mu_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_atomic_mu_i:0
    +inputs['batch_tot_nat'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat:0
    +inputs['batch_tot_nat_real'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat_real:0
    +inputs['bond_vector'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: serving_default_bond_vector:0
    +inputs['ind_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_i:0
    +inputs['ind_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_j:0
    +inputs['map_atoms_to_structure'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_map_atoms_to_structure:0
    +inputs['mu_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_mu_j:0
    +inputs['n_struct_total'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_n_struct_total:0


  The given SavedModel SignatureDef contains the following output(s):
    outputs['atomic_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall:0
    outputs['total_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall:1
    outputs['total_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:2
    outputs['virial'] tensor_info:
        dtype: DT_DOUBLE
        shape: (6)
        name: StatefulPartitionedCall:3
  Method name is: tensorflow/serving/predict
 */
void PairGRACE::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
//    tagint *tag = atom->tag;

    // number of atoms in cell
    int nlocal = atom->nlocal;
//    int n_atoms_extended = atom->nlocal + atom->nghost;
    int n_fake_atoms = 0; // no fake atoms needed. fake bonds will be 1e6
    int n_real_neighbours;
    int n_fake_neighbours;

    int newton_pair = force->newton_pair;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    data_timer.start();
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    // atomic_mu_i: per-atom species type + padding with type[0]
    std::vector<int32_t> atomic_mu_i_vector(nlocal + n_fake_atoms, element_type_mapping[type[0]]);
    for (i = 0; i < nlocal; ++i)
        atomic_mu_i_vector[i] = element_type_mapping[type[i]];

    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "atomic_mu_i" + ":0",
                        cppflow::tensor(atomic_mu_i_vector, {nlocal + n_fake_atoms}));

    //    map_atoms_to_structure
    if(has_map_atoms_to_structure_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "map_atoms_to_structure" + ":0",
                            cppflow::tensor(std::vector<int32_t>(nlocal + n_fake_atoms, 0), {nlocal + n_fake_atoms}));
    }

    // batch_nat = number of extened atoms + padding
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_tot_nat" + ":0",
                        cppflow::tensor(std::vector<int32_t>{nlocal + n_fake_atoms}, {}));

    // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_tot_nat_real" + ":0",
                        cppflow::tensor(std::vector<int32_t>{nlocal}, {}));

    // ind_i, ind_j: bonds
    //determine the maximum number of neighbours (within cutoff)
    std::vector<int> actual_jnum(inum, 0);
    double cutoff_sq = cutoff * cutoff;
#ifdef PACE_TP_OMP
#pragma omp parallel for default(none) shared(inum, ilist, jlist, cutoff_sq, numneigh, firstneigh, actual_jnum, x) \
                         private(i, jnum, jj, j, delx, dely, delz)
#endif
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        double xtmp = x[i][0];
        double ytmp = x[i][1];
        double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        int cur_actual_jnum = 0;
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                cur_actual_jnum += 1;
            }
        }
        actual_jnum[ii] = cur_actual_jnum;
    }
    n_real_neighbours = std::accumulate(actual_jnum.begin(), actual_jnum.end(), 0);

    std::vector<int> actual_jnum_shift(actual_jnum.size(), 0);
    std::partial_sum(actual_jnum.begin(), actual_jnum.end() - 1, actual_jnum_shift.begin() + 1);

    if (do_padding) {
        if (n_real_neighbours > tot_neighbours) {
            n_fake_neighbours = static_cast<int>(std::round(n_real_neighbours * neigh_padding_fraction));
            n_fake_neighbours = std::max(n_fake_neighbours, 1);
            tot_neighbours = n_real_neighbours + n_fake_neighbours; // add fake neighbours
            if (pad_verbose)
                utils::logmesg(lmp,
                               "[GRACE] Neighbours padding: new num. of neighbours = {} (+{:.3f}% fake neighbours)\n",
                               tot_neighbours, 100. * (double) n_fake_neighbours / n_real_neighbours);
        }
    } else {
        tot_neighbours = n_real_neighbours;
    }

    std::vector<int32_t> ind_i_vector(tot_neighbours);
    std::vector<int32_t> ind_j_vector(tot_neighbours);
    std::vector<int32_t> mu_i_vector(tot_neighbours);
    std::vector<int32_t> mu_j_vector(tot_neighbours);
    std::vector<double> bond_vector(tot_neighbours *3, 1e6);

#ifdef PACE_TP_OMP
#pragma omp parallel for default(none)  \
    shared(inum, ilist, firstneigh, numneigh, actual_jnum_shift, cutoff_sq, type, x, \
    ind_i_vector, ind_j_vector, mu_i_vector, mu_j_vector, bond_vector) \
    private(i, jlist, jnum, jj, j, delx, dely, delz)
#endif
    for (ii = 0; ii < inum; ++ii) {
        i = ilist[ii];
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        int tot_ind = actual_jnum_shift[ii];
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            const double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                ind_i_vector[tot_ind] = i;
                // remap j to j_local
                int j_local = atom->map(atom->tag[j]);
                ind_j_vector[tot_ind] = j_local;
                mu_i_vector[tot_ind] = element_type_mapping[type[i]];
                mu_j_vector[tot_ind] = element_type_mapping[type[j]];
                double bondx = atom->x[j][0] - atom->x[i][0];
                double bondy = atom->x[j][1] - atom->x[i][1];
                double bondz = atom->x[j][2] - atom->x[i][2];


                bond_vector[3 * tot_ind + 0] = bondx;
                bond_vector[3 * tot_ind + 1] = bondy;
                bond_vector[3 * tot_ind + 2] = bondz;
                ++tot_ind;
            }
        }
    }

    // add fake bonds
    int fake_atom_ind = nlocal + n_fake_atoms - 1;
#ifdef PACE_TP_OMP
#pragma omp parallel for default(none) shared(n_real_neighbours, tot_neighbours, fake_atom_ind, \
        ind_i_vector, ind_j_vector, mu_i_vector, mu_j_vector)
#endif
    for (int tot_ind = n_real_neighbours; tot_ind < tot_neighbours; tot_ind++) {
        ind_i_vector[tot_ind] = fake_atom_ind; // fake atom ind
        ind_j_vector[tot_ind] = fake_atom_ind; // fake atom ind
        mu_i_vector[tot_ind] = 0;
        mu_j_vector[tot_ind] = 0;
    }

    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_i" + ":0",
                        cppflow::tensor(ind_i_vector, {tot_neighbours}));
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_j" + ":0",
                        cppflow::tensor(ind_j_vector, {tot_neighbours}));


    // mu_i, mu_j: bonds
    if (has_mu_i_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
                            cppflow::tensor(mu_i_vector, {tot_neighbours}));
    }
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                        cppflow::tensor(mu_j_vector, {tot_neighbours}));


    // num_struc: 1
    if (has_nstruct_total_op) {
        inputs.emplace_back(DEFAULT_INPUT_PREFIX + "n_struct_total" + ":0",
                            cppflow::tensor(std::vector<int32_t>{1}, {}));
    }

    // vector_offsets: 1
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "bond_vector" + ":0",
                        cppflow::tensor(bond_vector, {tot_neighbours, 3}));


    data_timer.stop();
    tp_timer.start();
    //CALL MODEL
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            {
                    "StatefulPartitionedCall:0", // atomic_energy [nat,1]
                    "StatefulPartitionedCall:1", // total_energy [-1, 1]
                    "StatefulPartitionedCall:2", // total_f [n_at, 3]
                    "StatefulPartitionedCall:3", // virial [6]
            }
    );
    tp_timer.stop();
//    std::cout << "Ave.timing: " << (double) tp_timer.as_microseconds() / nlocal << " mcs/at" << std::endl;

    data_timer.start();
    auto &e_out = output[0]; // atomic_energy
    auto e_tens = e_out.get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));


//    auto &te_out = output[1]; // total_energy


    auto &total_f_out = output[2]; // total_f
    auto total_f_tens = total_f_out.get_tensor();
    const double *total_f_data = static_cast<double *>(TF_TensorData(total_f_tens.get()));


    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        const int itype = type[i];
        double fx = total_f_data[ii * 3 + 0];
        double fy = total_f_data[ii * 3 + 1];
        double fz = total_f_data[ii * 3 + 2];

        f[i][0] += scale[itype][itype] * fx;
        f[i][1] += scale[itype][itype] * fy;
        f[i][2] += scale[itype][itype] * fz;

        // tally energy contribution
        if (eflag_either) {
            // evdwl = energy of atom I
            evdwl = scale[itype][itype] * e_data[i];
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    } // end for(ii)


    // virial order, seems to be OK
    if (vflag_global) {
        auto &v_out = output[3]; // virial
        auto v_tens = v_out.get_tensor();
        const double *v_data = static_cast<double *>(TF_TensorData(v_tens.get()));

//            ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely, -delz);
        virial[0] += v_data[0];
        virial[1] += v_data[1];
        virial[2] += v_data[2];
        virial[3] += v_data[3];
        virial[4] += v_data[4];
        virial[5] += v_data[5];
    }

    data_timer.stop();
    // end modifications YL
}

//#endif // #ifdef GRACE