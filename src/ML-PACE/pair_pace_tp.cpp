//
// Created by Yury Lysogorskiy on 01.12.23.
//
#ifdef PACE_TP
#include "pair_pace_tp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"


#include <cstring>
#include <exception>
#include <numeric>
#include "yaml-cpp/yaml.h"

#include "utils_pace.h"

#define PACE_TP_OMP 1

const std::string DEFAULT_INPUT_PREFIX = "serving_default_";


namespace LAMMPS_NS {

    struct ACETPImpl {
        ACETPImpl() : model(nullptr) {}

        ~ACETPImpl() {
            delete model;
        }

        cppflow::model *model;
    };
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairPACETensorPotential::PairPACETensorPotential(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new ACETPImpl;

    scale = nullptr;

    chunksize = 4096;

    data_timer.init();
    tp_timer.init();
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairPACETensorPotential::~PairPACETensorPotential() {
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
void PairPACETensorPotential::allocate() {
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
void PairPACETensorPotential::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style pace/tp", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "ACE potentials require 'metal' units");

    if (comm->me == 0) {
        utils::logmesg(lmp, "ACE/TensorPotential\n");
    }

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else if (strcmp(arg[iarg], "padding") == 0) {
            neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);

            iarg += 2;
        } else
            error->all(FLERR, "Unknown pair_style pace keyword: {}", arg[iarg]);
    }

    do_padding = (neigh_padding_fraction > 0);
    if (do_padding && comm->me == 0)
        utils::logmesg(lmp, "Neighbour padding is ON, padding fraction: {}\n", neigh_padding_fraction);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACETensorPotential::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);

    auto potential_path = utils::get_potential_file_path(arg[2]);

    //load potential file
    delete aceimpl->model;
    //load potential file
    if (comm->me == 0) utils::logmesg(lmp, "Loading {}\n", potential_path);
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

    if (comm->me == 0) {
        utils::logmesg(lmp, "Model loaded\n");
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
            if (comm->me == 0) utils::logmesg(lmp, "Skipping LAMMPS atom type #{}(NULL)\n", i);
        } else {
            int atomic_number = PACE::AtomicNumberByName(elemname);
            if (atomic_number == -1) error->all(FLERR, "'{}' is not a valid element\n", elemname);
            int mu = elements_to_index_map.at(elemname);
            if (mu != -1) {
                if (comm->me == 0)
                    utils::logmesg(lmp, "Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
                                   elemname, mu);
                map[i] = mu;
                // set up LAMMPS atom type to ACE species  mapping for ace evaluator
                element_type_mapping[i] = mu;
            } else {
                error->all(FLERR, "Element {} is not supported by ACE-potential from file {}", elemname,
                           potential_path);
            }
        }
    }

    // initialize scale factor
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) scale[i][j] = 1.0;
    }

}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACETensorPotential::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style pace/tp requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style pace/tp requires newton pair on");

    // request a full neighbor list
    neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACETensorPotential::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];
    //TODO: need info on cutoff
    return cutoff;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairPACETensorPotential::extract(const char *str, int &dim) {
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}


/* ---------------------------------------------------------------------- */

void PairPACETensorPotential::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;

    // number of atoms in cell
    int nlocal = atom->nlocal;
    int n_atoms_extended = atom->nlocal + atom->nghost;
    int n_fake_atoms = (do_padding ? 1 : 0);
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
    std::vector<int32_t> atomic_mu_i_vector(n_atoms_extended + n_fake_atoms, element_type_mapping[type[0]]);
    for (i = 0; i < n_atoms_extended; ++i)
        atomic_mu_i_vector[i] = element_type_mapping[type[i]];
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "atomic_mu_i" + ":0",
                        cppflow::tensor(atomic_mu_i_vector, {n_atoms_extended + n_fake_atoms}));

    // batch_nat = number of extened atoms + padding
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nat" + ":0",
                        cppflow::tensor(std::vector<int32_t>{n_atoms_extended + n_fake_atoms}, {}));

    // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nreal_atoms_per_structure" + ":0",
                        cppflow::tensor(std::vector<int32_t>{n_atoms_extended + n_fake_atoms}, {1}));

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
            utils::logmesg(lmp,
                           "Neighbours padding: new num. of neighbours = {} (+{:.3f}% fake neighbours)\n",
                           tot_neighbours, 100. * (double) n_fake_neighbours / n_real_neighbours);
        }
    } else {
        tot_neighbours = n_real_neighbours;
    }

    std::vector<int32_t> ind_i_vector(tot_neighbours);
    std::vector<int32_t> ind_j_vector(tot_neighbours);
    std::vector<int32_t> mu_i_vector(tot_neighbours);
    std::vector<int32_t> mu_j_vector(tot_neighbours);
#ifdef PACE_TP_OMP
#pragma omp parallel for default(none)  \
    shared(inum, ilist, firstneigh, numneigh, actual_jnum_shift, cutoff_sq, type, x, \
    ind_i_vector, ind_j_vector, mu_i_vector, mu_j_vector) \
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
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            const double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                ind_i_vector[tot_ind] = i;
                ind_j_vector[tot_ind] = j;
                mu_i_vector[tot_ind] = element_type_mapping[type[i]];
                mu_j_vector[tot_ind] = element_type_mapping[type[j]];
                ++tot_ind;
            }
        }
    }

    // add fake bonds
    int fake_atom_ind = n_atoms_extended + n_fake_atoms - 1;
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
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
                        cppflow::tensor(mu_i_vector, {tot_neighbours}));
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                        cppflow::tensor(mu_j_vector, {tot_neighbours}));


    // num_struc: 1
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "num_struc" + ":0",
                        cppflow::tensor(std::vector<int32_t>{1}, {}));


    // positions: [n_extened_atoms+n_fake_atoms, 3]
    std::vector<double> positions_vec(3 * (n_atoms_extended + n_fake_atoms), 0.);
    for (int tot_ind = 0, i = 0; i < n_atoms_extended; ++i) {
        for (int k = 0; k < 3; ++k, ++tot_ind)
            positions_vec[tot_ind] = x[i][k];
    }

    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "positions" + ":0",
                        cppflow::tensor(positions_vec, {n_atoms_extended + n_fake_atoms, 3}));

    data_timer.stop();

    tp_timer.start();
    //CALL MODEL
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            {
                    "StatefulPartitionedCall:0", // energy [nat]
                    "StatefulPartitionedCall:1", // force-per-bond [n_neigh, 3]
            }
    );
    tp_timer.stop();

    data_timer.start();
    auto &e_out = output[0];
    auto &f_out = output[1];

    auto e_tens = e_out.get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));

    auto f_tens = f_out.get_tensor();
    const double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));
    // parallel loop over atoms, does not work because atomic addition on f[j][0..2]+=f_ij, does not work properly...
    // this section remains serial
//#ifdef PACE_TP_OMP
//    #pragma omp parallel for default(none) shared(inum, ilist, firstneigh, numneigh, x, actual_jnum_shift, cutoff_sq, \
//                e_data, f_data, f, nlocal, newton_pair, type)     private(i, jlist, jnum, jj, delx, dely, delz, j, fij) \
//                reduction (+: f[:nlocal][:3])
//#endif
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        const int itype = type[i];

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        int tot_ind = actual_jnum_shift[ii] * 3;
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;

            double rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoff_sq) {
                fij[0] = -scale[itype][itype] * f_data[tot_ind];
                fij[1] = -scale[itype][itype] * f_data[tot_ind + 1];
                fij[2] = -scale[itype][itype] * f_data[tot_ind + 2];

                tot_ind += 3;

                f[i][0] += fij[0];
                f[i][1] += fij[1];
                f[i][2] += fij[2];
                //for OpenMP it should be atomic (critical is too slow), but it does no work with f[j][0..2] array
                f[j][0] -= fij[0];
                f[j][1] -= fij[1];
                f[j][2] -= fij[2];

                // tally per-atom virial contribution, OpenMP critical !
                if (vflag_either) {
                    ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely, -delz);
                }
            }
        }

        // tally energy contribution
        if (eflag_either) {
            // evdwl = energy of atom I
            evdwl = scale[itype][itype] * e_data[i];
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    } // end for(ii)

    if (vflag_fdotr) virial_fdotr_compute();
    data_timer.stop();
    // end modifications YL
}
#endif // #ifdef PACE_TP