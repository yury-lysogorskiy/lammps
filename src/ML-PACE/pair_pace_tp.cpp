//
// Created by lysogy36 on 01.12.23.
//

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

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else
            error->all(FLERR, "Unknown pair_style pace keyword: {}", arg[iarg]);
    }

    if (comm->me == 0) {
        utils::logmesg(lmp, "ACE/TensorPotential");
    }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACETensorPotential::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);

    auto potential_file_name = utils::get_potential_file_path(arg[2]);

    //load potential file
    delete aceimpl->model;
    //load potential file
    if (comm->me == 0) utils::logmesg(lmp, "Loading {}\n", potential_file_name);
    // TODO: load cppflow model
    aceimpl->model = new cppflow::model(potential_file_name);

    if (comm->me == 0) {
        utils::logmesg(lmp, "Model loaded\n");
    }

    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    // TODO: aceimpl->ace->element_type_mapping.init(atom->ntypes + 1);

    const int n = atom->ntypes;
    //TODO: elements to species-type map
//    for (int i = 1; i <= n; i++) {
//        char *elemname = arg[2 + i];
//        if (strcmp(elemname, "NULL") == 0) {
//            // species_type=-1 value will not reach ACE Evaluator::compute_atom,
//            // but if it will ,then error will be thrown there
//            aceimpl->ace->element_type_mapping(i) = -1;
//            map[i] = -1;
//            if (comm->me == 0) utils::logmesg(lmp, "Skipping LAMMPS atom type #{}(NULL)\n", i);
//        } else {
//            int atomic_number = AtomicNumberByName_pace(elemname);
//            if (atomic_number == -1) error->all(FLERR, "'{}' is not a valid element\n", elemname);
//            SPECIES_TYPE mu = aceimpl->basis_set->get_species_index_by_name(elemname);
//            if (mu != -1) {
//                if (comm->me == 0)
//                    utils::logmesg(lmp, "Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n", i,
//                                   elemname, mu);
//                map[i] = mu;
//                // set up LAMMPS atom type to ACE species  mapping for ace evaluator
//                aceimpl->ace->element_type_mapping(i) = mu;
//            } else {
//                error->all(FLERR, "Element {} is not supported by ACE-potential from file {}", elemname,
//                           potential_file_name);
//            }
//        }
//    }

    // initialize scale factor
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) scale[i][j] = 1.0;
    }

//    aceimpl->ace->set_basis(*aceimpl->basis_set, 1);
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
//    n_nodes = atom->nlocal + atom->nghost;
    int n_atoms_extended = atom->nlocal + atom->nghost;


    int newton_pair = force->newton_pair;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;


    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    // atomic_mu_i: per-atom species type (+ padding?)
    std::vector<int32_t> atomic_mu_i_vector(type, type + n_atoms_extended);
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "atomic_mu_i" + ":0",
                        cppflow::tensor(atomic_mu_i_vector, {n_atoms_extended}));

    // batch_nat: number of extened atoms (+ padding?)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nat" + ":0",
                        cppflow::tensor(std::vector<int32_t>{n_atoms_extended}, {}));

    // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nreal_atoms_per_structure" + ":0",
                        cppflow::tensor(std::vector<int32_t>{n_atoms_extended}, {1}));

    // ind_i, ind_j: bonds
    //determine the maximum number of neighbours (within cutoff)
    int tot_number_of_neighbours = 0;
    double cutoff_sq = cutoff * cutoff;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
//        double xtmp = atom->x[i][0];
//        double ytmp = atom->x[i][1];
//        double ztmp = atom->x[i][2];
//        jlist = list->firstneigh[i];
        jnum = list->numneigh[i];
        tot_number_of_neighbours += jnum;
//        for (jj=0; jj<jnum; ++jj) {
//            j = jlist[jj];
//            j &= NEIGHMASK;
//            delx = xtmp - atom->x[j][0];
//            dely = ytmp - atom->x[j][1];
//            delz = ztmp - atom->x[j][2];
//            double rsq = delx * delx + dely * dely + delz * delz;
//            if (rsq < cutoff_sq) {
//                tot_number_of_neighbours += 1;
//            }
//        }
    }

    std::vector<int32_t> ind_i_vector(tot_number_of_neighbours);
    std::vector<int32_t> ind_j_vector(tot_number_of_neighbours);
    std::vector<int32_t> mu_i_vector(tot_number_of_neighbours);
    std::vector<int32_t> mu_j_vector(tot_number_of_neighbours);
    int tot_ind = 0;
    for (ii = 0; ii < inum; ++ii) {
        i = ilist[ii];
        jlist = list->firstneigh[i];
        jnum = list->numneigh[i];
        for (jj = 0; jj < jnum; ++jj, ++tot_ind) {
            j = jlist[jj];
            ind_i_vector[tot_ind] = i;
            ind_j_vector[tot_ind] = j;
            mu_i_vector[tot_ind] = type[i];
            mu_j_vector[tot_ind] = type[j];
        }
    }
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_i" + ":0",
                        cppflow::tensor(ind_i_vector, {tot_number_of_neighbours}));
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_j" + ":0",
                        cppflow::tensor(ind_j_vector, {tot_number_of_neighbours}));


    // mu_i, mu_j: bonds
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
                        cppflow::tensor(mu_i_vector, {tot_number_of_neighbours}));
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                        cppflow::tensor(mu_j_vector, {tot_number_of_neighbours}));


    // mu_ij: bonds - stub
    std::vector<int32_t> mu_ij_vector(tot_number_of_neighbours, 0);
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_ij" + ":0",
                        cppflow::tensor(mu_ij_vector, {tot_number_of_neighbours}));

    // num_struc: 1
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "num_struc" + ":0", cppflow::tensor(std::vector<int32_t>{1}, {}));


    // positions: [n_extened_atoms, 3]
    std::vector<double> positions_vec(3 * n_atoms_extended);
    for (i = 0, tot_ind = 0; i < n_atoms_extended; ++i) {
        for (int k = 0; k < 3; ++k, ++tot_ind)
            positions_vec[tot_ind] = x[i][k];
    }
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "positions" + ":0",
                        cppflow::tensor(positions_vec, {n_atoms_extended, 3}));


    // slice_mu_ij: MOCK
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "slice_mu_ij" + ":0", cppflow::tensor(std::vector<int32_t>{0}, {1}));

    //vector_offsets: (-1,3)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "vector_offsets" + ":0",
                        cppflow::tensor(std::vector<double>(3 * tot_number_of_neighbours, 0),
                                        {tot_number_of_neighbours, 3}));

    //CALL MODEL
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            {
                    "StatefulPartitionedCall:0",
                    "StatefulPartitionedCall:1",
            }
    );

    auto &e_out = output[0];
    auto &f_out = output[1];

    auto e_tens = e_out.get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));

    auto f_tens = f_out.get_tensor();
    const double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));

    //loop over atoms
    tot_ind = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        jlist = list->firstneigh[i];
        jnum = list->numneigh[i];

        const int itype = type[i];

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        // 'compute_atom' will update the `aceimpl->ace->e_atom` and `aceimpl->ace->neighbours_forces(jj, alpha)` arrays

        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;

//            fij[0] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 0);
//            fij[1] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 1);
//            fij[2] = scale[itype][itype] * aceimpl->ace->neighbours_forces(jj, 2);
            fij[0] = -f_data[tot_ind];
            fij[1] = -f_data[tot_ind+1];
            fij[2] = -f_data[tot_ind+2];

            tot_ind += 3;

            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution
            if (vflag_either)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely,
                             -delz);
        }

        // tally energy contribution
        if (eflag_either) {
            // evdwl = energy of atom I
            evdwl = scale[itype][itype] * e_data[i];
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }

    if (vflag_fdotr) virial_fdotr_compute();

    // end modifications YL
}