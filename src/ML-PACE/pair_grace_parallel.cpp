//
// Created by Yury Lysogorskiy on 01.12.23.
//
#ifndef NO_GRACE_TF
#define GRACE_PRINT_DEBUG

#include "pair_grace_parallel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include <unistd.h>

#include <numeric>
#include "yaml-cpp/yaml.h"


#include "utils_pace.h"
#include "utils_grace.h"
// CppFlow headers
#include <cppflow/ops.h>
#include <cppflow/model.h>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>
#include <iomanip>


namespace GRACEParallel {

    class GracePaddingDimension {
        public:
            int current_padded_size = 0;
            int num_of_reductions = 0;

            // Settings
            double padding_fraction = 0.01;           // How much to grow
            double reduction_threshold_fraction = 0.2; // Waste limit before reduction
            int max_reductions = 10;                  // -1 for unlimited
            bool enabled = true;
            bool verbose = false;

            GracePaddingDimension() = default;

            /**
             * Updates the padded size based on the current real count.
             * @return The new (or existing) padded size.
             */
            int update(int real_count) {
                if (!enabled) {
                    current_padded_size = real_count;
                    return current_padded_size;
                }

                // 1. Try to find a previously used size that fits
                auto it = padding_history.upper_bound(real_count);

                if (it == padding_history.end()) {
                    // GROW: No suitable size found in history
                    int extra = static_cast<int>(std::round(real_count * padding_fraction));
                    current_padded_size = real_count + std::max(extra, 1);
                    padding_history.insert(current_padded_size);
                    was_updated = true;
                } else {
                    // REUSE: Found a size in history
                    current_padded_size = *it;
                    was_updated = false;

                    // REDUCE: Check if the smallest historical size is too wasteful
                    if (it == padding_history.begin() && can_reduce()) {
                        double waste = (real_count > 0) ?
                                       static_cast<double>(current_padded_size - real_count) / real_count : 0;

                        if (waste > reduction_threshold_fraction) {
                            current_padded_size = real_count;
                            padding_history.insert(current_padded_size);
                            num_of_reductions++;
                            was_updated = true;
                        }
                    }
                }
                return current_padded_size;
            }

            bool last_update_triggered_resize() const { return was_updated; }

            void reset() {
                padding_history.clear();
                num_of_reductions = 0;
                current_padded_size = 0;
            }

        private:
            std::set<int> padding_history;
            bool was_updated = false;

            bool can_reduce() const {
                return (max_reductions == -1 || num_of_reductions < max_reductions);
            }
    };

}    // namespace GRACEParallel

namespace LAMMPS_NS {
    struct GRACEImpl {
        GRACEImpl() : model(nullptr) {}

        ~GRACEImpl() {
            delete model;
        }

        cppflow::model *model;

        // Use the new clean class for different quantities
        GRACEParallel::GracePaddingDimension atom_padding;
        GRACEParallel::GracePaddingDimension neighbor_padding;

        GRACEParallel::GracePaddingDimension atom_padding_1;
        GRACEParallel::GracePaddingDimension neighbor_padding_1;

        // Buffers for reuse across timesteps
        std::vector<int32_t> atomic_mu_i_vector;
        std::vector<int32_t> atomic_mu_i_1_vector;
        std::vector<int32_t> ind_i_vector;
        std::vector<int32_t> ind_j_vector;
        std::vector<int32_t> mu_i_vector;
        std::vector<int32_t> mu_j_vector;
        std::vector<double> bond_vector;
        std::vector<int32_t> ind_i_vector_1;
        std::vector<int32_t> ind_j_vector_1;
    };
}



using namespace LAMMPS_NS;
using namespace GRACEParallel;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairGRACEParallel::PairGRACEParallel(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new GRACEImpl;

    scale = nullptr;

    chunksize = 4096;

    data_timer.init();
    tp_timer.init();

    no_virial_fdotr_compute = 1;
}


/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairGRACEParallel::~PairGRACEParallel() {
    if (copymode) return;

    delete aceimpl;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }
    auto data_t = static_cast<double>(data_timer.as_microseconds());
    auto tp_t = static_cast<double>(tp_timer.as_microseconds());

    utils::logmesg(lmp,
                   "[GRACE:debug, proc #{:d}]: Data preparation timer: {:g} mcs, graph execution time: {:g} mcs, data preparation time fraction: {:.2f} %\n",
                   comm->me, data_t, tp_t, (data_t / (data_t + tp_t) * 1e2));
}

/* ---------------------------------------------------------------------- */
void PairGRACEParallel::allocate() {
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
void PairGRACEParallel::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace/parallel", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "GRACE potentials require 'metal' units");

    auto tf_version = TF_Version();
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", tf_version);

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "chunksize") == 0) {
            chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
        } else if (strcmp(arg[iarg], "padding") == 0) {
            auto padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            aceimpl->atom_padding.padding_fraction = padding_fraction;
            iarg += 2;

        } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
            aceimpl->atom_padding.verbose = true;
            iarg += 1;
        } else  if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
            auto max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            aceimpl->atom_padding.max_reductions = max_number_of_reduction;
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Maximum number of recompilation during padding reduction: {}\n",
                               max_number_of_reduction);
        } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
            auto reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            aceimpl->atom_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Reducing padding fraction: {}\n", reducing_neigh_padding_fraction);
        } else
            error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
    }

    aceimpl->atom_padding.enabled = (aceimpl->atom_padding.padding_fraction > 0);

    if(aceimpl->atom_padding.enabled && comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Neighbour padding is ON, padding fraction: {}, max padding fraction before reduction: {}, max number of reduction(s): {}\n",
                       aceimpl->atom_padding.padding_fraction,
                       aceimpl->atom_padding.reduction_threshold_fraction,
                       aceimpl->atom_padding.max_reductions);



    centroidstressflag = CENTROID_AVAIL;

    aceimpl->neighbor_padding = aceimpl->atom_padding;
    aceimpl->atom_padding_1 = aceimpl->atom_padding;
    aceimpl->neighbor_padding_1 = aceimpl->atom_padding;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGRACEParallel::coeff(int narg, char **arg) {

    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);
    auto potential_path = std::string(arg[2]);

    //load potential file
    delete aceimpl->model;
    //load potential file
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);
    // load cppflow model
    const std::vector<uint8_t> config_bytes = {
        0x32, 0x05,             // Field 6 (GPUOptions), Length 5
        0x82, 0x01, 0x02,       // Field 16 (Experimental), Length 2
        0x18, 0x00              // Field 3 (TF32 Enabled), Value 0 (False)
      };

    aceimpl->model = new cppflow::model(potential_path, config_bytes);
#ifdef GRACE_PRINT_DEBUG
    aceimpl->model->print_signatures();
#endif
    // read elements from metadata.yaml
    YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
    auto elements_yaml = metadata_yaml["chemical_symbols"];
    elements_name = elements_yaml.as<std::vector<std::string>>();
    nelements = static_cast<int>(elements_name.size());
    for (int mu = 0; mu < nelements; mu++) {
        elements_to_index_map[elements_name.at(mu)] = mu;
    }
    cutoff = metadata_yaml["cutoff"].as<double>();

    if (metadata_yaml["cutoff_matrix"]) {
        cutoff_matrix = metadata_yaml["cutoff_matrix"].as<vector<vector<double>>>();
        // assert square size of matrix
        if (cutoff_matrix.size() != nelements)
            error->all(FLERR,
                       "[GRACE] cutoff_matrix is provided, but it's size ({}) is not equal to number of elements ({})\n",
                       cutoff_matrix.size(), nelements);

        for (const auto &v: cutoff_matrix)
            if (v.size() != nelements)
                error->all(FLERR,
                           "[GRACE] cutoff_matrix is provided, but it's row size ({}) is not equal to number of elements ({})\n",
                           v.size(), nelements);

        if (comm->me == 0)
            utils::logmesg(lmp, "[GRACE] Custom cutoff matrix is loaded\n");
        is_custom_cutoffs = true;
    }
    if (comm->me == 0) {
        utils::logmesg(lmp, "[GRACE] Model loaded\n");
    }


    // read args that map atom types to PACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    const int ntypes = atom->ntypes;
    element_type_mapping.resize(ntypes + 1);
    // elements to species-type map
    for (int i = 1; i <= ntypes; i++) {
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
    for (int i = 1; i <= ntypes; i++) {
        for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;
    }

    if (is_custom_cutoffs) {
        // matrix of size [ntypes+1][ntypes+1]
        cutoff_matrix_per_lammps_type.resize(ntypes + 1, vector<double>(ntypes + 1));
        double min_cutoff = 1e99, max_cutoff = 0;
        for (int i = 1; i <= ntypes; i++) {
            for (int j = 1; j <= ntypes; j++) {
                auto val = cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
                cutoff_matrix_per_lammps_type[i][j] = val;
                if (val < min_cutoff) min_cutoff = val;
                if (val > max_cutoff) max_cutoff = val;
            }
        }

        if (comm->me == 0)
            utils::logmesg(lmp, "[GRACE] Custom cutoffs: min={}, max={}\n", min_cutoff, max_cutoff);

    }

    if (aceimpl->model->has_signature("parallel_compute")) {
        this->compute_function_name = "parallel_compute";
    } else  {
        throw std::runtime_error("No 'parallel_compute' function found in SavedModel");
    }
    this->DEFAULT_INPUT_PREFIX = this->compute_function_name+"_";

    //
    has_map_atoms_to_structure_op = aceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX+"map_atoms_to_structure");

    // if (comm->me == 0)
    //     utils::logmesg(lmp, "[GRACE/DEBUG] has_map_atoms_to_structure_op={}\n", has_map_atoms_to_structure_op);

    has_nstruct_total_op = aceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX+"n_struct_total");

    // if (comm->me == 0)
    //     utils::logmesg(lmp, "[GRACE/DEBUG] has_nstruct_total_op={}\n", has_nstruct_total_op);

    has_mu_i_op = aceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX+"mu_i");

    // if (comm->me == 0)
    //     utils::logmesg(lmp, "[GRACE/DEBUG] has_mu_i_op={}\n", has_mu_i_op);

    has_batch_tot_nat = aceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX+"batch_tot_nat");

    // if (comm->me == 0)
    //     utils::logmesg(lmp, "[GRACE/DEBUG] has_batch_tot_nat={}\n", has_batch_tot_nat);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACEParallel::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");

    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);

    // request atom map (maybe?)
    if (atom->map_style == Atom::MAP_NONE) {
        atom->map_init();
        atom->map_set();
    }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGRACEParallel::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    //cutoff from the basis set's radial functions settings
    scale[j][i] = scale[i][j];

    double rcut = (is_custom_cutoffs) ?  cutoff_matrix_per_lammps_type[i][j] : cutoff;

    if (comm->cutghostuser < 2*rcut)
        comm->cutghostuser = 2*rcut;

    return rcut;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairGRACEParallel::extract(const char *str, int &dim) {
    dim = 2;
    if (strcmp(str, "scale") == 0) return (void *) scale;
    return nullptr;
}


std::vector<int> PairGRACEParallel::get_atom_shell_mapping() const {
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    // Initialize all to -1 (OUT)
    std::vector<int> atom_type_map(nlocal + nghost, -1);

    // 1. Mark Real Atoms
    // Real atoms are ALWAYS indices 0 to nlocal-1.
    // We don't need ilist for this; ilist is for neighbor access.
    for (int i = 0; i < nlocal; ++i) {
        atom_type_map[i] = 0;
    }

    // 2. Mark Shell 1 (Neighbors of Real atoms that are ghosts)
    for (int ii = 0; ii < inum; ++ii) {
        int i = ilist[ii];

        // Only process neighbors of REAL atoms to find Shell 1
        if (i >= nlocal) continue;

        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; ++jj) {
            int j = jlist[jj] & NEIGHMASK;
            if (j >= nlocal && atom_type_map[j] == -1) {
                atom_type_map[j] = 1; // Mark as Shell 1
            }
        }
    }

    // 3. Mark Shell 2 (Neighbors of Shell 1 atoms that are neither Real nor Shell 1)
    for (int ii = nlocal; ii < nlocal+nghost; ++ii) {
        int i = ilist[ii];

        // Only process neighbors of SHELL 1 atoms to find Shell 2
        if (atom_type_map[i] == -1)
            atom_type_map[i] = 2;

    }
    return atom_type_map;
}


void PairGRACEParallel::print_atomic_neighbours(const std::vector<int>& atom_type_map) const {
    std::stringstream ss;
    double **x = atom->x;
    tagint *tag = atom->tag;

    int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    ss << "\n[GRACE-DEBUG, Proc #" << comm->me << "]: Mapping, Positions, and Neighbors\n";
    ss << std::left << std::setw(6)  << "ii"   << " | "
            << std::setw(8)  << "Tag"   << " | "
            << std::setw(12) << "map(tag[i])" << " | " // <--- New Column Header
            << std::setw(8)  << "Shell" << " | "
            << std::setw(8)  << "ilist[ii]" << " | "
            << std::setw(25) << "Positions (x, y, z)" << " | "
            << "Neighbors (Local Indices)\n";
    ss << std::string(125, '-') << "\n";

    int nall = atom->nlocal + atom->nghost;

    for (int ii = 0; ii < nall; ++ii) {
        // Safety check: only access ilist within its bounds
        int i = ilist[ii];

        // If we are beyond inum, we should decide if we still want to print.
        // For this debug table, we'll use 'ii' as the index if we aren't using ilist mapping.
        int active_idx = (i != -1) ? i : ii;

        // 1. Label the shell (mapped to active_idx)
        std::string s_label;
        if (atom_type_map[active_idx] == 0)      s_label = "REAL";
        else if (atom_type_map[active_idx] == 1) s_label = "SHELL1";
        else if (atom_type_map[active_idx] == 2) s_label = "SHELL2";
        else                                     s_label = "OUT";

        // 2. Perform the map lookup
        int mapped_idx = atom->map(tag[active_idx]);

        // 3. Format row
        ss << std::left << std::setw(6)  << ii << " | "
                << std::setw(8)  << tag[active_idx] << " | "
                << std::setw(12) << mapped_idx << " | " // <--- The map result
                << std::setw(8)  << s_label << " | "
                << std::setw(8)  << (i != -1 ? std::to_string(i) : "-") << " | "
                << std::fixed << std::setprecision(3)
                << std::setw(7) << x[active_idx][0] << " "
                << std::setw(7) << x[active_idx][1] << " "
                << std::setw(7) << x[active_idx][2] << " | ";

        // 4. Neighbors column (based on center active_idx)
        if (firstneigh[active_idx] != nullptr) {
            int jnum = numneigh[active_idx];
            int* jlist = firstneigh[active_idx];
            for (int jj = 0; jj < jnum; ++jj) {
                int j = jlist[jj] & NEIGHMASK;
                ss << j << (jj == jnum - 1 ? "" : ",");
            }
        } else {
            ss << "no-list";
        }
        ss << "\n";
    }

    utils::logmesg(lmp, "{}\n", ss.str());
}

void print_remap(std::vector<int> index_remap, int me, LAMMPS *lmp) {
    std::stringstream ss_remap;
    ss_remap << "[GRACE-DEBUG-REMAP] Proc " << me << " Index Remap (LAMMPS->GRACE): ";
    for (size_t k = 0; k < index_remap.size(); ++k) {
        ss_remap << k << "->" << index_remap[k];
        if (k < index_remap.size() - 1) ss_remap << ", ";
    }
    ss_remap << "\n";
    utils::logmesg(lmp, "{}", ss_remap.str());
}

/* ---------------------------------------------------------------------- */
/**
signature_def['parallel_compute']:
  The given SavedModel SignatureDef contains the following input(s):
    inputs['atomic_mu_i'] tensor_info:  real  +  shell 1 + PAD1 (atom_padding)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_atomic_mu_i:0
    inputs['atomic_mu_i_1'] tensor_info: real  + PAD2 (atom_padding_1)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_atomic_mu_i_1:0
    inputs['batch_tot_nat_real'] tensor_info: =len(atomic_mu_i), placeholder
        dtype: DT_INT32
        shape: ()
        name: parallel_compute_batch_tot_nat_real:0
    inputs['bond_vector'] tensor_info: (r-r, r-s1, s1-r, s1-s1, s1-s2, NO s2-s1, s2-s2) + PAD3 (neighbor_padding)
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: parallel_compute_bond_vector:0
    inputs['ind_i'] tensor_info: (r-r, r-s1, s1-r, s1-s1, s1-s2,) + PAD3 (neighbor_padding)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_ind_i:0
    inputs['ind_i_1'] tensor_info: (r-r, r-s1) + PAD4 (neighbor_padding_1)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_ind_i_1:0
    inputs['ind_j'] tensor_info: (r-r, r-s1, s1-r, s1-s1, s1-s2,) + PAD3 (neighbor_padding)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_ind_j:0
    inputs['ind_j_1'] tensor_info: (r-r, r-s1) + PAD4 (neighbor_padding_1)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_ind_j_1:0
    inputs['mu_i'] tensor_info: (r-r, r-s1, s1-r, s1-s1, s1-s2,) + PAD3 (neighbor_padding)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_mu_i:0
    inputs['mu_j'] tensor_info: (r-r, r-s1, s1-r, s1-s1, s1-s2,) + PAD3 (neighbor_padding)
        dtype: DT_INT32
        shape: (-1)
        name: parallel_compute_mu_j:0
  The given SavedModel SignatureDef contains the following output(s):
    outputs['atomic_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall_1:0
    outputs['z_pair_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall_1:1
 */
void PairGRACEParallel::compute(int eflag, int vflag) {
    int i, j, ii, jj, inum, jnum;
    double delx, dely, delz, evdwl, my_cut_sq, xtmp, ytmp, ztmp;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;

    ev_init(eflag, vflag);

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int type_i;
    tagint *tag = atom->tag;

    // number of atoms in cell
    int nlocal = atom->nlocal;
    if (nlocal==0) return;

    int nall = nlocal + atom->nghost;

    int newton_pair = force->newton_pair;

    // inum: length of the neighborlists list
    inum = list->inum;

    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;

    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;

    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;
#ifdef GRACE_PRINT_DEBUG
    usleep(comm->me * 1000);
#endif
    data_timer.start();
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
    auto parallel_compute_inputs_sig = aceimpl->model->signatures.at(this->compute_function_name).inputs;
    // ------------------------------------------------------------------
    // 1. Identify Atom Shells
    // ------------------------------------------------------------------
    std::vector<int> atom_shell_map = get_atom_shell_mapping();

    // Count shells for tensor sizing
    int nshell1 = 0;
    int nshell2 = 0;
    // Note: iterating over nall to count correctly based on mapping logic
    for(int k=0; k<nall; k++) {
        if(atom_shell_map[k] == 1) nshell1++;
        else if(atom_shell_map[k] == 2) nshell2++;
    }

#ifdef GRACE_PRINT_DEBUG
    MPI_Barrier(world);
    usleep(comm->me * 1000);
    print_atomic_neighbours(atom_shell_map);
#endif

  // ------------------------------------------------------------------
    // 2. Prepare Atomic Species (atomic_mu_i) with WIDENED GAP
    // ------------------------------------------------------------------

    // Step A: Determine padding needed for the LOCAL subset first
    // This dictates how big the "Gap" must be.
    int tot_atoms_1 = aceimpl->atom_padding_1.update(nlocal);

    // The gap must be large enough to hold all padding for atomic_mu_i_1
    // gap_size = (padded_local_size) - (actual_real_atoms)
    // We enforce at least 1 fake atom if you want a guaranteed gap.
    int gap_size_atom = tot_atoms_1 - nlocal;
    // if (gap_size < 1) gap_size = 1;

    // Step B: Determine total size (Real + Gap + Shell1 + FinalPadding)
    // We treat the Gap as "real" occupied space now.
    int num_real_plus_gap_plus_ghosts = nlocal + gap_size_atom + nshell1;
    int tot_atoms = aceimpl->atom_padding.update(num_real_plus_gap_plus_ghosts);

    int fake_atom_type = 0;

    // We use the START of the gap as the canonical "fake atom index" for neighbor lists
    int fake_atom_ind_global = nlocal;

    aceimpl->atomic_mu_i_vector.assign(tot_atoms, fake_atom_type);

    // 1. Fill Real Atoms [0 ... nlocal-1]
    for (i = 0; i < nlocal; ++i) {
        aceimpl->atomic_mu_i_vector[i] = element_type_mapping[type[i]];
    }

    // 2. The Gap [nlocal ... nlocal + gap_size - 1]
    // Is already filled with fake_atom_type by assign().
    // We don't need to do anything here.

    // 3. Fill Ghost Atoms [nlocal + gap_size ... ]
    int ghost_start_idx = nlocal + gap_size_atom;
    int ghost_ptr = ghost_start_idx;

    // Remap vector: maps LAMMPS index -> GRACE index
    std::vector<int> index_remap(nall, -1);
    // for (int k = 0; k < nlocal; ++k) index_remap[k] = k; // Real maps 1:1

    for (int ii = nlocal; ii < nall; ++ii) {
         int k = ilist[ii];
         if (atom_shell_map[k] == 1 || atom_shell_map[k] == 2) {
             if (ghost_ptr < tot_atoms) {
                aceimpl->atomic_mu_i_vector[ghost_ptr] = element_type_mapping[type[k]];
             }
             // index_remap[k] = ghost_ptr;
             ghost_ptr++;
         }
         // Shell 2 usually maps to fake or stays -1 depending on needs
    }

    // Input 1: Main Atoms
    inputs.emplace_back( parallel_compute_inputs_sig.at("atomic_mu_i").name,
                       cppflow::tensor(aceimpl->atomic_mu_i_vector, {tot_atoms}));

    // Input 2: Local Atoms Subset
    // Now safe to copy prefix because the Gap is wide enough
    aceimpl->atomic_mu_i_1_vector.assign(aceimpl->atomic_mu_i_vector.begin(),
                                         aceimpl->atomic_mu_i_vector.begin() + tot_atoms_1);

    inputs.emplace_back(parallel_compute_inputs_sig.at( "atomic_mu_i_1").name,
                        cppflow::tensor(aceimpl->atomic_mu_i_1_vector, {tot_atoms_1}));


    // batch_tot_nat_real: number of extened atoms (w/o padding)
    inputs.emplace_back(parallel_compute_inputs_sig.at("batch_tot_nat_real").name,
                        cppflow::tensor(std::vector<int32_t>{nlocal}, {})); //nlocal

    // ind_i, ind_j: bonds
    // ind_i and ind_j - those pairs, where ind_i not in shell 2, but in shell 1 or real; MUST BE COMPACTIFIED!!! reindexed
    // ind_i_1 and ind_j_1 - slice of ind_i and ind_j, where ind_i is in real shell

   // ------------------------------------------------------------------
    // 4. Neighbor List Construction
    // Layout: [Real Bonds] [Gap (Fake)] [Ghost Bonds] [Final Padding]
    // ------------------------------------------------------------------

    // First Pass: Count Neighbors
    int n_real_bonds = 0;  // Bonds where i is Real
    int n_ghost_bonds = 0; // Bonds where i is Shell 1

    double cutoff_sq_default = cutoff * cutoff;
    double rsq;

    ghost_start_idx = nlocal + gap_size_atom;
    ghost_ptr = ghost_start_idx;

    for (ii = 0; ii < nall; ii++) {
        i = ilist[ii];
        int shell = atom_shell_map[i];
        if (shell != 0 && shell != 1) continue; // Skip Shell 2 or OUT

        type_i = type[i];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        int valid_neighs = 0;
        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj] & NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;

            my_cut_sq = cutoff_sq_default;
            if (is_custom_cutoffs) {
                my_cut_sq = cutoff_matrix_per_lammps_type[type_i][type[j]];
                my_cut_sq *= my_cut_sq;
            }

            if (rsq < my_cut_sq) {
                valid_neighs++;
                if (index_remap[i]==-1) {
                    if (shell==0) //identity map
                        index_remap[i] = i;
                    else  // shell1
                        index_remap[i] = ghost_ptr++;
                }
            }
        }

        if (shell == 0) n_real_bonds += valid_neighs;
        else n_ghost_bonds += valid_neighs;
    }

    // fill rest of the remap
    for (ii = nlocal; ii < nall; ii++) {
        i = ilist[ii];
        if (index_remap[i] == -1 )
            index_remap[i] = ghost_ptr++;
    }
#ifdef GRACE_PRINT_DEBUG
    print_remap(index_remap, comm->me, lmp);
#endif

    // Determine Padding Sizes
    int tot_neighbours_1 = aceimpl->neighbor_padding_1.update(n_real_bonds); // Real + Gap
    int gap_size_bonds = tot_neighbours_1 - n_real_bonds;
    int tot_neighbours = aceimpl->neighbor_padding.update(n_real_bonds + gap_size_bonds + n_ghost_bonds); // Total
#ifdef GRACE_PRINT_DEBUG
    // --- DEBUG PRINT START ---
    utils::logmesg(lmp, "[GRACE-DEBUG-PAD] Proc {}: n_real_bonds={}, tot_neighbours_1={} (gap_size_bonds={}), n_ghost_bonds={}, tot_neighbours={}\n",
                   comm->me, n_real_bonds, tot_neighbours_1, gap_size_bonds, n_ghost_bonds, tot_neighbours);
    // --- DEBUG PRINT END ---
#endif

    // Resize Vectors
    aceimpl->ind_i_vector.resize(tot_neighbours);
    aceimpl->ind_j_vector.resize(tot_neighbours);
    aceimpl->mu_i_vector.resize(tot_neighbours);
    aceimpl->mu_j_vector.resize(tot_neighbours);
    aceimpl->bond_vector.resize(3 * tot_neighbours);

  // Second Pass: Fill Vectors
    // We maintain two insertion pointers
    int ptr_real = 0;
    int ptr_ghost = tot_neighbours_1; // Starts after Gap

    // Pre-fill Gap with Fake Data
    for(int k=n_real_bonds; k < tot_neighbours_1; ++k) {
        aceimpl->ind_i_vector[k] = fake_atom_ind_global;
        aceimpl->ind_j_vector[k] = fake_atom_ind_global;
        aceimpl->mu_i_vector[k]  = fake_atom_type;
        aceimpl->mu_j_vector[k]  = fake_atom_type;
        aceimpl->bond_vector[3*k+0] = 1e6;
        aceimpl->bond_vector[3*k+1] = 1e6;
        aceimpl->bond_vector[3*k+2] = 1e6;
    }

    // Pre-fill Final Padding with Fake Data
    for(int k=tot_neighbours_1 + n_ghost_bonds; k < tot_neighbours; ++k) {
        aceimpl->ind_i_vector[k] = fake_atom_ind_global;
        aceimpl->ind_j_vector[k] = fake_atom_ind_global;
        aceimpl->mu_i_vector[k]  = fake_atom_type;
        aceimpl->mu_j_vector[k]  = fake_atom_type;
        aceimpl->bond_vector[3*k+0] = 1e6;
        aceimpl->bond_vector[3*k+1] = 1e6;
        aceimpl->bond_vector[3*k+2] = 1e6;
    }

    for (ii = 0; ii < nall; ii++) {
        i = ilist[ii];
        int shell = atom_shell_map[i];
        if (shell != 0 && shell != 1) continue;

        // Decide where to write based on shell
        int k;
        if (shell == 0) k = ptr_real;
        else k = ptr_ghost;

        type_i = type[i];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        // Remapped index for i
        // Real: i, Ghost: index_remap[i] (which is > nlocal)
        int i_remapped = index_remap[i];
        int mu_i_val = element_type_mapping[type_i];

        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj] & NEIGHMASK;
            delx = x[j][0]-xtmp;
            dely = x[j][1]-ytmp;
            delz = x[j][2]-ztmp;

            my_cut_sq = cutoff_sq_default;
            if (is_custom_cutoffs) {
                my_cut_sq = cutoff_matrix_per_lammps_type[type_i][type[j]];
                my_cut_sq *= my_cut_sq;
            }
            rsq = delx * delx + dely * dely + delz * delz;

            if (rsq < my_cut_sq) {
                // Remapped index for j
                // Note: j can be > nall in some LAMMPS configs? usually no.
                // We assume j is within range. If j is ghost, we map it.
                // If j is Shell 2 (not in index_remap), we default to fake or map conservatively?
                // Logic: j MUST be real or shell 1 for valid forces usually, but forces from shell2 neighbors affect Energy.
                // We need to map j consistently.
                int  j_remapped = index_remap[j];


                aceimpl->ind_i_vector[k] = i_remapped;
                aceimpl->ind_j_vector[k] = j_remapped;
                aceimpl->mu_i_vector[k] = mu_i_val;
                aceimpl->mu_j_vector[k] = element_type_mapping[type[j]];

                // Bond vector: j - i
                aceimpl->bond_vector[3*k+0] = delx;
                aceimpl->bond_vector[3*k+1] = dely;
                aceimpl->bond_vector[3*k+2] = delz;

                k++;
            }
        }

        // Update pointers
        if (shell == 0) ptr_real = k;
        else ptr_ghost = k;
    }

    inputs.emplace_back(parallel_compute_inputs_sig.at("bond_vector").name,
                        cppflow::tensor(aceimpl->bond_vector, {tot_neighbours, 3}));

    inputs.emplace_back(parallel_compute_inputs_sig.at("ind_i").name,
                        cppflow::tensor(aceimpl->ind_i_vector, {tot_neighbours}));
    inputs.emplace_back(parallel_compute_inputs_sig.at("ind_j").name,
                        cppflow::tensor(aceimpl->ind_j_vector, {tot_neighbours}));

    if (has_mu_i_op) {
        inputs.emplace_back(parallel_compute_inputs_sig.at("mu_i").name,
                            cppflow::tensor(aceimpl->mu_i_vector, {tot_neighbours}));
    }
    inputs.emplace_back(parallel_compute_inputs_sig.at("mu_j").name,
                        cppflow::tensor(aceimpl->mu_j_vector, {tot_neighbours}));

    // ------------------------------------------------------------------
    // 5. Prepare Subsets (_1 vectors) for Neighbors
    // ind_i_1 is just the prefix: Real Bonds + Gap
    // ------------------------------------------------------------------

    // We construct these by copying the prefix of the main vectors
    //TODO: maybe avoid creating vectors and just copy tot_neighbours_1 elements from ind_i_vector, ind_j_vector
    aceimpl->ind_i_vector_1.assign(aceimpl->ind_i_vector.begin(),
                                   aceimpl->ind_i_vector.begin() + tot_neighbours_1);
    aceimpl->ind_j_vector_1.assign(aceimpl->ind_j_vector.begin(),
                                   aceimpl->ind_j_vector.begin() + tot_neighbours_1);

    inputs.emplace_back(parallel_compute_inputs_sig.at("ind_i_1").name,
                        cppflow::tensor(aceimpl->ind_i_vector_1, {tot_neighbours_1}));
    inputs.emplace_back(parallel_compute_inputs_sig.at("ind_j_1").name,
                        cppflow::tensor(aceimpl->ind_j_vector_1, {tot_neighbours_1}));

#ifdef GRACE_PRINT_DEBUG
    print_tf_inputs(inputs, comm->me, lmp, true);
#endif

    data_timer.stop();
    tp_timer.start();
    auto parallel_compute_outputs_sig = aceimpl->model->signatures.at("parallel_compute").outputs;
    vector<string> output_names = {
        parallel_compute_outputs_sig.at("atomic_energy").name,// "StatefulPartitionedCall_1:0", // atomic_energy [nat,1]
        parallel_compute_outputs_sig.at("z_pair_f").name  // "StatefulPartitionedCall_1:1", // pair_f [n_bonds, 3]
    };


    //CALL MODEL
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            output_names
    );
    tp_timer.stop();
//    std::cout << "Ave.timing: " << (double) tp_timer.as_microseconds() / nlocal << " mcs/at" << std::endl;

    data_timer.start();
    auto &e_out = output[0]; // atomic_energy
    auto e_tens = e_out.get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));



    // pair forces
    auto &f_out = output[1];
    auto f_tens = f_out.get_tensor();
    const double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));

#ifdef GRACE_PRINT_DEBUG
    MPI_Barrier(world);
    usleep(comm->me*1000);
    // print_f_data(f_data,comm->me, lmp, aceimpl->ind_i_vector, aceimpl->ind_j_vector,  atom->tag);
    MPI_Barrier(world);
#endif

    // ------------------------------------------------------------------
    // 6. Force & Energy Tally (Pass 3)
    // Relies on deterministic order of ii/jj loops to match Tensor indices
    // ------------------------------------------------------------------


    // 1. Reset pointers EXACTLY as in the Fill Pass
    // ptr_real starts at 0
    // ptr_ghost starts AFTER the Real block + Gap
    ptr_real = 0;
    ptr_ghost =tot_neighbours_1; // Must match 'ghost_start_idx' from Fill Pass


    // loop only over LOCAL REAL atoms + shell 1 !!
    for (ii = 0; ii < nall; ++ii) {
        i = ilist[ii];
        int shell = atom_shell_map[i]; //

        // 2. Skip atoms we didn't send to TF (Shell 2/OUT)
        if (shell != 0 && shell != 1) continue;

        // 3. Select the correct pointer based on Shell type
        // This handles the "Jump over Gap" logic automatically
        int &curr_ptr = (shell == 0) ? ptr_real : ptr_ghost;

        type_i = type[i];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        // 4. Calculate Energy Index
        // Real atoms: i maps to i.
        // Ghost atoms: i maps to index_remap[i] (which handles the atomic gap)

        // Tally Energy (only for Real atoms usually)
        if (eflag_either && shell == 0) {
            evdwl = scale[type_i][type_i] * e_data[ii];
            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }

        for (jj = 0; jj < jnum; ++jj) {
            j = jlist[jj] & NEIGHMASK;

            // 5. Reproduce Filtering Logic
            // CRITICAL: You must use the EXACT same cutoff check as the Fill Pass.
            // If you process a neighbor here that you skipped in Fill (or vice versa),
            // the pointers will desync and you will read the wrong force.
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];

            my_cut_sq = cutoff_sq_default;
            if (is_custom_cutoffs) {
                my_cut_sq = cutoff_matrix_per_lammps_type[type_i][type[j]];
                my_cut_sq *= my_cut_sq;
            }
            rsq = delx * delx + dely * dely + delz * delz;

            if (rsq < my_cut_sq) {
                // 6. Read Force from the current pointer location
                int k = curr_ptr;

                double fx = f_data[3*k + 0];
                double fy = f_data[3*k + 1];
                double fz = f_data[3*k + 2];

                // 7. Apply Forces to LAMMPS array
                double s = scale[type_i][type_i];
                fij[0] = -s * fx;
                fij[1] = -s * fy;
                fij[2] = -s * fz;
#ifdef GRACE_PRINT_DEBUG
                // Print pair-specific force mapping
                utils::logmesg(lmp, "[GRACE-FORCE] Proc {}: Pair ({}-{}) | Tag ({}-{}) | BondIdx {} | F_tf: [{}, {}, {}]\n",
                                   comm->me, i, j, atom->tag[i], atom->tag[j], curr_ptr, fx, fy, fz);
#endif
                f[i][0] += fij[0];
                f[i][1] += fij[1];
                f[i][2] += fij[2];
                f[j][0] -= fij[0];
                f[j][1] -= fij[1];
                f[j][2] -= fij[2];

                if (vflag_either) {
                    ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], delx, dely, delz);
                    // (Optional: Add centroid stress logic here)

                    // Centroid Stress
                    if (cvflag_atom) {
                        const double fx = fij[0];
                        const double fy = fij[1];
                        const double fz = fij[2];

                        cvatom[i][0] += 0.5 * delx * fx; // xx
                        cvatom[i][1] += 0.5 * dely * fy; // yy
                        cvatom[i][2] += 0.5 * delz * fz; // zz
                        cvatom[i][3] += 0.5 * delx * fy; // xy
                        cvatom[i][4] += 0.5 * delx * fz; // xz
                        cvatom[i][5] += 0.5 * dely * fz; // yz
                        cvatom[i][6] += 0.5 * dely * fx; // yx
                        cvatom[i][7] += 0.5 * delz * fx; // zx
                        cvatom[i][8] += 0.5 * delz * fy; // zy


                        cvatom[j][0] += 0.5 * delx * fx; // xx
                        cvatom[j][1] += 0.5 * dely * fy; // yy
                        cvatom[j][2] += 0.5 * delz * fz; // zz
                        cvatom[j][3] += 0.5 * delx * fy; // xy
                        cvatom[j][4] += 0.5 * delx * fz; // xz
                        cvatom[j][5] += 0.5 * dely * fz; // yz
                        cvatom[j][6] += 0.5 * dely * fx; // yx
                        cvatom[j][7] += 0.5 * delz * fx; // zx
                        cvatom[j][8] += 0.5 * delz * fy; // zy
                    }
                }

                // 8. Increment the active pointer
                curr_ptr++; //it is reference
            }
        } // loop over neighbours
    } // loop over atoms -i

    if (vflag_fdotr) virial_fdotr_compute();



    data_timer.stop();
    // end modifications YL
}

#endif //#ifndef NO_GRACE_TF
