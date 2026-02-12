#ifndef NO_GRACE_TF
// #define GRACE_DEBUG
#define GRACE_PROFILE

#include "pair_grace_2layer_parallel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"

#include <cstring>
#include <algorithm>
#include <tuple>
#include "yaml-cpp/yaml.h"

#include "utils_pace.h"
#include "utils_grace.h"

// CppFlow headers
#include <unistd.h>
#include <cppflow/ops.h>
#include <cppflow/model.h>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>

namespace LAMMPS_NS {
    struct GRACE2LayerImpl {
        GRACE2LayerImpl() : model(nullptr) {}
        ~GRACE2LayerImpl() { delete model; }
        cppflow::model *model;

        GRACE::GracePaddingDimension all_atoms_padding;
        GRACE::GracePaddingDimension real_atoms_padding;
        GRACE::GracePaddingDimension neighbor_padding;

        std::vector<int32_t> mu_i;
        std::vector<int32_t> mu_j;
        std::vector<int32_t> ind_i;
        std::vector<int32_t> ind_j;
        std::vector<double> bond_vector;
        std::vector<int32_t> atomic_mu_i; // for real + ghost atoms + pad
        std::vector<int32_t> atomic_mu_i_1; // for real atoms + pad

        std::vector<double> grad_bond_vector; // Temporary storage for summing forces

        int n_all_atoms_padded = 0;
        int n_local_atoms_padded = 0;
        int n_neighbours_padded = 0;
        int nlocal_bonds = 0;
    };
}

using namespace LAMMPS_NS;

PairGRACE2LayerParallel::PairGRACE2LayerParallel(LAMMPS *lmp) : Pair(lmp) {
    single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    aceimpl = new GRACE2LayerImpl;
    scale = nullptr;
    
#ifdef GRACE_PROFILE
    total_timer.init();
    data_timer.init();
    tp_timer.init();
    comm_timer.init();
    model1_timer.init();
    model2_timer.init();
    model3_timer.init();
#endif

    no_virial_fdotr_compute = 1;
    chunksize = 4096;
    nelements = 0;
}

PairGRACE2LayerParallel::~PairGRACE2LayerParallel() {
    if (copymode) return;
    delete aceimpl;
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }
}

void PairGRACE2LayerParallel::allocate() {
    allocated = 1;
    int n = atom->ntypes + 1;
    memory->create(setflag, n, n, "pair:setflag");
    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(scale, n, n, "pair:scale");
    map = new int[n];
}

void PairGRACE2LayerParallel::settings(int narg, char **arg) {
    if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace", error);

    // ACE potentials are parameterized in metal units
    if (strcmp("metal", update->unit_style) != 0)
        error->all(FLERR, "GRACE potentials require 'metal' units");

    auto tf_version = TF_Version();
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", tf_version);

    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "padding") == 0) {
            neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;

        } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
            pad_verbose = true;
            iarg += 1;
        } else if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
            max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Maximum number of recompilation during padding reduction: {}\n",
                               max_number_of_reduction);
        } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
            reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
            iarg += 2;
            if (comm->me == 0)
                utils::logmesg(lmp, "[GRACE] Reducing padding fraction: {}\n", reducing_neigh_padding_fraction);
        } else
            error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
    }

    do_padding = (neigh_padding_fraction > 0);
    if(do_padding) {
        if (comm->me == 0)
            utils::logmesg(lmp, "[GRACE] Neighbour padding is ON, padding fraction: {}, max padding fraction before reduction: {}, max number of reduction(s): {}\n",
                           neigh_padding_fraction, reducing_neigh_padding_fraction, max_number_of_reduction);
    }

    // Apply to padding helpers
    aceimpl->all_atoms_padding.enabled = do_padding;
    aceimpl->all_atoms_padding.padding_fraction = neigh_padding_fraction;
    aceimpl->all_atoms_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
    aceimpl->all_atoms_padding.max_reductions = max_number_of_reduction;
    aceimpl->all_atoms_padding.verbose = pad_verbose;

    aceimpl->real_atoms_padding.enabled = do_padding;
    aceimpl->real_atoms_padding.padding_fraction = neigh_padding_fraction;
    aceimpl->real_atoms_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
    aceimpl->real_atoms_padding.max_reductions = max_number_of_reduction;
    aceimpl->real_atoms_padding.verbose = pad_verbose;

    aceimpl->neighbor_padding.enabled = do_padding;
    aceimpl->neighbor_padding.padding_fraction = neigh_padding_fraction;
    aceimpl->neighbor_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
    aceimpl->neighbor_padding.max_reductions = max_number_of_reduction;
    aceimpl->neighbor_padding.verbose = pad_verbose;

    if (!pair_forces && comm->nprocs > 1) {
        pair_forces = true;
        if (comm->me == 0)
            utils::logmesg(lmp,
                           "[GRACE] ENFORCE pair-force mode to ON, because number of processes {} is more than one.\n",
                           comm->nprocs);
    }

    
    // mark as available centroid stress flag
    centroidstressflag = CENTROID_AVAIL;
    
}

void PairGRACE2LayerParallel::coeff(int narg, char **arg) {
    if (!allocated) allocate();

    map_element2type(narg - 3, arg + 3);
    auto potential_path = std::string(arg[2]);

    delete aceimpl->model;
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);

    const std::vector<uint8_t> config_bytes = { 0x32, 0x05, 0x82, 0x01, 0x02, 0x18, 0x00 };
    aceimpl->model = new cppflow::model(potential_path, config_bytes);

    YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
    elements_name = metadata_yaml["chemical_symbols"].as<std::vector<std::string>>();
    nelements = static_cast<int>(elements_name.size());
    for (int mu = 0; mu < nelements; mu++) {
        elements_to_index_map[elements_name.at(mu)] = mu;
    }
    cutoff = metadata_yaml["cutoff"].as<double>();

    // Parse parallel communication features
    if (metadata_yaml["parallel_communication"]) {
        YAML_PACE::Node parallel_comm = metadata_yaml["parallel_communication"];
        for (YAML_PACE::const_iterator it = parallel_comm.begin(); it != parallel_comm.end(); ++it) {
            std::string key = it->first.as<std::string>();
            YAML_PACE::Node node = it->second;
            
            std::vector<int64_t> shape;
            int size = 1;
            
            if (node["shape"]) {
                for (const auto& dim : node["shape"]) {
                    int d = dim.as<int>();
                    shape.push_back(d);
                    size *= d;
                }
            } else {
                 // Fallback if shape is directly the value (old format support if needed, though user specified new format)
                 // Assuming new format as per request
                 error->all(FLERR, "[GRACE] Feature '{}' missing 'shape' in metadata.yaml", key);
            }
            
            bool non_local; // Default to true
            if (node["non_local"]) {
                non_local = node["non_local"].as<bool>();
            } else {
                error->all(FLERR, "[GRACE] Feature '{}' missing 'non_local' in metadata.yaml", key);
            }

            feature_shapes[key] = shape;
            feature_sizes[key] = size;
            feature_is_local[key] = !non_local;

            if (comm->me == 0) {
                std::string shape_str;
                for (size_t i = 0; i < shape.size(); ++i) {
                    shape_str += std::to_string(shape[i]);
                    if (i < shape.size() - 1) shape_str += ", ";
                }
                utils::logmesg(lmp, "[GRACE] Feature '{}' shape: [{}], size: {}, non_local: {}\n", key, shape_str, size, non_local);
            }
        }
    } else {
        error->all(FLERR, "[GRACE] metadata.yaml missing 'parallel_communication' section");
    }

    // Calculate total communication size
    comm_forward = 0;
    comm_reverse = 0;
    for (const auto& [key, size] : feature_sizes) {
        comm_forward += size;
        if (!feature_is_local[key]) {
            comm_reverse += size;
        }
    }

    const int ntypes = atom->ntypes;
    element_type_mapping.resize(ntypes + 1);
    for (int i = 1; i <= ntypes; i++) {
        char *elemname = arg[2 + i];
        if (strcmp(elemname, "NULL") == 0) {
            element_type_mapping[i] = -1;
            map[i] = -1;
        } else {
            int mu = elements_to_index_map.at(elemname);
            map[i] = mu;
            element_type_mapping[i] = mu;
        }
    }

    for (int i = 1; i <= ntypes; i++)
        for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;

    if (!aceimpl->model->has_signature(forward_layer_1_name)) error->all(FLERR, "Model missing forward_layer_1");
    if (!aceimpl->model->has_signature(backward_layer_2_name)) error->all(FLERR, "Model missing backward_layer_2");
    if (!aceimpl->model->has_signature(backward_layer_1_name)) error->all(FLERR, "Model missing backward_layer_1");

    if (metadata_yaml["cutoff_matrix"]) {
        cutoff_matrix = metadata_yaml["cutoff_matrix"].as<vector<vector<double>>>();
        is_custom_cutoffs = true;
        cutoff_matrix_per_lammps_type.resize(ntypes + 1, vector<double>(ntypes + 1));
        for (int i = 1; i <= ntypes; i++) {
            for (int j = 1; j <= ntypes; j++) {
                cutoff_matrix_per_lammps_type[i][j] = cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
            }
        }
    }

    if (comm->me == 0) {
#ifdef GRACE_DEBUG
        utils::logmesg(lmp, "[GRACE-DEBUG] Available signatures:\n");
        for (auto const& [name, sig] : aceimpl->model->signatures) {
            utils::logmesg(lmp, "  Signature: {}\n", name);
            utils::logmesg(lmp, "    Inputs:\n");
            for (auto const& [in_key, in_val] : sig.inputs) utils::logmesg(lmp, "      {} -> {}\n", in_key, in_val.name);
            utils::logmesg(lmp, "    Outputs:\n");
            for (auto const& [out_key, out_val] : sig.outputs) utils::logmesg(lmp, "      {} -> {}\n", out_key, out_val.name);
        }
#endif
    }
}

void PairGRACE2LayerParallel::init_style() {
    if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
    if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");
    neighbor->add_request(this, NeighConst::REQ_FULL);
    
    // Initialize atom map for ghost-to-parent lookup in force tallying
    if (atom->map_style == Atom::MAP_NONE) {
        atom->map_init();
        atom->map_set();
    }
}


double PairGRACE2LayerParallel::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    return cutoff;
}

void PairGRACE2LayerParallel::compute(int eflag, int vflag) {

#ifdef GRACE_PROFILE
    total_timer.init();
    data_timer.init();
    tp_timer.init();
    comm_timer.init();
    model1_timer.init();
    model2_timer.init();
    model3_timer.init();

    total_timer.start();
    data_timer.start();
#endif
    ev_init(eflag, vflag);

    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;

    aceimpl->n_all_atoms_padded = aceimpl->all_atoms_padding.update(nall);
    aceimpl->n_local_atoms_padded =aceimpl->real_atoms_padding.update(nlocal);
    if (comm->me == 0 && pad_verbose && aceimpl->all_atoms_padding.last_update_triggered_resize()) {
        utils::logmesg(lmp, "[GRACE] Atom padding resize: all atoms {} -> {}, real atoms {} -> {}\n",
            nall, aceimpl->n_all_atoms_padded,
            nlocal, aceimpl->n_local_atoms_padded
            );
    }
    
    int n_bonds = 0;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;
    double **x = atom->x;
    int *type = atom->type;
    double cutoff_sq = cutoff * cutoff;

    for (int ii = 0; ii < nlocal; ii++) {
        int i = ilist[ii];
        int type_i = type[i];
        double xtmp = x[i][0];
        double ytmp = x[i][1];
        double ztmp = x[i][2];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
            int j = jlist[jj] & NEIGHMASK;
            double dx = xtmp - x[j][0];
            double dy = ytmp - x[j][1];
            double dz = ztmp - x[j][2];
            if (is_custom_cutoffs) {
                int type_j = type[j];
                double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
                cutoff_sq = cur_cutoff * cur_cutoff;
            }
            if (dx*dx + dy*dy + dz*dz < cutoff_sq) n_bonds++;
        }
    }
    aceimpl->nlocal_bonds = n_bonds;
    aceimpl->n_neighbours_padded = aceimpl->neighbor_padding.update(n_bonds);
    if (comm->me == 0 && pad_verbose && aceimpl->neighbor_padding.last_update_triggered_resize()) {
        utils::logmesg(lmp, "[GRACE] Neighbor padding resize: {} -> {}\n", n_bonds, aceimpl->n_neighbours_padded);
    }

    // Resize features and gradients based on parsed sizes and locality
    for (const auto& [key, size] : feature_sizes) {
        int n_padded = feature_is_local[key] ? aceimpl->n_local_atoms_padded : aceimpl->n_all_atoms_padded;
        features[key].resize(n_padded * size);
        gradients[key].assign(n_padded * size, 0.0);
    }

    aceimpl->mu_i.resize(aceimpl->n_neighbours_padded);
    aceimpl->mu_j.resize(aceimpl->n_neighbours_padded);
    aceimpl->ind_i.resize(aceimpl->n_neighbours_padded);
    aceimpl->ind_j.resize(aceimpl->n_neighbours_padded);
    aceimpl->bond_vector.resize(3 * aceimpl->n_neighbours_padded, 1e6);
    grad_bv_L2.resize(3 * aceimpl->n_neighbours_padded);

    aceimpl->atomic_mu_i.assign(aceimpl->n_all_atoms_padded, element_type_mapping[type[0]]);
    for (int i = 0; i < nall; i++) aceimpl->atomic_mu_i[i] = element_type_mapping[type[i]];

    aceimpl->atomic_mu_i_1.assign(aceimpl->n_local_atoms_padded, element_type_mapping[type[0]]);
    for (int i = 0; i < nlocal; i++) aceimpl->atomic_mu_i_1[i] = element_type_mapping[type[i]];
    
    int tot_ind = 0;
    for (int ii = 0; ii < nlocal; ii++) {
        int i = ilist[ii];
        int type_i = type[i];
        double xtmp = x[i][0]; double ytmp = x[i][1]; double ztmp = x[i][2];
        int *jlist = firstneigh[i]; int jnum = numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
            int j = jlist[jj] & NEIGHMASK;
            double dx = x[j][0] - xtmp;
            double dy = x[j][1] - ytmp;
            double dz = x[j][2] - ztmp;
            if (is_custom_cutoffs) {
                int type_j = type[j];
                double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
                cutoff_sq = cur_cutoff * cur_cutoff;
            }
            if (dx*dx + dy*dy + dz*dz < cutoff_sq) {
                aceimpl->ind_i[tot_ind] = i;
                aceimpl->ind_j[tot_ind] = j;
                aceimpl->mu_i[tot_ind] = element_type_mapping[type[i]];
                aceimpl->mu_j[tot_ind] = element_type_mapping[type[j]];
                aceimpl->bond_vector[3*tot_ind + 0] = dx;
                aceimpl->bond_vector[3*tot_ind + 1] = dy;
                aceimpl->bond_vector[3*tot_ind + 2] = dz;
                tot_ind++;
            }
        }
    }
    // int fake_at = aceimpl->n_all_atoms_padded - 1;
    int fake_at = aceimpl->n_local_atoms_padded - 1;
    for (int k = n_bonds; k < aceimpl->n_neighbours_padded; k++) {
        aceimpl->ind_i[k] = aceimpl->ind_j[k] = fake_at;
        aceimpl->mu_i[k] = aceimpl->mu_j[k] = 0;
    }

#ifdef GRACE_PROFILE
    data_timer.stop();
    model1_timer.start();
#endif
    run_forward_layer_1();
#ifdef GRACE_PROFILE
    model1_timer.stop();

    comm_timer.start();
#endif
    comm->forward_comm(this);
#ifdef GRACE_PROFILE
    comm_timer.stop();

    model2_timer.start();
#endif
    run_backward_layer_2(eflag, vflag);
#ifdef GRACE_PROFILE
    model2_timer.stop();

    comm_timer.start();
#endif
    comm->reverse_comm(this);
#ifdef GRACE_PROFILE
    comm_timer.stop();

    data_timer.start();
#endif
    
    // Zero ghost/padding gradients after reverse_comm has accumulated them to parents
    // This ensures ghosts don't contribute again in backward_layer_1
    for (auto& [key, grad_vec] : gradients) {
        int size = feature_sizes[key];
        int n_padded = feature_is_local[key] ? aceimpl->n_local_atoms_padded : aceimpl->n_all_atoms_padded;
        if (n_padded > atom->nlocal) {
            std::fill_n(&grad_vec[atom->nlocal * size], (n_padded - atom->nlocal) * size, 0.0);
        }
    }
    
#ifdef GRACE_PROFILE
    data_timer.stop();
    model3_timer.start();
#endif
    run_backward_layer_1();
#ifdef GRACE_PROFILE
    model3_timer.stop();
    data_timer.start();
#endif

    
#ifdef GRACE_DEBUG
    if (comm->me == 0) {
        int n_bonds_print = std::min(aceimpl->nlocal_bonds, 8);
        utils::logmesg(lmp, "[GRACE-DEBUG] Sum (grad L2 + grad L1):\n");
        for (int k = 0; k < n_bonds_print; k++) {
            double fx_tot = aceimpl->grad_bond_vector[3 * k + 0];
            double fy_tot = aceimpl->grad_bond_vector[3 * k + 1];
            double fz_tot = aceimpl->grad_bond_vector[3 * k + 2];
            utils::logmesg(lmp, "{}[{:12.8f} {:12.8f} {:12.8f}]{}\n", 
                          (k == 0 ? "[" : " "), fx_tot, fy_tot, fz_tot, 
                          (k == n_bonds_print - 1 ? "]" : ""));
        }
        
        // Tally forces for debug printing (using atom->map for ghost-to-parent lookup)
        std::vector<std::vector<double>> f_par(atom->nlocal, std::vector<double>(3, 0.0));
        for (int k = 0; k < aceimpl->nlocal_bonds; k++) {
            int i = aceimpl->ind_i[k];
            int j = aceimpl->ind_j[k];
            
            // Map ghost j to parent local atom using LAMMPS atom map
            if (j >= nlocal) {
                j = atom->map(atom->tag[j]);
            }
            
            double sc = scale[type[i]][type[i]];
            double fx = sc * aceimpl->grad_bond_vector[3 * k + 0];
            double fy = sc * aceimpl->grad_bond_vector[3 * k + 1];
            double fz = sc * aceimpl->grad_bond_vector[3 * k + 2];
            if (i < atom->nlocal) {
                f_par[i][0] += fx; f_par[i][1] += fy; f_par[i][2] += fz;
            }
            if (j >= 0 && j < atom->nlocal) {
                f_par[j][0] -= fx; f_par[j][1] -= fy; f_par[j][2] -= fz;
            }
        }

        int n_atoms_print = std::min(atom->nlocal, 8);
        utils::logmesg(lmp, "[GRACE-DEBUG] Computed Atomic Forces from Parallel Grads:\n");
        for (int i = 0; i < n_atoms_print; i++) {
            utils::logmesg(lmp, "{}[{:12.8f} {:12.8f} {:12.8f}]{}\n", 
                          (i == 0 ? "[" : " "), f_par[i][0], f_par[i][1], f_par[i][2],
                          (i == n_atoms_print - 1 ? "]" : ""));
        }
    }
#endif

    // Final force tally
    // For ghosts:
    //   - Same-proc PBC ghosts: atom->map(tag) returns local index, use it
    //   - Cross-proc ghosts: atom->map(tag) returns -1, apply to ghost f[j]
    //     and LAMMPS reverse_comm will communicate forces back to owner
    double **f = atom->f;
    int newton_pair = force->newton_pair;
    
    for (int k = 0; k < aceimpl->nlocal_bonds; k++) {
        int i = aceimpl->ind_i[k];
        int j = aceimpl->ind_j[k];
        int j_orig = j;  // Keep original for virial
        
        // Try to map ghost j to its parent local atom
        if (j >= nlocal) {
            int j_local = atom->map(atom->tag[j]);
            if (j_local >= 0 && j_local < nlocal) {
                // Same-proc PBC ghost: use the local parent atom
                j = j_local;
            }
            // Else: cross-proc ghost, keep j as is and apply to ghost atom
            // LAMMPS reverse_comm will communicate the force back to owner
        }
        
        double sc = scale[type[i]][type[i]];
        double fx = sc * aceimpl->grad_bond_vector[3 * k + 0];
        double fy = sc * aceimpl->grad_bond_vector[3 * k + 1];
        double fz = sc * aceimpl->grad_bond_vector[3 * k + 2];
        f[i][0] += fx; f[i][1] += fy; f[i][2] += fz;
        f[j][0] -= fx; f[j][1] -= fy; f[j][2] -= fz;
        if (evflag) {
            double dx = -aceimpl->bond_vector[3*k+0];
            double dy = -aceimpl->bond_vector[3*k+1];
            double dz = -aceimpl->bond_vector[3*k+2];
            ev_tally_xyz(i, j_orig, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, dx, dy, dz);

            if (cvflag_atom) {
                cvatom[i][0] += 0.5 * dx * fx; // xx
                cvatom[i][1] += 0.5 * dy * fy; // yy
                cvatom[i][2] += 0.5 * dz * fz; // zz
                cvatom[i][3] += 0.5 * dx * fy; // xy
                cvatom[i][4] += 0.5 * dx * fz; // xz
                cvatom[i][5] += 0.5 * dy * fz; // yz
                cvatom[i][6] += 0.5 * dy * fx; // yx
                cvatom[i][7] += 0.5 * dz * fx; // zx
                cvatom[i][8] += 0.5 * dz * fy; // zy

                cvatom[j][0] += 0.5 * dx * fx; // xx
                cvatom[j][1] += 0.5 * dy * fy; // yy
                cvatom[j][2] += 0.5 * dz * fz; // zz
                cvatom[j][3] += 0.5 * dx * fy; // xy
                cvatom[j][4] += 0.5 * dx * fz; // xz
                cvatom[j][5] += 0.5 * dy * fz; // yz
                cvatom[j][6] += 0.5 * dy * fx; // yx
                cvatom[j][7] += 0.5 * dz * fx; // zx
                cvatom[j][8] += 0.5 * dz * fy; // zy
            }
        }
    }



    if (vflag_fdotr) virial_fdotr_compute();
#ifdef GRACE_PROFILE
    data_timer.stop();
    total_timer.stop();
#endif

#ifdef GRACE_PROFILE
    // if (comm->me == 0) { // PRINT FOR ALL RANKS
    {
        double d_t = data_timer.as_microseconds();
        double c_t = comm_timer.as_microseconds();
        double m1_t = model1_timer.as_microseconds();
        double m2_t = model2_timer.as_microseconds();
        double m3_t = model3_timer.as_microseconds();
        double total_t = total_timer.as_microseconds();

        auto pct = [&](double t) { return (total_t > 0) ? (t / total_t * 100.0) : 0.0; };

        // Use fprintf to stderr to ensure immediate output from all ranks
        fprintf(stderr, "[GRACE-PROFILE] [Rank %d] Timings (mcs): Data: %.1f (%.1f%%) | Comm: %.1f (%.1f%%) | M1: %.1f (%.1f%%) | M2: %.1f (%.1f%%) | M3: %.1f (%.1f%%) | Total: %.1f\n",
                comm->me, d_t, pct(d_t), c_t, pct(c_t), m1_t, pct(m1_t), m2_t, pct(m2_t), m3_t, pct(m3_t), total_t);
        
        fprintf(stderr, "[GRACE-PROFILE] [Rank %d] Array sizes: Atoms: %d (padded: %d) | Neighbors: %d (padded: %d)\n",
                comm->me, nall, aceimpl->n_all_atoms_padded, aceimpl->nlocal_bonds, aceimpl->n_neighbours_padded);
    }
#endif
}

int PairGRACE2LayerParallel::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {
    int m = 0;
    
    // Hoist map lookups out of the atom loop for performance.
    // We pre-fetch base pointers and sizes so lookups happen once per routine, 
    // rather than n * num_features times.
    struct FeaturePtr {
        double* data;
        int size;
    };
    std::vector<FeaturePtr> active_features;
    for (const auto& [key, size] : feature_sizes) {
        if (!feature_is_local[key]) {
            active_features.push_back({features[key].data(), size});
        }
    }

    for (int i = 0; i < n; i++) {
        int j = list[i];
        for (const auto& feat : active_features) {
            std::copy_n(&feat.data[j * feat.size], feat.size, &buf[m]);
            m += feat.size;
        }
    }
    return m;
}

void PairGRACE2LayerParallel::unpack_forward_comm(int n, int first, double *buf) {
    int m = 0; int last = first + n;
    
    // Hoist map lookups out of the atom loop for performance.
    struct FeaturePtr {
        double* data;
        int size;
    };
    std::vector<FeaturePtr> active_features;
    for (const auto& [key, size] : feature_sizes) {
        if (!feature_is_local[key]) {
            active_features.push_back({features[key].data(), size});
        }
    }

    for (int i = first; i < last; i++) {
        for (const auto& feat : active_features) {
            std::copy_n(&buf[m], feat.size, &feat.data[i * feat.size]);
            m += feat.size;
        }
    }
}

int PairGRACE2LayerParallel::pack_reverse_comm(int n, int first, double *buf) {
    int m = 0; int last = first + n;

    // Hoist map lookups out of the atom loop for performance.
    struct FeaturePtr {
        double* data;
        int size;
    };
    std::vector<FeaturePtr> active_features;
    for (const auto& [key, size] : feature_sizes) {
        if (!feature_is_local[key]) {
            active_features.push_back({gradients[key].data(), size});
        }
    }

    for (int i = first; i < last; i++) {
        for (const auto& feat : active_features) {
            std::copy_n(&feat.data[i * feat.size], feat.size, &buf[m]);
            m += feat.size;
        }
    }
    return m;
}

void PairGRACE2LayerParallel::unpack_reverse_comm(int n, int *list, double *buf) {
    int m = 0;

    // Hoist map lookups out of the atom loop for performance.
    struct FeaturePtr {
        double* data;
        int size;
    };
    std::vector<FeaturePtr> active_features;
    for (const auto& [key, size] : feature_sizes) {
        if (!feature_is_local[key]) {
            active_features.push_back({gradients[key].data(), size});
        }
    }

    for (int i = 0; i < n; i++) {
        int j = list[i];
        for (const auto& feat : active_features) {
            for (int k = 0; k < feat.size; k++) feat.data[j * feat.size + k] += buf[m++];
        }
    }
}


void PairGRACE2LayerParallel::run_forward_layer_1() {
#ifdef GRACE_DEBUG
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE-DEBUG] Entering run_forward_layer_1\n");
#endif
    auto sig = aceimpl->model->signatures.at(forward_layer_1_name);
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
    
    auto add_input = [&](const std::string& key, const cppflow::tensor& t) {
        if (sig.inputs.count(key)) {
            inputs.emplace_back(sig.inputs.at(key).name, t);
        } else {
            error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature", key);
        }
    };

    add_input("atomic_mu_i", cppflow::tensor(aceimpl->atomic_mu_i, {aceimpl->n_all_atoms_padded}));
    add_input("atomic_mu_i_1", cppflow::tensor(aceimpl->atomic_mu_i_1, {aceimpl->n_local_atoms_padded}));
    add_input("ind_i", cppflow::tensor(aceimpl->ind_i, {aceimpl->n_neighbours_padded}));
    add_input("ind_j", cppflow::tensor(aceimpl->ind_j, {aceimpl->n_neighbours_padded}));
    add_input("mu_i", cppflow::tensor(aceimpl->mu_i, {aceimpl->n_neighbours_padded}));
    add_input("mu_j", cppflow::tensor(aceimpl->mu_j, {aceimpl->n_neighbours_padded}));
    add_input("bond_vector", cppflow::tensor(aceimpl->bond_vector, {aceimpl->n_neighbours_padded, 3}));
    add_input("batch_tot_nat_real", cppflow::tensor(std::vector<int32_t>{atom->nlocal}, {}));

#ifdef GRACE_DEBUG
    if (comm->me == 0) print_tensors(forward_layer_1_name, inputs);
#endif
    
    std::vector<std::string> out_names;
    std::vector<std::string> ordered_keys;
    
    // We need to request outputs for all keys in feature_shapes
    for (const auto& [key, shape] : feature_shapes) {
        if (sig.outputs.count(key)) {
            out_names.push_back(sig.outputs.at(key).name);
            ordered_keys.push_back(key);
        } else {
            error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature '{}' outputs", key, forward_layer_1_name);
        }
    }

    auto outputs = aceimpl->model->operator()(inputs, out_names);
    
#ifdef GRACE_DEBUG
    if (comm->me == 0) {
        std::vector<std::tuple<std::string, cppflow::tensor>> out_tensors;
        for (size_t i=0; i<outputs.size(); ++i) out_tensors.push_back({out_names[i], outputs[i]});
        print_tensors(forward_layer_1_name, out_tensors, "Output");
    }
#endif
    
    for (size_t i = 0; i < outputs.size(); ++i) {
        std::string key = ordered_keys[i];
        auto tens = outputs[i].get_tensor();
        const double *data = static_cast<double *>(TF_TensorData(tens.get()));
        int size = feature_sizes[key];
        // Copy exactly what we got (n_local_atoms_padded * size)
        std::copy_n(data, aceimpl->n_local_atoms_padded * size, &features[key][0]);
    }
}

void PairGRACE2LayerParallel::run_backward_layer_2(int eflag, int vflag) {
#ifdef GRACE_DEBUG
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE-DEBUG] Entering run_backward_layer_2\n");
#endif
    auto sig = aceimpl->model->signatures.at(backward_layer_2_name);
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    auto add_input = [&](const std::string& key, const cppflow::tensor& t) {
        if (sig.inputs.count(key)) {
            inputs.emplace_back(sig.inputs.at(key).name, t);
        } else {
            error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature", key);
        }
    };

    add_input("atomic_mu_i", cppflow::tensor(aceimpl->atomic_mu_i, {aceimpl->n_all_atoms_padded}));
    add_input("atomic_mu_i_1", cppflow::tensor(aceimpl->atomic_mu_i_1, {aceimpl->n_local_atoms_padded}));
    add_input("ind_i", cppflow::tensor(aceimpl->ind_i, {aceimpl->n_neighbours_padded}));
    add_input("ind_j", cppflow::tensor(aceimpl->ind_j, {aceimpl->n_neighbours_padded}));
    add_input("mu_i", cppflow::tensor(aceimpl->mu_i, {aceimpl->n_neighbours_padded}));
    add_input("mu_j", cppflow::tensor(aceimpl->mu_j, {aceimpl->n_neighbours_padded}));
    add_input("bond_vector", cppflow::tensor(aceimpl->bond_vector, {aceimpl->n_neighbours_padded, 3}));
    add_input("batch_tot_nat_real", cppflow::tensor(std::vector<int32_t>{atom->nlocal}, {}));
    
    for (const auto& [key, shape] : feature_shapes) {
        if (feature_is_local[key]) {
             std::vector<int64_t> local_shape = {aceimpl->n_local_atoms_padded};
             local_shape.insert(local_shape.end(), shape.begin(), shape.end());
             add_input(key, cppflow::tensor(features[key], local_shape)); 
        } else {
             std::vector<int64_t> full_shape = {aceimpl->n_all_atoms_padded};
             full_shape.insert(full_shape.end(), shape.begin(), shape.end());
             add_input(key, cppflow::tensor(features[key], full_shape));
        }
    }

#ifdef GRACE_DEBUG
    if (comm->me == 0) print_tensors(backward_layer_2_name, inputs);
#endif
    
    std::vector<std::string> out_names;
    std::vector<std::string> ordered_keys;
    
    // Energy
    if (sig.outputs.count(ENERGY_KEY)) {
        out_names.push_back(sig.outputs.at(ENERGY_KEY).name);
        ordered_keys.push_back(ENERGY_KEY);
    } else error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature '{}' outputs", ENERGY_KEY, backward_layer_2_name);

    // Gradients for features
    for (const auto& [key, shape] : feature_shapes) {
        std::string grad_key = "grad_" + key;
        if (sig.outputs.count(grad_key)) {
            out_names.push_back(sig.outputs.at(grad_key).name);
            ordered_keys.push_back(grad_key);
        } else {
             error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature '{}' outputs", grad_key, backward_layer_2_name);
        }
    }

    // Grad bond vector
    if (sig.outputs.count(GRAD_BOND_KEY)) {
        out_names.push_back(sig.outputs.at(GRAD_BOND_KEY).name);
        ordered_keys.push_back(GRAD_BOND_KEY);
    } else error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature '{}' outputs", GRAD_BOND_KEY, backward_layer_2_name);
    
    auto outputs = aceimpl->model->operator()(inputs, out_names);

#ifdef GRACE_DEBUG
    if (comm->me == 0) {
        std::vector<std::tuple<std::string, cppflow::tensor>> out_tensors;
        for (size_t i=0; i<outputs.size(); ++i) out_tensors.push_back({out_names[i], outputs[i]});
        print_tensors(backward_layer_2_name, out_tensors, "Output");
    }
#endif

    int out_idx = 0;
    
    // Energy
    auto e_tens = outputs[out_idx++].get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));
    if (eflag_either) {
        for (int i = 0; i < atom->nlocal; i++)
            ev_tally_full(i, 2.0 * e_data[i], 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    // Gradients
    for (const auto& [key, shape] : feature_shapes) {
        std::string grad_key = "grad_" + key;
        // We iterate in the same order as we pushed to ordered_keys
        // ordered_keys[out_idx] should be grad_key
        
        auto g_tens = outputs[out_idx++].get_tensor();
        const double *g_data = static_cast<double *>(TF_TensorData(g_tens.get()));
        int size = feature_sizes[key];
        int n_padded = feature_is_local[key] ? aceimpl->n_local_atoms_padded : aceimpl->n_all_atoms_padded;
        // Store in gradients map using the original key (without "grad_")
        std::copy_n(g_data, n_padded * size, &gradients[key][0]);
    }

    // Grad bond vector
    auto gf_tens = outputs[out_idx++].get_tensor();
    const double *gf_data = static_cast<double *>(TF_TensorData(gf_tens.get()));

    aceimpl->grad_bond_vector.resize(3 * aceimpl->n_neighbours_padded);
    std::copy_n(gf_data, aceimpl->n_neighbours_padded * 3, aceimpl->grad_bond_vector.data());
    grad_bv_L2 = aceimpl->grad_bond_vector;
}

void PairGRACE2LayerParallel::run_backward_layer_1() {
#ifdef GRACE_DEBUG
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE-DEBUG] Entering run_backward_layer_1\n");
#endif
    auto sig = aceimpl->model->signatures.at(backward_layer_1_name);
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
    
    auto add_input = [&](const std::string& key, const cppflow::tensor& t) {
        if (sig.inputs.count(key)) {
            inputs.emplace_back(sig.inputs.at(key).name, t);
        } else {
            error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature", key);
        }
    };

    add_input("atomic_mu_i", cppflow::tensor(aceimpl->atomic_mu_i, {aceimpl->n_all_atoms_padded}));
    add_input("atomic_mu_i_1", cppflow::tensor(aceimpl->atomic_mu_i_1, {aceimpl->n_local_atoms_padded}));
    add_input("ind_i", cppflow::tensor(aceimpl->ind_i, {aceimpl->n_neighbours_padded}));
    add_input("ind_j", cppflow::tensor(aceimpl->ind_j, {aceimpl->n_neighbours_padded}));
    add_input("mu_i", cppflow::tensor(aceimpl->mu_i, {aceimpl->n_neighbours_padded}));
    add_input("mu_j", cppflow::tensor(aceimpl->mu_j, {aceimpl->n_neighbours_padded}));
    add_input("bond_vector", cppflow::tensor(aceimpl->bond_vector, {aceimpl->n_neighbours_padded, 3}));
    add_input("batch_tot_nat_real", cppflow::tensor(std::vector<int32_t>{atom->nlocal}, {}));
    
    for (const auto& [key, shape] : feature_shapes) {
        std::string grad_key = "grad_" + key;
        std::vector<int64_t> local_shape = {aceimpl->n_local_atoms_padded};
        local_shape.insert(local_shape.end(), shape.begin(), shape.end());

        if (feature_is_local[key]) {
            add_input(grad_key, cppflow::tensor(gradients[key], local_shape));
        } else {
            // Temporarily resize to local padded size to avoid intermediate copy
            // The model expects a tensor of size n_local_atoms_padded
            size_t original_size = gradients[key].size();
            gradients[key].resize(aceimpl->n_local_atoms_padded * feature_sizes[key]);
            add_input(grad_key, cppflow::tensor(gradients[key], local_shape));
            gradients[key].resize(original_size);
        }
    }

#ifdef GRACE_DEBUG
    if (comm->me == 0) print_tensors(backward_layer_1_name, inputs);
#endif

    std::vector<std::string> out_names;
    if (sig.outputs.count(GRAD_BOND_KEY)) out_names.push_back(sig.outputs.at(GRAD_BOND_KEY).name);
    else error->all(FLERR, "[GRACE-ERROR] Key '{}' not found in signature '{}' outputs", GRAD_BOND_KEY, backward_layer_1_name);

    auto outputs = aceimpl->model->operator()(inputs, out_names);

#ifdef GRACE_DEBUG
    if (comm->me == 0) {
        std::vector<std::tuple<std::string, cppflow::tensor>> out_tensors;
        for (size_t i=0; i<outputs.size(); ++i) out_tensors.push_back({out_names[i], outputs[i]});
        print_tensors(backward_layer_1_name, out_tensors, "Output");
    }
#endif
    
    auto gf_tens = outputs[0].get_tensor();
    const double *gf_data = static_cast<double *>(TF_TensorData(gf_tens.get()));
    // Combine with L2 gradients
    for (int k = 0; k < aceimpl->n_neighbours_padded * 3; k++)
        aceimpl->grad_bond_vector[k] += gf_data[k];
}

void PairGRACE2LayerParallel::print_tensors(const std::string& name, const std::vector<std::tuple<std::string, cppflow::tensor>>& tensors, const std::string& type_prefix) {
    if (type_prefix == "Input") utils::logmesg(lmp, "[GRACE-DEBUG] Calling model in {}\n", name);
    for (const auto& tensor : tensors) {
        const auto& key = std::get<0>(tensor);
        const auto& t = std::get<1>(tensor);
        auto shape = t.shape().get_data<int64_t>();
        std::string shape_str = "[";
        for (size_t i = 0; i < shape.size(); ++i) shape_str += std::to_string(shape[i]) + (i == shape.size() - 1 ? "" : ", ");
        shape_str += "]";
        utils::logmesg(lmp, "  {}: {} | Shape: {}\n", type_prefix, key, shape_str);
        
        auto dtype = TF_TensorType(t.get_tensor().get());
        int total_elements = static_cast<int>(TF_TensorElementCount(t.get_tensor().get()));
        
        if (shape.size() >= 1) {
            int n_total_atoms = static_cast<int>(shape[0]);
            int elements_per_atom = total_elements / n_total_atoms;

            auto print_atom = [&](int i) {
                utils::logmesg(lmp, "    Atom {}:", i);
                auto log_val = [&](int k) {
                    if (dtype == TF_DOUBLE) {
                        const double *data = static_cast<const double *>(TF_TensorData(t.get_tensor().get()));
                        utils::logmesg(lmp, " {:8.4f}", data[i * elements_per_atom + k]);
                    } else if (dtype == TF_INT32) {
                        const int32_t *data = static_cast<const int32_t *>(TF_TensorData(t.get_tensor().get()));
                        utils::logmesg(lmp, " {}", data[i * elements_per_atom + k]);
                    }
                };

                if (elements_per_atom <= 10) {
                    for (int k = 0; k < elements_per_atom; k++) log_val(k);
                } else {
                    for (int k = 0; k < 5; k++) log_val(k);
                    utils::logmesg(lmp, " ...");
                    for (int k = elements_per_atom - 5; k < elements_per_atom; k++) log_val(k);
                }
                utils::logmesg(lmp, "\n");
            };

            if (n_total_atoms <= 10) {
                for (int i = 0; i < n_total_atoms; i++) print_atom(i);
            } else {
                for (int i = 0; i < 5; i++) print_atom(i);
                utils::logmesg(lmp, "    ...\n");
                for (int i = n_total_atoms - 5; i < n_total_atoms; i++) print_atom(i);
            }
        } else {
            // Scalar
             if (dtype == TF_INT32) {
                const int32_t *data = static_cast<int32_t *>(TF_TensorData(t.get_tensor().get()));
                utils::logmesg(lmp, "    Value: {}\n", data[0]);
             } else if (dtype == TF_DOUBLE) {
                const double *data = static_cast<double *>(TF_TensorData(t.get_tensor().get()));
                utils::logmesg(lmp, "    Value: {:8.4f}\n", data[0]);
             }
        }
    }
}



void *PairGRACE2LayerParallel::extract(const char *str, int &dim) {
    dim = 2; if (strcmp(str, "scale") == 0) return (void *) scale; return nullptr;
}

#endif
