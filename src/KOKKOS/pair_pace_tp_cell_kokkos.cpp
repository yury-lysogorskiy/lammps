//
// Created by Lysogorskiy Yry on 14.01.24.
//

#include "pair_pace_tp_cell_kokkos.h"


#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

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

// Also used as reference: https://github.com/ACEsuit/lammps/blob/mace/src/KOKKOS/pair_mace_kokkos.cpp

/* ---------------------------------------------------------------------- */
using namespace LAMMPS_NS;

template<class DeviceType>
PairPACETensorPotentialCellKokkos<DeviceType>::PairPACETensorPotentialCellKokkos(LAMMPS *lmp) : PairPACETensorPotentialCell(lmp)
{
    respa_enable = 0;

    kokkosable = 1;
    atomKK = (AtomKokkos *) atom;
    execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
    datamask_read = EMPTY_MASK;
    datamask_modify = EMPTY_MASK;

    host_flag = (execution_space == Host);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairPACETensorPotentialCellKokkos<DeviceType>::~PairPACETensorPotentialCellKokkos()
{
    if (copymode) return;

}
/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairPACETensorPotentialCellKokkos<DeviceType>::grow(int totnatom, int totneigh)
{
    if ((int)d_atomic_mu_i_vector.extent(0) < totnatom) {
        MemKK::realloc_kokkos(d_atomic_mu_i_vector, "pace/tp_cell/kk:d_atomic_mu_i_vector", totnatom);
        MemKK::realloc_kokkos(d_positions, "pace/tp_cell/kk:d_positions", totnatom,3);


        MemKK::realloc_kokkos(d_atom_energy, "pace/tp_cell/kk:d_atom_energy", totnatom);
    }

    if ((int)d_ind_i_vector.extent(0) < totneigh ) {
        MemKK::realloc_kokkos(d_ind_i_vector, "pace/tp_cell/kk:d_ind_i_vector", totneigh);
        MemKK::realloc_kokkos(d_ind_j_vector, "pace/tp_cell/kk:d_ind_j_vector", totneigh);
        MemKK::realloc_kokkos(d_mu_i_vector, "pace/tp_cell/kk:d_mu_i_vector", totneigh);
        MemKK::realloc_kokkos(d_mu_j_vector, "pace/tp_cell/kk:d_mu_j_vector", totneigh);

        MemKK::realloc_kokkos(d_vector_offset, "pace/tp_cell/kk:d_vector_offset", totneigh,3);
        MemKK::realloc_kokkos(d_fpair, "pace/tp_cell/kk:d_fpair", totneigh,3);
    }
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACETensorPotentialCellKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
    eflag = eflag_in;
    vflag = vflag_in;

    if (neighflag == FULL) no_virial_fdotr_compute = 1;
    ev_init(eflag,vflag,0);

    // if (eflag_atom) {
    //     // memoryKK->destroy_kokkos(k_eatom,eatom);
    //     // memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    //     // d_eatom = k_eatom.view<DeviceType>();
    // }

    atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK|TAG_MASK);
    printf("LABEL 1\n");

    NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
    auto d_numneigh = k_list->d_numneigh;
    auto d_neighbors = k_list->d_neighbors;
    auto d_ilist = k_list->d_ilist;

    if (atom->nlocal != list->inum)
        error->all(FLERR, "ERROR: nlocal != inum.");

    if (eflag_atom || vflag_atom)
        error->all(FLERR, "ERROR: mace/kokkos eflag_atom and/or vflag_atom not implemented.");

    int nlocal = atom->nlocal;
    int n_fake_atoms = (do_padding ? 1 : 0);
    int n_real_neighbours=0;
    int n_fake_neighbours=0;

    // atom map
    auto map_style = atom->map_style;
    auto k_map_array = atomKK->k_map_array;
    auto k_map_hash = atomKK->k_map_hash;
    k_map_array.template sync<DeviceType>();

    auto x = atomKK->k_x.view<DeviceType>();
    auto f = atomKK->k_f.view<DeviceType>();
    auto tag = atomKK->k_tag.view<DeviceType>();
    auto type = atomKK->k_type.view<DeviceType>();

    data_timer.start();
    //count number of real atoms
    if((int)d_actual_jnum.extent(0)<nlocal) {
        printf("LABEL grow d_actual_jnum, d_actual_jnum_shift to nlocal=%d\n",nlocal);
        MemKK::realloc_kokkos(d_actual_jnum, "pace/tp_cell/kk:d_actual_jnum", nlocal);
        MemKK::realloc_kokkos(d_actual_jnum_shift, "pace/tp_cell/kk:d_actual_jnum_shift", nlocal);
    }
    printf("LABEL 2\n");

    // d_actual_jnum = Kokkos::View<int64_t*,DeviceType>("d_actual_jnum", nlocal); // grow ?
    double cutoff_sq = cutoff * cutoff;
    {
        Kokkos::parallel_for("PairPACETensorPotentialCellKokkos: Fill d_actual_jnum.", nlocal, KOKKOS_LAMBDA (const int ii) {
            const int i = d_ilist(ii);
            const double xtmp = x(i,0);
            const double ytmp = x(i,1);
            const double ztmp = x(i,2);
            for (int jj=0; jj<d_numneigh(i); ++jj) {
              int j = d_neighbors(i,jj);
              j &= NEIGHMASK;
              const double delx = xtmp - x(j,0);
              const double dely = ytmp - x(j,1);
              const double delz = ztmp - x(j,2);
              const double rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutoff_sq) {
                d_actual_jnum(ii) += 1;
              }
            }
        });
    }
    printf("LABEL 3\n");
    // count number of real neighbours and build shift vector d_actual_jnum_shift
    // d_actual_jnum_shift = Kokkos::View<int64_t*,DeviceType>("d_actual_jnum_shift", nlocal); // grow ?
    {
        n_real_neighbours = 0;
        Kokkos::parallel_scan("PairPACETensorPotentialCellKokkos: Sum real neighbours", nlocal,
            KOKKOS_LAMBDA(const int ii, int& n_neigh, bool is_final) {
                if(is_final) d_actual_jnum_shift(ii) = n_neigh;
                n_neigh += d_actual_jnum(ii);
            },
        n_real_neighbours);
    }

    printf("LABEL 4\n");
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

    printf("LABEL 5\n");
    grow(nlocal+n_fake_atoms, tot_neighbours);

    printf("LABEL 6\n");
    // atomic_mu_i: per-atom species type + padding with type[0]
    {
        Kokkos::parallel_for("PairPACETensorPotentialCellKokkos:d_atomic_mu_i_vector", nlocal, KOKKOS_LAMBDA(const int ii) {
            const int i = d_ilist(ii);
            d_atomic_mu_i_vector(i) =  d_map(type(i));
        });
    }

    printf("LABEL 7\n");
    // fill bond-related vectors: d_ind_i_vector, d_ind_j_vector, d_mu_i_vector, d_mu_j_vector, d_vector_offset
    {
        Kokkos::parallel_for("PairPACETensorPotentialCellKokkos: Fill bond-related vectors", nlocal, KOKKOS_LAMBDA(const int ii) {
            const int i = d_ilist(ii);
            const double xtmp = x(i,0);
            const double ytmp = x(i,1);
            const double ztmp = x(i,2);
            int tot_ind = d_actual_jnum_shift(ii);
            for (int jj=0; jj<d_numneigh(i); ++jj) {
              int j = d_neighbors(i,jj);
              j &= NEIGHMASK;
              const double delx = xtmp - x(j,0);
              const double dely = ytmp - x(j,1);
              const double delz = ztmp - x(j,2);
              const double rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutoff_sq) {
                    d_ind_i_vector(tot_ind) = i;
                    // remap j to j_local
                    int j_local = AtomKokkos::map_kokkos<DeviceType>(tag(j),map_style,k_map_array,k_map_hash);
                    d_ind_j_vector(tot_ind) = j_local;

                    d_mu_i_vector(tot_ind) = d_map(type(i));
                    d_mu_j_vector(tot_ind) = d_map(type(j));

                    double shiftx = x(j,0) - x(j_local,0);
                    double shifty = x(j,1) - x(j_local,1);
                    double shiftz = x(j,2) - x(j_local,2);

                    d_vector_offset(tot_ind,0) = shiftx;
                    d_vector_offset(tot_ind,1) = shifty;
                    d_vector_offset(tot_ind,2) = shiftz;
                    ++tot_ind;
              }
            }
        });
        // add fake bonds
        int fake_atom_ind = nlocal + n_fake_atoms - 1;
        Kokkos::parallel_for("PairPACETensorPotentialCellKokkos: Fake fill bond-related vectors", n_fake_neighbours, KOKKOS_LAMBDA(const int ii) {
            d_ind_i_vector(n_real_neighbours+ii) = fake_atom_ind;
            d_ind_j_vector(n_real_neighbours+ii) = fake_atom_ind;
            d_mu_i_vector(n_real_neighbours+ii) = 0;
            d_mu_j_vector(n_real_neighbours+ii) = 0;
        });
    }

    printf("LABEL 8\n");
    {
        Kokkos::parallel_for("PairMACEKokkos: Fill d_positions.", nlocal, KOKKOS_LAMBDA (const int i) {
           d_positions(i,0) = x(i,0);
           d_positions(i,1) = x(i,1);
           d_positions(i,2) = x(i,2);
        });

    }

    // Currently, we will transfer data to HOST
    // this will decrease performance due to GPU-CPU data transfer

    printf("LABEL copy d2h\n");
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    auto h_atomic_mu_i_vector = Kokkos::create_mirror_view(d_atomic_mu_i_vector);
    Kokkos::deep_copy(h_atomic_mu_i_vector, d_atomic_mu_i_vector);
    cppflow::tensor atomic_mu_i_vector_tens(TF_DOUBLE, h_atomic_mu_i_vector.data(),
             (nlocal + n_fake_atoms)*sizeof(int), {nlocal + n_fake_atoms});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "atomic_mu_i" + ":0", atomic_mu_i_vector_tens);

    // batch_nat = number of extened atoms + padding
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nat" + ":0",
                        cppflow::tensor(std::vector<int32_t>{nlocal + n_fake_atoms}, {}));

    // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "batch_nreal_atoms_per_structure" + ":0",
                        cppflow::tensor(std::vector<int32_t>{nlocal + n_fake_atoms}, {1}));


    auto h_ind_i_vector = Kokkos::create_mirror_view(d_ind_i_vector);
    Kokkos::deep_copy(h_ind_i_vector, d_ind_i_vector);
    cppflow::tensor ind_i_tens(TF_DOUBLE, h_ind_i_vector.data(), tot_neighbours*sizeof(int), {tot_neighbours});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_i" + ":0", ind_i_tens);

    auto h_ind_j_vector = Kokkos::create_mirror_view(d_ind_j_vector);
    Kokkos::deep_copy(h_ind_j_vector, d_ind_j_vector);
    cppflow::tensor ind_j_tens(TF_DOUBLE, h_ind_j_vector.data(), tot_neighbours*sizeof(int), {tot_neighbours});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "ind_j" + ":0", ind_j_tens);

    auto h_mu_i_vector = Kokkos::create_mirror_view(d_mu_i_vector);
    Kokkos::deep_copy(h_mu_i_vector, d_mu_i_vector);
    cppflow::tensor mu_i_tens(TF_DOUBLE, h_mu_i_vector.data(), tot_neighbours*sizeof(int), {tot_neighbours});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_i" + ":0", mu_i_tens);

    auto h_mu_j_vector = Kokkos::create_mirror_view(d_mu_j_vector);
    Kokkos::deep_copy(h_mu_j_vector, d_mu_j_vector);
    cppflow::tensor mu_j_tens(TF_DOUBLE, h_mu_j_vector.data(), tot_neighbours*sizeof(int), {tot_neighbours});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "mu_j" + ":0", mu_j_tens);

    // num_struc: 1
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "num_struc" + ":0",
                        cppflow::tensor(std::vector<int32_t>{1}, {}));


    // vector_offsets: 1
    auto h_vector_offset = Kokkos::create_mirror_view(d_vector_offset);
    Kokkos::deep_copy(h_vector_offset, d_vector_offset);
    cppflow::tensor vector_offset_tens(TF_DOUBLE, h_vector_offset.data(), tot_neighbours*3*sizeof(int), {tot_neighbours, 3});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "vector_offsets" + ":0", vector_offset_tens);


    // positions: [n_extened_atoms+n_fake_atoms, 3]
    auto h_positions = Kokkos::create_mirror_view(d_positions);
    Kokkos::deep_copy(h_positions, d_positions);
    cppflow::tensor pos_tens(TF_DOUBLE, h_positions.data(),
            (nlocal + n_fake_atoms)*3*sizeof(double), {nlocal + n_fake_atoms,3});
    inputs.emplace_back(DEFAULT_INPUT_PREFIX + "positions" + ":0", pos_tens);
    data_timer.stop();

    tp_timer.start();
    //CALL MODEL
    printf("LABEL CALL MODEL\n");
    std::vector<cppflow::tensor> output = aceimpl->model->operator()(
            inputs,
            {
                    "StatefulPartitionedCall:0", // energy [nat]
                    "StatefulPartitionedCall:1", // force-per-bond [n_neigh, 3]
            }
    );
    tp_timer.stop();




    if (eflag_global) {
        auto &e_out = output[0];
        auto e_tens = e_out.get_tensor();
        double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));

        auto h_atom_energy = Kokkos::View<double*,Kokkos::LayoutRight,LMPHostType,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>(e_data, nlocal+n_fake_atoms);

        Kokkos::deep_copy(d_atom_energy, h_atom_energy);

        eng_vdwl = 0.0;
        Kokkos::parallel_reduce("PairPACETensorPotentialCellKokkos: Accumulate site energies.", nlocal,  KOKKOS_LAMBDA(const int ii, double &eng_vdwl) {
              const int i = d_ilist(ii);
              eng_vdwl += d_atom_energy(i);
        }, eng_vdwl);
    }

    auto &f_out = output[1];
    auto f_tens = f_out.get_tensor();
    double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));

    auto h_fpair = Kokkos::View<double*[3],Kokkos::LayoutRight,LMPHostType,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>(f_data, nlocal+n_fake_atoms,3);
    Kokkos::deep_copy(d_fpair, h_fpair);

}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACETensorPotentialCellKokkos<DeviceType>::coeff(int narg, char **arg)
{
    if (!allocated) allocate();
    PairPACETensorPotentialCell::coeff(narg,arg);

    // Set up element lists

    auto h_map = Kokkos::create_mirror_view(d_map);

    for (int i = 1; i <= atom->ntypes; i++)
        h_map(i) = map[i];

    Kokkos::deep_copy(d_map,h_map);
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairPACETensorPotentialCellKokkos<DeviceType>::init_style()
{
    PairPACETensorPotentialCell::init_style();
    auto request = neighbor->find_request(this);
    request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                             !std::is_same<DeviceType,LMPDeviceType>::value);
    request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
}

template<class DeviceType>
double PairPACETensorPotentialCellKokkos<DeviceType>::init_one(int i, int j)
{
    double cutone = PairPACETensorPotentialCell::init_one(i,j);
    k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
    k_cutsq.template modify<LMPHostType>();
    return cutone;
}

template<class DeviceType>
void PairPACETensorPotentialCellKokkos<DeviceType>::allocate()
{
    PairPACETensorPotentialCell::allocate();

    int n = atom->ntypes + 1;
    MemKK::realloc_kokkos(d_map, "pace/tp_cell/kk:map", n);

    MemKK::realloc_kokkos(k_cutsq, "pace/tp_cell/kk:cutsq", n, n);
    d_cutsq = k_cutsq.template view<DeviceType>();
}

namespace LAMMPS_NS {
    template class PairPACETensorPotentialCellKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
    template class PairPACETensorPotentialCellKokkos<LMPHostType>;
#endif
}
