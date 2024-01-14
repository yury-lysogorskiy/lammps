//
// Created by Lysogorskiy Yry on 14.01.24.
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/tp_cell/kk,PairPACETensorPotentialCellKokkos<LMPDeviceType>);
PairStyle(pace/tp_cell/kk/device,PairPACETensorPotentialCellKokkos<LMPDeviceType>);
PairStyle(pace/tp_cell/kk/host,PairPACETensorPotentialCellKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_TP_CELL_KOKKOS_H
#define LMP_PAIR_PACE_TP_CELL_KOKKOS_H

#include "kokkos_type.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"

#include "pair_pace_tp_cell.h"

namespace LAMMPS_NS {

    template<class DeviceType>
    class PairPACETensorPotentialCellKokkos : public PairPACETensorPotentialCell {

    public:

        typedef DeviceType device_type;
        typedef ArrayTypes<DeviceType> AT;
        PairPACETensorPotentialCellKokkos(class LAMMPS *);
        ~PairPACETensorPotentialCellKokkos() override;
        void compute(int, int) override;
        void coeff(int, char **) override;
        void init_style() override;
        double init_one(int, int) override;
        void allocate();

    protected:


        int host_flag;

        int eflag, vflag;

        int neighflag;

        typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
        tdual_fparams k_cutsq;
        typedef Kokkos::View<F_FLOAT**, DeviceType> t_fparams;
        t_fparams d_cutsq;

        typename AT::t_int_1d d_map;

        typename AT::t_int_1d d_atomic_mu_i_vector;
        typename AT::t_int_1d d_ind_i_vector;
        typename AT::t_int_1d d_ind_j_vector;
        typename AT::t_int_1d d_mu_i_vector;
        typename AT::t_int_1d d_mu_j_vector;

        typename AT::t_int_1d d_actual_jnum;
        typename AT::t_int_1d d_actual_jnum_shift;


        typedef Kokkos::View<double*[3], DeviceType> t_double_2d3;
        t_double_2d3 d_vector_offset;
        t_double_2d3 d_positions;
        t_double_2d3 d_fpair;

        typedef Kokkos::View<double*, DeviceType> t_double_1d;
        t_double_1d d_atom_energy;



        void grow(int, int);

    };
}    // namespace LAMMPS_NS

#endif
#endif