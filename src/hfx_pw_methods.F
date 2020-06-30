!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Test routines for HFX caclulations using PW
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!>     Made GAPW compatible 12.2019 (A. Bussy)
!>     Refactored from hfx_admm_utils (JGH)
!> \author MI
! **************************************************************************************************
MODULE hfx_pw_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        copy_info_type, cp_fm_cleanup_copy_general, cp_fm_create, cp_fm_finish_copy_general, &
        cp_fm_get_info, cp_fm_p_type, cp_fm_release, cp_fm_start_copy_general, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_type
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              twopi
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: &
        COMPLEXDATA1D, COMPLEXDATA3D, REALDATA3D, REALSPACE, RECIPROCALSPACE, pw_create, &
        pw_p_type, pw_release, pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: pw_hfx, pw_hfx_kp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_pw_methods'

CONTAINS

! **************************************************************************************************
!> \brief computes the Hartree-Fock energy brute force in a pw basis
!> \param qs_env ...
!> \param ehfx ...
!> \param hfx_section ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param irep ...
!> \par History
!>      12.2007 created [Joost VandeVondele]
!> \note
!>      only computes the HFX energy, no derivatives as yet
! **************************************************************************************************
   SUBROUTINE pw_hfx(qs_env, ehfx, hfx_section, poisson_env, auxbas_pw_pool, irep)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: ehfx
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER                                            :: irep

      CHARACTER(*), PARAMETER                            :: routineN = 'pw_hfx'

      INTEGER                                            :: blocksize, handle, iloc, iorb, &
                                                            iorb_block, ispin, iw, jloc, jorb, &
                                                            jorb_block, norb, potential_type
      LOGICAL                                            :: do_pw_hfx, explicit
      REAL(KIND=dp)                                      :: exchange_energy, fracval, omega, &
                                                            pair_energy, rcut, scaling
      REAL(KIND=dp), DIMENSION(3)                        :: kpp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_p_type)                                    :: pot_g, rho_g, rho_r
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: rho_i, rho_j
      TYPE(pw_type), POINTER                             :: greenfn
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ip_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL section_vals_val_get(hfx_section, "PW_HFX", l_val=do_pw_hfx, i_rep_section=irep)

      IF (do_pw_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fracval)
         CALL section_vals_val_get(hfx_section, "PW_HFX_BLOCKSIZE", i_val=blocksize)

         CALL get_qs_env(qs_env, mos=mo_array, pw_env=pw_env, dft_control=dft_control, &
                         cell=cell, particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)

         ! limit the blocksize by the number of orbitals
         CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=norb)
         blocksize = MAX(1, MIN(blocksize, norb))

         CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         ALLOCATE (rho_i(blocksize))
         ALLOCATE (rho_j(blocksize))

         NULLIFY (greenfn)
         CALL pw_pool_create_pw(auxbas_pw_pool, greenfn, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         ip_section => section_vals_get_subs_vals(hfx_section, "INTERACTION_POTENTIAL")
         CALL section_vals_get(ip_section, explicit=explicit)
         potential_type = do_potential_coulomb
         rcut = 0.0_dp
         omega = 0.0_dp
         kpp = 0.0_dp
         IF (explicit) THEN
            CALL section_vals_val_get(ip_section, "POTENTIAL_TYPE", i_val=potential_type)
            CALL section_vals_val_get(ip_section, "CUTOFF_RADIUS", r_val=rcut)
            CALL section_vals_val_get(ip_section, "OMEGA", r_val=omega)
         ENDIF
         grid => poisson_env%green_fft%influence_fn%pw_grid
         CALL set_greenfn(greenfn, grid, potential_type, kpp, kpp, rcut, omega)

         DO iorb_block = 1, blocksize
            NULLIFY (rho_i(iorb_block)%pw)
            CALL pw_create(rho_i(iorb_block)%pw, rho_r%pw%pw_grid, &
                           use_data=REALDATA3D, &
                           in_space=REALSPACE)
            NULLIFY (rho_j(iorb_block)%pw)
            CALL pw_create(rho_j(iorb_block)%pw, rho_r%pw%pw_grid, &
                           use_data=REALDATA3D, &
                           in_space=REALSPACE)
         ENDDO

         exchange_energy = 0.0_dp

         DO ispin = 1, SIZE(mo_array)
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)

            IF (mo_array(ispin)%mo_set%use_mo_coeff_b) THEN !fm->dbcsr
               CALL copy_dbcsr_to_fm(mo_coeff_b, mo_coeff) !fm->dbcsr
            ENDIF !fm->dbcsr

            CALL cp_fm_get_info(mo_coeff, ncol_global=norb)

            DO iorb_block = 1, norb, blocksize

               DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)

                  iloc = iorb - iorb_block + 1
                  CALL calculate_wavefunction(mo_coeff, iorb, rho_i(iloc), rho_g, &
                                              atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                              pw_env)

               ENDDO

               DO jorb_block = iorb_block, norb, blocksize

                  DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)

                     jloc = jorb - jorb_block + 1
                     CALL calculate_wavefunction(mo_coeff, jorb, rho_j(jloc), rho_g, &
                                                 atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                                 pw_env)

                  ENDDO

                  DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)
                     iloc = iorb - iorb_block + 1
                     DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)
                        jloc = jorb - jorb_block + 1
                        IF (jorb < iorb) CYCLE

                        ! compute the pair density
                        rho_r%pw%cr3d = rho_i(iloc)%pw%cr3d*rho_j(jloc)%pw%cr3d

                        ! go the g-space and compute hartree energy
                        CALL pw_transfer(rho_r%pw, rho_g%pw)
                        CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw, &
                                              greenfn=greenfn)

                        ! sum up to the full energy
                        scaling = fracval
                        IF (SIZE(mo_array) == 1) scaling = scaling*2.0_dp
                        IF (iorb /= jorb) scaling = scaling*2.0_dp

                        exchange_energy = exchange_energy - scaling*pair_energy

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO

         DO iorb_block = 1, blocksize
            CALL pw_release(rho_i(iorb_block)%pw)
            CALL pw_release(rho_j(iorb_block)%pw)
         ENDDO

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, greenfn)

         iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                   extension=".scfLog")

         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10))") &
               "HF_PW_HFX| PW exchange energy:", exchange_energy
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10),/)") &
               "HF_PW_HFX| Gaussian exchange energy:", ehfx
         ENDIF

         CALL cp_print_key_finished_output(iw, logger, hfx_section, "HF_INFO")

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pw_hfx

! **************************************************************************************************
!> \brief computes the Hartree-Fock energy brute force in a pw basis for k-point calculations
!> \param qs_env ...
!> \param ehfx ...
!> \param hfx_section ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param irep ...
!> \par History
!>      12.2007 created [Joost VandeVondele]
!>      06.2020 kpoints [JGH]
!> \note
!>      only computes the HFX energy, no derivatives as yet
! **************************************************************************************************
   SUBROUTINE pw_hfx_kp(qs_env, ehfx, hfx_section, poisson_env, auxbas_pw_pool, irep)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: ehfx
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER                                            :: irep

      CHARACTER(*), PARAMETER :: routineN = 'pw_hfx_kp', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: debug_level = 1

      INTEGER :: blocksize, handle, ic, ik, ik1, ik2, ikk, iloc, indx, iorb, iorb_block, iounit, &
         is, ispin, iw, jloc, jorb, jorb_block, nc, ncol, nkp, norb, nrow, nspin, potential_type
      LOGICAL                                            :: do_pw_hfx, explicit, full_grid, &
                                                            my_kpgrp, real_only
      REAL(KIND=dp)                                      :: ekpp, exchange_energy, fracval, omega, &
                                                            pair_energy, rcut, scaling, wk1, wk2
      REAL(KIND=dp), DIMENSION(3)                        :: kp1, kp2
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:)    :: info
      TYPE(cp_fm_p_type), DIMENSION(:, :, :), POINTER    :: fmwork
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fmdummy, mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_p_type)                                    :: pot_g, rho_g, rho_r
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:, :)      :: rho_i, rho_j
      TYPE(pw_type), POINTER                             :: greenfn
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ip_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL section_vals_val_get(hfx_section, "PW_HFX", l_val=do_pw_hfx, i_rep_section=irep)

      IF (do_pw_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fracval)
         CALL section_vals_val_get(hfx_section, "PW_HFX_BLOCKSIZE", i_val=blocksize)

         CALL get_qs_env(qs_env, mos=mo_array, pw_env=pw_env, dft_control=dft_control, &
                         kpoints=kpoint, cell=cell, particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
         para_env => kpoint%blacs_env_all%para_env

         kp => kpoint%kp_env(1)%kpoint_env
         nspin = SIZE(kp%mos, 2)
         nc = SIZE(kp%mos, 1)
         real_only = kpoint%use_real_wfn
         IF (real_only) THEN
            CPASSERT(nc == 1)
         END IF
         full_grid = kpoint%full_grid
         nkp = kpoint%nkp
         ALLOCATE (info(nspin*nkp*nc))
         NULLIFY (fmdummy)

         ! limit the blocksize by the number of orbitals
         CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=norb)
         blocksize = MAX(1, MIN(blocksize, norb))

         CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                use_data=COMPLEXDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         ALLOCATE (rho_i(nc, blocksize))
         ALLOCATE (rho_j(nc, blocksize))

         DO ic = 1, nc
            DO iorb_block = 1, blocksize
               NULLIFY (rho_i(ic, iorb_block)%pw)
               CALL pw_create(rho_i(ic, iorb_block)%pw, rho_r%pw%pw_grid, &
                              use_data=REALDATA3D, &
                              in_space=REALSPACE)
               NULLIFY (rho_j(ic, iorb_block)%pw)
               CALL pw_create(rho_j(ic, iorb_block)%pw, rho_r%pw%pw_grid, &
                              use_data=REALDATA3D, &
                              in_space=REALSPACE)
            ENDDO
         ENDDO

         exchange_energy = 0.0_dp

         IF (debug_level > 0) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(A)") "PW_HFX_KP (DEBUG): Start of subroutine"
               WRITE (iounit, "(A,T72,I9)") "PW_HFX_KP (DEBUG): Number of kpoints active ", nkp
               WRITE (iounit, "(A,T72,I9)") "PW_HFX_KP (DEBUG): Number of spins ", nspin
               WRITE (iounit, "(A,T80,L1)") "PW_HFX_KP (DEBUG): Real WFN ", real_only
               WRITE (iounit, "(A,T72,I9)") "PW_HFX_KP (DEBUG): Blocksize ", blocksize
               WRITE (iounit, "(A,T72,F9.5)") "PW_HFX_KP (DEBUG): HFX fraction", fracval
               WRITE (iounit, "(A)") "PW_HFX_KP (DEBUG): Start redistributing wavefunctions"
            END IF
         END IF

         CALL cp_fm_get_info(mo_coeff, matrix_struct=matrix_struct)
         ALLOCATE (fmwork(nc, nkp, nspin))
         DO is = 1, nspin
            DO ik = 1, nkp
               DO ic = 1, nc
                  CALL cp_fm_create(fmwork(ic, ik, is)%matrix, matrix_struct)
               END DO
            END DO
         END DO

         DO ispin = 1, nspin
            indx = 0
            DO ik = 1, nkp
               my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
               IF (my_kpgrp) THEN
                  ikk = ik - kpoint%kp_range(1) + 1
                  kp => kpoint%kp_env(ikk)%kpoint_env
               ELSE
                  NULLIFY (kp)
               END IF
               IF (my_kpgrp) THEN
                  DO ic = 1, nc
                     mo_set => kp%mos(ic, ispin)%mo_set
                     indx = indx + 1
                     CALL cp_fm_start_copy_general(mo_set%mo_coeff, fmwork(ic, ik, ispin)%matrix, &
                                                   para_env, info(indx))
                  END DO
               ELSE
                  DO ic = 1, nc
                     indx = indx + 1
                     CALL cp_fm_start_copy_general(fmdummy, fmwork(ic, ik, ispin)%matrix, &
                                                   para_env, info(indx))
                  END DO
               END IF
            END DO
         END DO
         ! Finish communication
         DO ispin = 1, nspin
            indx = 0
            DO ik = 1, nkp
               DO ic = 1, nc
                  indx = indx + 1
                  CALL cp_fm_finish_copy_general(fmwork(ic, ik, ispin)%matrix, info(indx))
               END DO
            END DO
         END DO
         ! Clean up communication
         DO ispin = 1, nspin
            indx = 0
            DO ik = 1, nkp
               my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
               IF (my_kpgrp) THEN
                  ikk = ik - kpoint%kp_range(1) + 1
                  kp => kpoint%kp_env(ikk)%kpoint_env
                  DO ic = 1, nc
                     mo_set => kp%mos(ic, ispin)%mo_set
                     indx = indx + 1
                     CALL cp_fm_cleanup_copy_general(mo_set%mo_coeff, info(indx))
                  END DO
               ELSE
                  ! calls with dummy arguments, so not included
                  ! therefore just increment counter by trip count
                  indx = indx + nc
               END IF
            END DO
         END DO
         DEALLOCATE (info)
         IF (debug_level > 0) THEN
            CALL cp_fm_get_info(mo_coeff, nrow_global=nrow, ncol_global=ncol)
            IF (iounit > 0) THEN
               wk1 = 8._dp*REAL(nspin*nc*nkp, KIND=dp)*REAL(nrow*ncol, KIND=dp)
               WRITE (iounit, "(A,T69,F12.0)") "PW_HFX_KP (DEBUG): Local memory needed for wavefunctions (bytes)", wk1
               WRITE (iounit, "(A)") "PW_HFX_KP (DEBUG): End redistributing wavefunctions"
            END IF
         END IF

         DO ispin = 1, nspin

            DO ik1 = 1, nkp
               kp1 = kpoint%xkp(1:3, ik1)*twopi
               wk1 = kpoint%wkp(ik1)
               DO ik2 = 1, nkp
                  kp2 = kpoint%xkp(1:3, ik2)*twopi
                  wk2 = kpoint%wkp(ik2)
                  IF (debug_level > 0) THEN
                     grid => poisson_env%green_fft%influence_fn%pw_grid
                     IF (iounit > 0) THEN
                        WRITE (iounit, "(A,T65,I2,2I7)") "PW_HFX_KP (DEBUG): spin, kp1, kp2", ispin, ik1, ik2
                        WRITE (iounit, "(A,T32,3F10.5,T71,F10.5)") "PW_HFX_KP (DEBUG): KP1", kp1, wk1
                        WRITE (iounit, "(A,T32,3F10.5,T71,F10.5)") "PW_HFX_KP (DEBUG): KP2", kp2, wk2
                        WRITE (iounit, "(A,T51,2I15)") "PW_HFX_KP (DEBUG): Number of PW grid points", grid%ngpts, grid%ngpts_cut
                        WRITE (iounit, "(A,T54,3I9)") "PW_HFX_KP (DEBUG): Number of RS grid points", grid%npts
                     END IF
                  END IF
                  ekpp = 0.0_dp

                  NULLIFY (greenfn)
                  CALL pw_pool_create_pw(auxbas_pw_pool, greenfn, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  ip_section => section_vals_get_subs_vals(hfx_section, "INTERACTION_POTENTIAL")
                  CALL section_vals_get(ip_section, explicit=explicit)
                  potential_type = do_potential_coulomb
                  rcut = 0.0_dp
                  omega = 0.0_dp
                  IF (explicit) THEN
                     CALL section_vals_val_get(ip_section, "POTENTIAL_TYPE", i_val=potential_type)
                     CALL section_vals_val_get(ip_section, "CUTOFF_RADIUS", r_val=rcut)
                     CALL section_vals_val_get(ip_section, "OMEGA", r_val=omega)
                  ENDIF
                  grid => poisson_env%green_fft%influence_fn%pw_grid
                  CALL set_greenfn(greenfn, grid, potential_type, kp1, kp2, rcut, omega)

                  mo_coeff => fmwork(1, 1, ispin)%matrix
                  CALL cp_fm_get_info(mo_coeff, ncol_global=norb)

                  DO iorb_block = 1, norb, blocksize

                     DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)
                        iloc = iorb - iorb_block + 1
                        DO ic = 1, nc
                           mo_coeff => fmwork(ic, ik1, ispin)%matrix
                           CALL calculate_wavefunction(mo_coeff, iorb, rho_i(ic, iloc), rho_g, &
                                                       atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                                       pw_env)
                        END DO
                     END DO

                     DO jorb_block = iorb_block, norb, blocksize

                        DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)
                           jloc = jorb - jorb_block + 1
                           DO ic = 1, nc
                              mo_coeff => fmwork(ic, ik2, ispin)%matrix
                              CALL calculate_wavefunction(mo_coeff, jorb, rho_j(ic, jloc), rho_g, &
                                                          atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                                          pw_env)
                           END DO
                        END DO

                        DO iorb = iorb_block, MIN(iorb_block + blocksize - 1, norb)
                           iloc = iorb - iorb_block + 1
                           DO jorb = jorb_block, MIN(jorb_block + blocksize - 1, norb)
                              jloc = jorb - jorb_block + 1
                              IF (jorb < iorb .AND. ik1 == ik2) CYCLE

                              ! compute the pair density
                              IF (real_only) THEN
                                 rho_r%pw%cc3d = CMPLX(rho_i(1, iloc)%pw%cr3d*rho_j(1, jloc)%pw%cr3d, 0.0_dp, KIND=dp)
                              ELSE
                                 rho_r%pw%cc3d = CMPLX( &
                                                 rho_i(1, iloc)%pw%cr3d*rho_j(1, jloc)%pw%cr3d- &
                                                 rho_i(2, iloc)%pw%cr3d*rho_j(2, jloc)%pw%cr3d, &
                                                 rho_i(1, iloc)%pw%cr3d*rho_j(2, jloc)%pw%cr3d+ &
                                                 rho_i(2, iloc)%pw%cr3d*rho_j(1, jloc)%pw%cr3d, &
                                                 KIND=dp)
                              END IF

                              ! go the g-space and compute hartree energy
                              CALL pw_transfer(rho_r%pw, rho_g%pw)
                              CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw, &
                                                    greenfn=greenfn)

                              ! sum up to the full energy
                              scaling = fracval*wk1*wk2
                              IF (nspin == 1) scaling = scaling*2.0_dp
                              IF (iorb /= jorb .AND. ik1 == ik2) scaling = scaling*2.0_dp

                              exchange_energy = exchange_energy - scaling*pair_energy
                              ekpp = ekpp - scaling*pair_energy
                              IF (debug_level > 1) THEN
                                 IF (iounit > 0) THEN
                                    WRITE (iounit, "(A,T50,2I5,T63,F18.12)") "PW_HFX_KP (DEBUG): Pair energy", &
                                       iorb, jorb, -scaling*pair_energy
                                    WRITE (iounit, "(A,T45,2F18.12)") "PW_HFX_KP (DEBUG): Pair density", &
                                       SUM(rho_r%pw%cc3d)*rho_r%pw%pw_grid%dvol
                                 END IF
                              END IF

                           ENDDO
                        ENDDO

                     ENDDO
                  ENDDO
                  IF (debug_level > 0) THEN
                     IF (iounit > 0) THEN
                        WRITE (iounit, "(A,T50,2I5,T63,F18.12)") "PW_HFX_KP (DEBUG): KP energy", &
                           ik1, ik2, ekpp
                     END IF
                  END IF
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, greenfn)
               ENDDO
            ENDDO
         ENDDO

         DO ic = 1, nc
            DO iorb_block = 1, blocksize
               CALL pw_release(rho_i(ic, iorb_block)%pw)
               CALL pw_release(rho_j(ic, iorb_block)%pw)
            ENDDO
         ENDDO

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)

         IF (ASSOCIATED(fmwork)) THEN
            DO is = 1, nspin
               DO ik = 1, nkp
                  DO ic = 1, nc
                     CALL cp_fm_release(fmwork(ic, ik, is)%matrix)
                  END DO
               END DO
            END DO
            DEALLOCATE (fmwork)
         END IF

         iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                   extension=".scfLog")

         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10))") &
               "HF_PW_HFX| PW exchange energy(kpoints):", exchange_energy
            WRITE (UNIT=iw, FMT="((T3,A,T61,F20.10),/)") &
               "HF_PW_HFX| Gaussian exchange energy(kpoints):", ehfx
         ENDIF

         CALL cp_print_key_finished_output(iw, logger, hfx_section, "HF_INFO")

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pw_hfx_kp

! **************************************************************************************************
!> \brief ...
!> \param greenfn ...
!> \param grid ...
!> \param potential_type ...
!> \param kp1 ...
!> \param kp2 ...
!> \param rcut ...
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE set_greenfn(greenfn, grid, potential_type, kp1, kp2, rcut, omega)
      TYPE(pw_type), POINTER                             :: greenfn
      TYPE(pw_grid_type), POINTER                        :: grid
      INTEGER, INTENT(IN)                                :: potential_type
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kp1, kp2
      REAL(KIND=dp), INTENT(IN)                          :: rcut, omega

      CHARACTER(*), PARAMETER :: routineN = 'set_greenfn', routineP = moduleN//':'//routineN

      INTEGER                                            :: ig
      REAL(KIND=dp)                                      :: g2, g3d, gg, gx, gy, gz, o2

      IF (potential_type == do_potential_coulomb) THEN
         DO ig = grid%first_gne0, grid%ngpts_cut_local
            gx = grid%g(1, ig) + kp1(1) - kp2(1)
            gy = grid%g(2, ig) + kp1(2) - kp2(2)
            gz = grid%g(3, ig) + kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            greenfn%cc(ig) = fourpi/g2
         END DO
         IF (grid%have_g0) THEN
            gx = kp1(1) - kp2(1)
            gy = kp1(2) - kp2(2)
            gz = kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            IF (g2 > 0.0_dp) THEN
               greenfn%cc(1) = fourpi/g2
            ELSE
               greenfn%cc(1) = 0.0_dp
            END IF
         END IF
      ELSEIF (potential_type == do_potential_truncated) THEN
         DO ig = grid%first_gne0, grid%ngpts_cut_local
            gx = grid%g(1, ig) + kp1(1) - kp2(1)
            gy = grid%g(2, ig) + kp1(2) - kp2(2)
            gz = grid%g(3, ig) + kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            gg = SQRT(g2)
            g3d = fourpi/g2
            greenfn%cc(ig) = g3d*(1.0_dp - COS(rcut*gg))
         END DO
         IF (grid%have_g0) THEN
            gx = kp1(1) - kp2(1)
            gy = kp1(2) - kp2(2)
            gz = kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            IF (g2 > 0.0_dp) THEN
               gg = SQRT(g2)
               g3d = fourpi/g2
               greenfn%cc(1) = g3d*(1.0_dp - COS(rcut*gg))
            ELSE
               greenfn%cc(1) = 0.5_dp*fourpi*rcut*rcut
            END IF
         END IF
      ELSEIF (potential_type == do_potential_short) THEN
         o2 = 0.0_dp
         IF (omega > 0.0_dp) o2 = 0.25_dp/(omega*omega)
         DO ig = grid%first_gne0, grid%ngpts_cut_local
            gx = grid%g(1, ig) + kp1(1) - kp2(1)
            gy = grid%g(2, ig) + kp1(2) - kp2(2)
            gz = grid%g(3, ig) + kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            gg = -o2*g2
            g3d = fourpi/g2
            greenfn%cc(ig) = g3d*(1.0_dp - EXP(gg))
         END DO
         IF (grid%have_g0) THEN
            gx = kp1(1) - kp2(1)
            gy = kp1(2) - kp2(2)
            gz = kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            IF (g2 > 0.0_dp) THEN
               gg = -o2*g2
               g3d = fourpi/g2
               greenfn%cc(1) = g3d*(1.0_dp - EXP(gg))
            ELSE
               greenfn%cc(1) = 0.0_dp
            END IF
         END IF
      ELSE
         CPWARN("PW_SCF: Potential type not supported, calculation uses Coulomb potential.")
         DO ig = grid%first_gne0, grid%ngpts_cut_local
            gx = grid%g(1, ig) + kp1(1) - kp2(1)
            gy = grid%g(2, ig) + kp1(2) - kp2(2)
            gz = grid%g(3, ig) + kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            greenfn%cc(ig) = fourpi/g2
         END DO
         IF (grid%have_g0) THEN
            gx = kp1(1) - kp2(1)
            gy = kp1(2) - kp2(2)
            gz = kp1(3) - kp2(3)
            g2 = gx*gx + gy*gy + gz*gz
            IF (g2 > 0.0_dp) THEN
               greenfn%cc(1) = fourpi/g2
            ELSE
               greenfn%cc(1) = 0.0_dp
            END IF
         END IF
      END IF
   END SUBROUTINE set_greenfn

END MODULE hfx_pw_methods
