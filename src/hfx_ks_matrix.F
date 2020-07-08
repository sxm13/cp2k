!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief HFX part of KS matrix
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!>     Made GAPW compatible 12.2019 (A. Bussy)
!>     refactoring 06-2020 [JGH]
!> \author MI
! **************************************************************************************************
MODULE hfx_ks_matrix
   USE admm_methods,                    ONLY: scale_dm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_pw_methods,                  ONLY: pw_hfx,&
                                              pw_hfx_kp
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: rt_prop_type
   USE virial_types,                    ONLY: virial_type
   USE xc_adiabatic_utils,              ONLY: rescale_xc_potential
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: hfx_ks_matrix_gamma, hfx_ks_matrix_kp, tddft_hfx_matrix

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ks_matrix'

CONTAINS

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param qs_env ...
!> \param matrix_ks ...
!> \param rho ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \par History
!>     refactoring 03-2011 [MI]
! **************************************************************************************************

   SUBROUTINE hfx_ks_matrix_gamma(qs_env, matrix_ks, rho, energy, calculate_forces, &
                                  just_energy, v_rspace_new, v_tau_rspace)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_new, v_tau_rspace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ks_matrix_gamma', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, img, irep, ispin, mspin, &
                                                            n_rep_hf, nimages, ns, nspins
      LOGICAL                                            :: distribute_fock_matrix, &
                                                            do_adiabatic_rescaling, &
                                                            hfx_treat_lsd_in_core, &
                                                            s_mstruct_changed, use_virial
      REAL(dp)                                           :: eh1, ehfx, ehfxrt, eold
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: hf_energy
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks_1d, matrix_ks_aux_fit, &
         matrix_ks_aux_fit_hfx, matrix_ks_aux_fit_im, matrix_ks_im, rho_ao_1d, rho_ao_resp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_ks_orb, rho_ao_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho_orb
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: adiabatic_rescaling_section, &
                                                            hfx_sections, input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, dft_control, force, hfx_sections, input, &
               para_env, poisson_env, pw_env, virial, matrix_ks_im, &
               matrix_ks_orb, rho_ao_orb, matrix_h, matrix_ks_aux_fit, &
               matrix_ks_aux_fit_im, matrix_ks_aux_fit_hfx)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      virial=virial, &
                      matrix_ks_im=matrix_ks_im, &
                      matrix_ks_aux_fit=matrix_ks_aux_fit, &
                      matrix_ks_aux_fit_im=matrix_ks_aux_fit_im, &
                      matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env=qs_env, mos_aux_fit=mo_array)
      ELSE
         CALL get_qs_env(qs_env=qs_env, mos=mo_array)
      ENDIF

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)
      adiabatic_rescaling_section => section_vals_get_subs_vals(input, "DFT%XC%ADIABATIC_RESCALING")
      CALL section_vals_get(adiabatic_rescaling_section, explicit=do_adiabatic_rescaling)

      ! *** Initialize the auxiliary ks matrix to zero if required
      IF (dft_control%do_admm) THEN
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
         END DO
      END IF
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_set(matrix_ks(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

      IF (calculate_forces) THEN
         !! initialize force array to zero
         CALL get_qs_env(qs_env=qs_env, force=force)
         DO ikind = 1, SIZE(force)
            force(ikind)%fock_4c(:, :) = 0.0_dp
         END DO
      END IF
      ALLOCATE (hf_energy(n_rep_hf))

      eold = 0.0_dp

      DO irep = 1, n_rep_hf
         ! Remember: Vhfx is added, energy is calclulated from total Vhfx,
         ! so energy of last iteration is correct

         IF (do_adiabatic_rescaling .AND. hfx_treat_lsd_in_core) &
            CPABORT("HFX_TREAT_LSD_IN_CORE not implemented for adiabatically rescaled hybrids")
         ! everything is calculated with adiabatic rescaling but the potential is not added in a first step
         distribute_fock_matrix = .NOT. do_adiabatic_rescaling

         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins

         ! fetch the correct matrices for normal HFX or ADMM
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env=qs_env, matrix_ks_aux_fit=matrix_ks_1d, &
                            rho_aux_fit=rho_orb)
            ns = SIZE(matrix_ks_1d)
            matrix_ks_orb(1:ns, 1:1) => matrix_ks_1d(1:ns)
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks_orb, rho=rho_orb)
         END IF
         CALL qs_rho_get(rho_struct=rho_orb, rho_ao_kp=rho_ao_orb)
         ! Finally the real hfx calulation
         ehfx = 0.0_dp

         IF (x_data(irep, 1)%do_hfx_ri .AND. calculate_forces) THEN
            CPABORT("RI forces not yet implemented in HFX")
         ENDIF

         IF (x_data(irep, 1)%do_hfx_ri) THEN
            CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_orb, ehfx, &
                                  mo_array, rho_ao_orb, &
                                  s_mstruct_changed, nspins, &
                                  x_data(irep, 1)%general_parameter%fraction)
            IF (dft_control%do_admm) THEN
               !for ADMMS, we need the exchange matrix k(d) for both spins
               DO ispin = 1, mspin
                  CALL dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_orb(ispin, 1)%matrix, &
                                  name="HF exch. part of matrix_ks_aux_fit for ADMMS")
               ENDDO
            END IF
         ELSE

            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_orb, eh1, rho_ao_orb, hfx_sections, &
                                          para_env, s_mstruct_changed, irep, distribute_fock_matrix, &
                                          ispin=ispin)
               ehfx = ehfx + eh1
            END DO
         ENDIF

         IF (calculate_forces .AND. .NOT. do_adiabatic_rescaling) THEN
            !Scale auxiliary density matrix for ADMMP (see Merlot2014) with gsi(ispin) to scale force
            IF (dft_control%do_admm) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.FALSE.)
            END IF
            IF (ASSOCIATED(qs_env%mp2_env)) THEN
               CALL get_qs_env(qs_env, matrix_p_mp2=rho_ao_resp)
            ELSE
               NULLIFY (rho_ao_resp)
            END IF
            CALL derivatives_four_center(qs_env, rho_ao_orb, rho_ao_resp, hfx_sections, &
                                         para_env, irep, use_virial)
            !Scale auxiliary density matrix for ADMMP back with 1/gsi(ispin)
            IF (dft_control%do_admm) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.TRUE.)
            END IF
         END IF

         !! If required, the calculation of the forces will be done later with adiabatic rescaling
         IF (do_adiabatic_rescaling) hf_energy(irep) = ehfx

         ! special case RTP/EMD we have a full complex density and HFX has a contribution from the imaginary part
         ehfxrt = 0.0_dp
         IF (qs_env%run_rtp) THEN

            CALL get_qs_env(qs_env=qs_env, rtp=rtp)
            DO ispin = 1, nspins
               CALL dbcsr_set(matrix_ks_im(ispin)%matrix, 0.0_dp)
            END DO
            IF (dft_control%do_admm) THEN
               ! matrix_ks_orb => matrix_ks_aux_fit_im
               ns = SIZE(matrix_ks_aux_fit_im)
               matrix_ks_orb(1:ns, 1:1) => matrix_ks_aux_fit_im(1:ns)
               DO ispin = 1, nspins
                  CALL dbcsr_set(matrix_ks_aux_fit_im(ispin)%matrix, 0.0_dp)
               END DO
            ELSE
               ! matrix_ks_orb => matrix_ks_im
               ns = SIZE(matrix_ks_im)
               matrix_ks_orb(1:ns, 1:1) => matrix_ks_im(1:ns)
            END IF

            CALL qs_rho_get(rho_orb, rho_ao_im=rho_ao_1d)
            ns = SIZE(rho_ao_1d)
            rho_ao_orb(1:ns, 1:1) => rho_ao_1d(1:ns)

            ehfxrt = 0.0_dp

            IF (x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_orb, ehfx, &
                                     mo_array, rho_ao_orb, &
                                     .FALSE., nspins, &
                                     x_data(irep, 1)%general_parameter%fraction)
               IF (dft_control%do_admm) THEN
                  !for ADMMS, we need the exchange matrix k(d) for both spins
                  DO ispin = 1, mspin
                     CALL dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_orb(ispin, 1)%matrix, &
                                     name="HF exch. part of matrix_ks_aux_fit for ADMMS")
                  ENDDO
               END IF

            ELSE
               DO ispin = 1, mspin
                  CALL integrate_four_center(qs_env, x_data, matrix_ks_orb, eh1, rho_ao_orb, hfx_sections, &
                                             para_env, .FALSE., irep, distribute_fock_matrix, &
                                             ispin=ispin)
                  ehfxrt = ehfxrt + eh1
               END DO
            ENDIF

            IF (calculate_forces .AND. .NOT. do_adiabatic_rescaling) THEN
               NULLIFY (rho_ao_resp)
               CALL derivatives_four_center(qs_env, rho_ao_orb, rho_ao_resp, hfx_sections, &
                                            para_env, irep, use_virial)
            END IF

            !! If required, the calculation of the forces will be done later with adiabatic rescaling
            IF (do_adiabatic_rescaling) hf_energy(irep) = ehfx + ehfxrt
         END IF

         IF (.NOT. qs_env%run_rtp) THEN
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                            poisson_env=poisson_env)
            eh1 = ehfx - eold
            CALL pw_hfx(qs_env, eh1, hfx_sections, poisson_env, auxbas_pw_pool, irep)
            eold = ehfx
         END IF

      END DO

      ! *** Set the total HFX energy
      energy%ex = ehfx + ehfxrt

      ! *** Add Core-Hamiltonian-Matrix ***
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_add(matrix_ks(ispin, img)%matrix, matrix_h(1, img)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
      END DO
      IF (use_virial .AND. calculate_forces) THEN
         virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         virial%pv_calculate = .FALSE.
      ENDIF

      !! If we perform adiabatic rescaling we are now able to rescale the xc-potential
      IF (do_adiabatic_rescaling) THEN
         CALL rescale_xc_potential(qs_env, matrix_ks, rho, energy, v_rspace_new, v_tau_rspace, &
                                   hf_energy, just_energy, calculate_forces, use_virial)
      END IF ! do_adiabatic_rescaling

      !update the hfx aux_fit matrixIF (dft_control%do_admm) THEN
      IF (dft_control%do_admm) THEN
         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, &
                           0.0_dp, 1.0_dp)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_ks_matrix_gamma

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian (kpoint version)
!>
!> \param qs_env ...
!> \param matrix_ks ...
!> \param energy ...
!> \param calculate_forces ...
!> \par History
!>     refactoring 06-2020 [JGH]
! **************************************************************************************************

   SUBROUTINE hfx_ks_matrix_kp(qs_env, matrix_ks, energy, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ks_matrix_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, img, irep, ispin, mspin, &
                                                            n_rep_hf, nimages, ns, nspins
      LOGICAL                                            :: do_adiabatic_rescaling, &
                                                            hfx_treat_lsd_in_core, distribute_fock_matrix, &
                                                            s_mstruct_changed, use_virial
      REAL(dp)                                           :: eh1, ehfx, eold
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: hf_energy
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_1d, matrix_ks_aux_fit, &
                                                            matrix_ks_aux_fit_hfx, &
                                                            matrix_ks_aux_fit_im, matrix_ks_im, &
                                                            rho_ao_resp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_ks_orb, rho_ao_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho_orb
      TYPE(section_vals_type), POINTER                   :: adiabatic_rescaling_section, &
                                                            hfx_sections, input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, dft_control, force, hfx_sections, input, &
               para_env, poisson_env, pw_env, virial, matrix_ks_im, &
               matrix_ks_orb, rho_ao_orb, matrix_h, matrix_ks_aux_fit, &
               matrix_ks_aux_fit_im, matrix_ks_aux_fit_hfx)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      virial=virial, &
                      matrix_ks_im=matrix_ks_im, &
                      matrix_ks_aux_fit=matrix_ks_aux_fit, &
                      matrix_ks_aux_fit_im=matrix_ks_aux_fit_im, &
                      matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      ! No RTP
      CPASSERT(.NOT. qs_env%run_rtp)

      ! MOs are in kpoint_env ?
      ! get them here redistributed?
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env=qs_env, mos_aux_fit=mo_array)
      ELSE
         CALL get_qs_env(qs_env=qs_env, mos=mo_array)
      ENDIF

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)
      ! no adiabatic rescaling with kpoints
      adiabatic_rescaling_section => section_vals_get_subs_vals(input, "DFT%XC%ADIABATIC_RESCALING")
      CALL section_vals_get(adiabatic_rescaling_section, explicit=do_adiabatic_rescaling)
      CPASSERT(.NOT. do_adiabatic_rescaling)

      ! probably correct?
      distribute_fock_matrix = .TRUE.

      ! *** Initialize the auxiliary ks matrix to zero if required
      IF (dft_control%do_admm) THEN
         CPABORT("ADMM KP")
         DO ispin = 1, nspins
            ! admm matrices need image index
            CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
         END DO
      END IF
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_set(matrix_ks(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

      IF (calculate_forces) THEN
         !! initialize force array to zero
         CALL get_qs_env(qs_env=qs_env, force=force)
         DO ikind = 1, SIZE(force)
            force(ikind)%fock_4c(:, :) = 0.0_dp
         END DO
      END IF
      ALLOCATE (hf_energy(n_rep_hf))

      eold = 0.0_dp

      DO irep = 1, n_rep_hf
         ! Remember: Vhfx is added, energy is calclulated from total Vhfx,
         ! so energy of last iteration is correct

         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins

         ! fetch the correct matrices for normal HFX or ADMM
         IF (dft_control%do_admm) THEN
            CPABORT("ADMM KP")
            ! ADMM matrices need additional index
            CALL get_qs_env(qs_env=qs_env, matrix_ks_aux_fit=matrix_ks_1d, &
                            rho_aux_fit=rho_orb)
            ns = SIZE(matrix_ks_1d)
            matrix_ks_orb(1:ns, 1:1) => matrix_ks_1d(1:ns)
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks_orb, rho=rho_orb)
         END IF
         CALL qs_rho_get(rho_struct=rho_orb, rho_ao_kp=rho_ao_orb)
         ! Finally the real hfx calulation
         ehfx = 0.0_dp

         IF (x_data(irep, 1)%do_hfx_ri .AND. calculate_forces) THEN
            CPABORT("RI HFX for KP not yet implemented")
            CPABORT("RI forces not yet implemented in HFX")
         ENDIF

         IF (x_data(irep, 1)%do_hfx_ri) THEN
            CPABORT("RI HFX for KP not yet implemented")
            CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_orb, ehfx, &
                                  mo_array, rho_ao_orb, &
                                  s_mstruct_changed, nspins, &
                                  x_data(irep, 1)%general_parameter%fraction)
            IF (dft_control%do_admm) THEN
               CPABORT("ADMM KP")
               !for ADMMS, we need the exchange matrix k(d) for both spins
               DO ispin = 1, mspin
                  CALL dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_orb(ispin, 1)%matrix, &
                                  name="HF exch. part of matrix_ks_aux_fit for ADMMS")
               ENDDO
            END IF
         ELSE
            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_orb, eh1, rho_ao_orb, hfx_sections, &
                                          para_env, s_mstruct_changed, irep, distribute_fock_matrix, &
                                          ispin=ispin)
               ehfx = ehfx + eh1
            END DO
         ENDIF

         IF (calculate_forces) THEN
            !Scale auxiliary density matrix for ADMMP (see Merlot2014) with gsi(ispin) to scale force
            IF (dft_control%do_admm) THEN
               CPABORT("ADMM KP")
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.FALSE.)
            END IF
            IF (ASSOCIATED(qs_env%mp2_env)) THEN
               CPABORT("MP2 KP")
               CALL get_qs_env(qs_env, matrix_p_mp2=rho_ao_resp)
            ELSE
               NULLIFY (rho_ao_resp)
            END IF
!deb        CALL derivatives_four_center(qs_env, rho_ao_orb, rho_ao_resp, hfx_sections, &
!deb                                     para_env, irep, use_virial)
            !Scale auxiliary density matrix for ADMMP back with 1/gsi(ispin)
            IF (dft_control%do_admm) THEN
               CPABORT("ADMM KP")
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.TRUE.)
            END IF
         END IF

         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
         eh1 = ehfx - eold
         CALL pw_hfx_kp(qs_env, eh1, hfx_sections, poisson_env, auxbas_pw_pool, irep)
         eold = ehfx

      END DO

      ! *** Set the total HFX energy
      energy%ex = ehfx

      ! *** Add Core-Hamiltonian-Matrix ***
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_add(matrix_ks(ispin, img)%matrix, matrix_h(1, img)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
      END DO
      IF (use_virial .AND. calculate_forces) THEN
         virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         virial%pv_calculate = .FALSE.
      ENDIF

      !update the hfx aux_fit matrixIF (dft_control%do_admm) THEN
      IF (dft_control%do_admm) THEN
         CPABORT("ADMM KP")
         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, &
                           0.0_dp, 1.0_dp)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_ks_matrix_kp

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param matrix_ks Kohn-Sham matrix (updated on exit)
!> \param rho_ao    electron density expressed in terms of atomic orbitals
!> \param qs_env    Quickstep environment
!> \note
!>     Simplified version of subroutine hfx_ks_matrix()
! **************************************************************************************************
   SUBROUTINE tddft_hfx_matrix(matrix_ks, rho_ao, qs_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, rho_ao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_hfx_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, irep, ispin, mspin, n_rep_hf, &
                                                            nspins
      LOGICAL                                            :: distribute_fock_matrix, do_hfx, &
                                                            hfx_treat_lsd_in_core, &
                                                            s_mstruct_changed
      REAL(KIND=dp)                                      :: eh1, ehfx
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, hfx_sections, input, para_env, matrix_ks_kp, rho_ao_kp)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      para_env=para_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      IF (do_hfx) THEN
         CPASSERT(dft_control%nimages == 1)
         nspins = dft_control%nspins

         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
         CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)

         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
         distribute_fock_matrix = .TRUE.

         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins

         matrix_ks_kp(1:nspins, 1:1) => matrix_ks(1:nspins)
         rho_ao_kp(1:nspins, 1:1) => rho_ao(1:nspins)

         DO irep = 1, n_rep_hf
            ! the real hfx calulation
            ehfx = 0.0_dp
            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_kp, eh1, rho_ao_kp, hfx_sections, para_env, &
                                          s_mstruct_changed, irep, distribute_fock_matrix, ispin=ispin)
               ehfx = ehfx + eh1
            END DO
         END DO
         energy%ex = ehfx
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddft_hfx_matrix

END MODULE hfx_ks_matrix
