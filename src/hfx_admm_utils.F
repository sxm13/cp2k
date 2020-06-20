!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for hfx and admm methods
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!>     Made GAPW compatible 12.2019 (A. Bussy)
!> \author MI
! **************************************************************************************************
MODULE hfx_admm_utils
   USE admm_dm_types,                   ONLY: admm_dm_create,&
                                              admm_dm_type
   USE admm_types,                      ONLY: admm_env_create,&
                                              admm_gapw_type,&
                                              admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE external_potential_types,        ONLY: copy_potential
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_bee, do_admm_aux_exch_func_default, do_admm_aux_exch_func_none, &
        do_admm_aux_exch_func_opt, do_admm_aux_exch_func_pbex, do_potential_coulomb, &
        do_potential_long, do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_short, &
        do_potential_truncated, xc_funct_no_shortcut
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              init_gapw_basis_set,&
                                              init_gapw_nlcc,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_local_rho_types,              ONLY: local_rho_set_create
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_rho_atom_methods,             ONLY: init_rho_atom
   USE qs_rho_types,                    ONLY: qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: hfx_admm_init, create_admm_xc_section

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_admm_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param do_mp2_hfx ...
! **************************************************************************************************
   SUBROUTINE hfx_admm_init(qs_env, do_mp2_hfx)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in), OPTIONAL                      :: do_mp2_hfx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_admm_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_rep_hf, natoms, nspins
      LOGICAL                                            :: do_hfx, my_do_mp2_hfx, s_mstruct_changed
      TYPE(admm_dm_type), POINTER                        :: admm_dm
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_aux_fit
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input, xc_section

      my_do_mp2_hfx = .FALSE.
      IF (PRESENT(do_mp2_hfx)) my_do_mp2_hfx = do_mp2_hfx

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, admm_dm, hfx_sections, mos, mos_aux_fit, para_env, &
               particle_set, xc_section, matrix_s_aux_fit, &
               matrix_s_aux_fit_vs_orb, rho, rho_aux_fit, ks_env, dft_control, &
               input)

      CALL get_qs_env(qs_env, &
                      mos_aux_fit=mos_aux_fit, &
                      mos=mos, &
                      admm_env=admm_env, &
                      admm_dm=admm_dm, &
                      matrix_s_aux_fit=matrix_s_aux_fit, &
                      matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, &
                      para_env=para_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      rho=rho, &
                      rho_aux_fit=rho_aux_fit, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      input=input)

      nspins = dft_control%nspins

      IF (my_do_mp2_hfx) THEN
         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      ELSE
         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      END IF

      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      !! ** ADMM can only be used with HFX
      IF (.NOT. do_hfx) &
         CPABORT("Wavefunction fitting requested without Hartree-Fock.")

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      IF (n_rep_hf > 1) &
         CPABORT("ADMM can handle only one HF section.")

      IF (.NOT. ASSOCIATED(admm_env)) THEN
         ! setup admm environment
         CALL get_qs_env(qs_env, input=input, particle_set=particle_set)
         natoms = SIZE(particle_set, 1)
         CALL admm_env_create(admm_env, dft_control%admm_control, mos, mos_aux_fit, &
                              para_env, natoms)
         CALL set_qs_env(qs_env, admm_env=admm_env)
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         CALL create_admm_xc_section(qs_env=qs_env, xc_section=xc_section, &
                                     admm_env=admm_env)

         ! Initialize the GAPW data types
         IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) &
            CALL init_admm_gapw(qs_env)

      ENDIF

      IF (admm_env%do_gapw .AND. dft_control%do_admm_dm) THEN
         CPABORT("GAPW ADMM not implemented for MCWEENY or NONE_DM purification.")
      END IF

      IF (admm_env%do_gapw .AND. qs_env%run_rtp) THEN
         CPABORT("GAPW ADMM not implemented for real-time propagation.")
      END IF

      IF (dft_control%do_admm_dm .AND. .NOT. ASSOCIATED(admm_dm)) THEN
         CALL admm_dm_create(admm_dm, dft_control%admm_control, nspins=nspins, natoms=natoms)
         CALL set_ks_env(ks_env, admm_dm=admm_dm)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE hfx_admm_init

! **************************************************************************************************
!> \brief Sets up the admm_gapw env
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE init_admm_gapw(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_admm_gapw', routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nkind
      TYPE(admm_gapw_type), POINTER                      :: admm_gapw_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis, aux_fit_soft_basis, &
                                                            orb_basis, soft_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: admm_kind_set, qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      NULLIFY (admm_kind_set, aux_fit_basis, atomic_kind_set, aux_fit_soft_basis, &
               dft_control, input, orb_basis, para_env, qs_kind_set, soft_basis)

      CALL get_qs_env(qs_env, admm_env=admm_env, &
                      atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, &
                      input=input, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set)

      admm_env%do_gapw = .TRUE.
      ALLOCATE (admm_env%admm_gapw_env)
      admm_gapw_env => admm_env%admm_gapw_env
      NULLIFY (admm_gapw_env%local_rho_set)
      NULLIFY (admm_gapw_env%admm_kind_set)
      NULLIFY (admm_gapw_env%task_list)

      !Create a new kind set for the ADMM stuff (paw_proj soft AUX_FIT basis, etc)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (admm_gapw_env%admm_kind_set(nkind))
      admm_kind_set => admm_gapw_env%admm_kind_set

      !In this new kind set, we want the AUX_FIT basis to be known as ORB, such that GAPW routines work
      DO ikind = 1, nkind
         !copying over simple data  of interest from qs_kind_set
         admm_kind_set(ikind)%name = qs_kind_set(ikind)%name
         admm_kind_set(ikind)%element_symbol = qs_kind_set(ikind)%element_symbol
         admm_kind_set(ikind)%natom = qs_kind_set(ikind)%natom
         admm_kind_set(ikind)%hard_radius = qs_kind_set(ikind)%hard_radius
         admm_kind_set(ikind)%max_rad_local = qs_kind_set(ikind)%max_rad_local
         admm_kind_set(ikind)%gpw_type_forced = qs_kind_set(ikind)%gpw_type_forced
         admm_kind_set(ikind)%ngrid_rad = qs_kind_set(ikind)%ngrid_rad
         admm_kind_set(ikind)%ngrid_ang = qs_kind_set(ikind)%ngrid_ang

         !copying potentials of interest from qs_kind_set
         IF (ASSOCIATED(qs_kind_set(ikind)%all_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%all_potential, admm_kind_set(ikind)%all_potential)
         END IF
         IF (ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%gth_potential, admm_kind_set(ikind)%gth_potential)
         END IF
         IF (ASSOCIATED(qs_kind_set(ikind)%sgp_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%sgp_potential, admm_kind_set(ikind)%sgp_potential)
         END IF

         NULLIFY (orb_basis)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=aux_fit_basis, basis_type="AUX_FIT")
         CALL copy_gto_basis_set(aux_fit_basis, orb_basis)
         CALL add_basis_set_to_container(admm_kind_set(ikind)%basis_sets, orb_basis, "ORB")
      END DO

      !Create the corresponding soft basis set (and projectors)
      CALL init_gapw_basis_set(admm_kind_set, dft_control%qs_control, input, &
                               modify_qs_control=.FALSE.)

      !Make sure the basis and the projectors are well initialized
      CALL init_interaction_radii(dft_control%qs_control, admm_kind_set)

      !We also init the atomic grids and harmonics
      CALL local_rho_set_create(admm_gapw_env%local_rho_set)
      CALL init_rho_atom(admm_gapw_env%local_rho_set%rho_atom_set, &
                         atomic_kind_set, admm_kind_set, dft_control, para_env)

      !Make sure that any NLCC potential is well initialized
      CALL init_gapw_nlcc(admm_kind_set)

      !Need to have access to the soft AUX_FIT basis from the qs_env => add it to the qs_kinds
      DO ikind = 1, nkind
         NULLIFY (aux_fit_soft_basis)
         CALL get_qs_kind(admm_kind_set(ikind), soft_basis_set=soft_basis)
         CALL copy_gto_basis_set(soft_basis, aux_fit_soft_basis)
         CALL add_basis_set_to_container(qs_kind_set(ikind)%basis_sets, aux_fit_soft_basis, "AUX_FIT_SOFT")
      END DO

   END SUBROUTINE init_admm_gapw

! **************************************************************************************************
!> \brief This routine modifies the xc section depending on the potential type
!>        used for the HF exchange and the resulting correction term. Currently
!>        three types of corrections are implemented:
!>
!>        coulomb:     Ex,hf = Ex,hf' + (PBEx-PBEx')
!>        shortrange:  Ex,hf = Ex,hf' + (XWPBEX-XWPBEX')
!>        truncated:   Ex,hf = Ex,hf' + ( (XWPBEX0-PBE_HOLE_TC_LR) -(XWPBEX0-PBE_HOLE_TC_LR)' )
!>
!>        with ' denoting the auxiliary basis set and
!>
!>          PBEx:           PBE exchange functional
!>          XWPBEX:         PBE exchange hole for short-range potential (erfc(omega*r)/r)
!>          XWPBEX0:        PBE exchange hole for standard coulomb potential
!>          PBE_HOLE_TC_LR: PBE exchange hole for longrange truncated coulomb potential
!>
!>
!> \param qs_env the qs environment
!> \param xc_section the original xc_section
!> \param admm_env the ADMM environment
!> \par History
!>      12.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE create_admm_xc_section(qs_env, xc_section, admm_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(admm_type), POINTER                           :: admm_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_admm_xc_section', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: debug_functional = .FALSE.

      CHARACTER(LEN=20)                                  :: name_x_func
      INTEGER                                            :: hfx_potential_type, ifun, iounit, nfun
      LOGICAL                                            :: funct_found
      REAL(dp)                                           :: cutoff_radius, hfx_fraction, omega, &
                                                            scale_coulomb, scale_longrange, scale_x
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: xc_fun, xc_fun_section

      logger => cp_get_default_logger()
      NULLIFY (admm_env%xc_section_aux, admm_env%xc_section_primary)

      CALL get_qs_env(qs_env)

      !! ** Duplicate existing xc-section
      CALL section_vals_duplicate(xc_section, admm_env%xc_section_aux)
      CALL section_vals_duplicate(xc_section, admm_env%xc_section_primary)
      !** Now modify the auxiliary basis
      !** First remove all functionals
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux, "XC_FUNCTIONAL")

      !* Overwrite possible shortcut
      CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                i_val=xc_funct_no_shortcut)

      !** Get number of Functionals in the list
      ifun = 0
      nfun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         nfun = nfun + 1
      END DO

      ifun = 0
      DO ifun = 1, nfun
         xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=1)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         CALL section_vals_remove_values(xc_fun)
      END DO

      hfx_potential_type = qs_env%x_data(1, 1)%potential_parameter%potential_type
      hfx_fraction = qs_env%x_data(1, 1)%general_parameter%fraction

      !in case of no admm exchange corr., no auxiliary exchange functional needed
      IF (admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
         hfx_fraction = 0.0_dp
      END IF

      ! default PBE Functional
      IF (admm_env%aux_exch_func == do_admm_aux_exch_func_default .OR. &
          admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN

         !! ** Add functionals evaluated with auxiliary basis
         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            CALL section_vals_val_set(xc_fun_section, "PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_C", &
                                      r_val=0.0_dp)
         CASE (do_potential_short)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=0.0_dp)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_truncated)
            cutoff_radius = qs_env%x_data(1, 1)%potential_parameter%cutoff_radius
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=0.0_dp)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction)
         CASE (do_potential_long)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_mix_cl)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            scale_coulomb = qs_env%x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = qs_env%x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_mix_cl_trunc)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            cutoff_radius = qs_env%x_data(1, 1)%potential_parameter%cutoff_radius
            scale_coulomb = qs_env%x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = qs_env%x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE DEFAULT
            CPABORT("")
         END SELECT

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_C", &
                                         r_val=0.0_dp)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_short)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=0.0_dp)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_long)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF
         CASE (do_potential_truncated)
            cutoff_radius = qs_env%x_data(1, 1)%potential_parameter%cutoff_radius
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=0.0_dp)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_mix_cl_trunc)
            cutoff_radius = qs_env%x_data(1, 1)%potential_parameter%cutoff_radius
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            scale_coulomb = qs_env%x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = qs_env%x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)

               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF
         CASE (do_potential_mix_cl)
            omega = qs_env%x_data(1, 1)%potential_parameter%omega
            scale_coulomb = qs_env%x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = qs_env%x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)

               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF

         END SELECT

         ! PBEX (always bare form), OPTX and Becke88 functional
      ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_opt .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_bee) THEN
         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
            name_x_func = 'PBE'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
            name_x_func = 'OPTX'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_bee) THEN
            name_x_func = 'BECKE88'
         END IF
         !primary basis
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                   l_val=.TRUE.)
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                   r_val=-hfx_fraction)

         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_C", r_val=0.0_dp)
         END IF

         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
            IF (admm_env%aux_exch_func_param) THEN
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A1", &
                                         r_val=admm_env%aux_x_param(1))
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A2", &
                                         r_val=admm_env%aux_x_param(2))
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%GAMMA", &
                                         r_val=admm_env%aux_x_param(3))
            END IF
         END IF

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible L")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            IF (xc_fun%section%name == TRIM(name_x_func)) THEN
               funct_found = .TRUE.
            END IF
         END DO
         IF (.NOT. funct_found) THEN
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=hfx_fraction)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_C", &
                                         r_val=0.0_dp)
            ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
               IF (admm_env%aux_exch_func_param) THEN
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A1", &
                                            r_val=admm_env%aux_x_param(1))
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A2", &
                                            r_val=admm_env%aux_x_param(2))
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%GAMMA", &
                                            r_val=admm_env%aux_x_param(3))
               END IF
            END IF

         ELSE
            CALL section_vals_val_get(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=scale_x)
            scale_x = scale_x + hfx_fraction
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=scale_x)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
               CPASSERT(.NOT. admm_env%aux_exch_func_param)
            END IF
         END IF

      END IF

      IF (debug_functional) THEN
         iounit = cp_logger_get_default_io_unit(logger)
         IF (iounit > 0) THEN
            WRITE (iounit, "(A)") " ADMM Primary Basis Set Functional"
         END IF
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT

            scale_x = -1000.0_dp
            IF (xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X", r_val=scale_x)
            END IF
            IF (xc_fun%section%name == "XWPBE") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X0", r_val=hfx_fraction)
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,2F10.3)") TRIM(xc_fun%section%name), scale_x, hfx_fraction
               END IF
            ELSE
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,F10.3)") TRIM(xc_fun%section%name), scale_x
               END IF
            END IF
         END DO

         IF (iounit > 0) THEN
            WRITE (iounit, "(A)") " Auxiliary Basis Set Functional"
         END IF
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux, "XC_FUNCTIONAL")
         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            scale_x = -1000.0_dp
            IF (xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X", r_val=scale_x)
            END IF
            IF (xc_fun%section%name == "XWPBE") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X0", r_val=hfx_fraction)
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,2F10.3)") TRIM(xc_fun%section%name), scale_x, hfx_fraction
               END IF
            ELSE
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,F10.3)") TRIM(xc_fun%section%name), scale_x
               END IF
            END IF
         END DO
      END IF

   END SUBROUTINE create_admm_xc_section

END MODULE hfx_admm_utils
