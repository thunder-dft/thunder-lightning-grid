! copyright info:
!
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Program Description
! ===========================================================================
!> This is the main driver for the LIGHTNING (grid) version of FIREBALL.
! ==========================================================================
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program lightning

! /GLOBAL
        use M_welcome

! /SYSTEM
        use M_species
        use M_configuraciones
        use M_neighbors
        use M_neighbors_PP
        use M_atom_functions

! /FDATA
        use M_Fdata_1c
        use M_Fdata_2c

! /GRID
        use M_grid
        use M_project_grid
        use M_fft_grid

! /ASSEMBLERS
        use M_assemble_2c
        use M_assemble_3c
        use M_assemble_ewald
        use M_assemble_usr
        use M_assemble_vxc
        use M_assemble_PP_2c
        use M_assemble_PP_3c

! /SOLVESH
        use M_kspace
        use M_density_matrix

! /SCF
        use M_charges

        implicit none

        interface
           subroutine Qmixer (t, iscf_iteration, sigma)
             use M_configuraciones
             implicit none
             integer, intent (in) :: iscf_iteration
             type(T_structure), target :: t
             real, intent (inout) :: sigma
           end subroutine Qmixer

           subroutine absorption (t)
             use M_species
             use M_configuraciones
             use M_atom_functions
             implicit none
             type(T_structure), target :: t           !< the structure to be used
           end subroutine absorption

           subroutine dos (t)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t           !< the structure to be used
           end subroutine dos

           subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
             use M_assemble_blocks
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t             ! the structure to be used
             real, intent (in) :: ebs                   ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc       ! short-range energies
           end subroutine writeout_energies

           subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t             ! the structure to be used
             real, intent (in) :: ebs                   ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc       ! short-range energies
           end subroutine writeout_xyz
        end interface

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iscf_iteration
        integer istructure, iseparate

        real sigma

        character (len = 25) :: slogfile

!       real rcbohr                      !< cutoff in Bohr radii

!       character (len = 11) buffer     !< buffer for generating wavefunction
!       character (len = 3) rcchar
!       character (len = 1), dimension (0:9) :: zchar

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        real timei, timef

! Energies
        real ebs                                     ! band-structure energy
        real uii_uee, uxcdcc                         ! short-range energies
!       real etot                                    ! total energy
!       real etot_per_atom                           ! total energy per atom

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        iseparate = 1
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome_lightning

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
        call read_Fdata_location
        call read_info
        call read_Fdata_1c
        call read_Fdata_2c
        call read_Fdata_3c

! Read in the wavefunctions and potentials
        call read_create
        call read_wavefunctions
        call read_napotentials

! Read parameters from structures.inp file
        write (ilogfile,*)
        write (ilogfile,*) ' Reading parameters from the structure input. '
        call read_parameters

        write (ilogfile, *)
        open (unit = 1, file = 'structures.inp', status = 'old')
        read (1, *) nstructures
        if (nstructures .gt. 999) then
          stop ' Cannot calculate more than 999 structures! '
        end if
        write (ilogfile, *) ' Number of structure calculating = ', nstructures
        allocate (structures (nstructures))

! Loop over all structures
! This loop can be made parallel if each subroutine in lightning
! is modified to take s as a parameter, rather than reference s directly
        do istructure = 1, nstructures
          s => structures(istructure)
          read (1, *) s%basisfile
          if (iseparate .eq. 1) then
            s%logfile = istructure + 1000
            s%inpfile = istructure + 2000
            slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
            slogfile = trim(slogfile)//'.log'
            open (unit = s%logfile, file = slogfile, status = 'replace')
          end if
          write (ilogfile,*)
          write (ilogfile, 100) slogfile
          write (s%logfile, *) ' Structure = ', istructure

          ! Read in the coordinates and parameters
          call read_positions (s)

          ! Set the charges
          call read_charges (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!           N E I G H B O R S   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
          call driver_neighbors (s)
          call driver_neighbors_PP (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!              A S S E M B L E   T H E   H A M I L T O N I A N
! ---------------------------------------------------------------------------
! ===========================================================================
! Assemble the Hamiltonian matrix:
          write (s%logfile, *)
          write (s%logfile, *) ' Calling two-center non-charge dependent assemblers. '
          call assemble_S (s)
          call assemble_T (s)
          call assemble_svnl (s)
          call assemble_vnl_2c (s)

          write (s%logfile,*) ' Calling three-center non-charge dependent assemblers. '
          call assemble_vnl_3c (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                  G R I D    I N I T I A L I Z A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
          write (s%logfile, *)
          write (s%logfile, *) ' Initialize the grid. '
          call initialize_grid (s)
          call initialize_project_grid (s)
          call initialize_density_matrix (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!              S C F   L O O P
! ---------------------------------------------------------------------------
! Note that self-consistency is preformed regardless of method used.
! But, in Harris, we just do a single scf cycle.
! ===========================================================================
          sigma = 999.0d0
          iscf_iteration = 1
          do while (sigma .gt. scf_tolerance_set .and.                       &
      &             iscf_iteration .le. max_scf_iterations_set)
            write (s%logfile, *)
            write (s%logfile, *) ' Begin scf step = ', iscf_iteration
            write (s%logfile, *) ' Calling two-center charge dependent assemblers. '
            call assemble_vna_2c (s)

            write (s%logfile, *) ' Calling three-center charge dependent assemblers. '
            call assemble_vna_3c (s)
            call assemble_vxc (s)

            write (s%logfile, *) ' Calling Kohn-Sham assemblers. '
            ! integrate vxc[rhoG] and vna[drhoG] on the grid
            call assemble_Gmatrix (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
            write (s%logfile, *) ' Calling kspace: '
            call cpu_time(timei)
            call driver_kspace (s, iscf_iteration)
            call cpu_time(timef)
            write (s%logfile, *)
            write (s%logfile, *) ' kspace time: ', timef - timei
            call density_matrix (s)
            if (iwriteout_density .eq. 1) call writeout_density (s)

            call calculate_charges (s)
            call Qmixer (s, iscf_iteration, sigma)
            if (iwriteout_charges .eq. 1) call writeout_charges (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                       T O T A L   E N E R G I E S
! ---------------------------------------------------------------------------
! ===========================================================================
            call calculate_ebs (s, ebs)
            call assemble_uee (s, uii_uee)
            call assemble_uxc (s, uxcdcc)
            call writeout_energies (s, ebs, uii_uee, uxcdcc)

            ! Evaluate new density on grid
            call project_density_grid (s)
            ! Solve Poisson equation to get new Hartree potential
            call Laplace_fft (s)

! End scf loop
            if (sigma .gt. 0.0d0) then
              iscf_iteration = iscf_iteration + 1
            else
              exit
            end if
          end do

! Write out the xyz file - we need this for examining the .xsf files
          if (iwriteout_xyz .eq. 1) call writeout_xyz (s, ebs, uii_uee, uxcdcc)

! Destroy Hamiltonian matrix elements storage
          call destroy_assemble_2c (s)
          call destroy_assemble_PP_2c (s)
          call destroy_assemble_vxc_McWEDA (s)

          call destroy_grid

! Calculate the absorption spectrum.
          if (iwriteout_abs .eq. 1) call absorption (s)

! Destroy final arrays
          call destroy_kspace (s)
          call destroy_denmat (s)
          call destroy_charges (s)
          call destroy_neighbors (s)
          call destroy_neighbors_PP (s)

          ! destroy neighbors last
          deallocate (s%xl) ! where to put this?

! Calculate the electronic density of states.
! We do this after destorying some arrays so that we can optimize the
! memory usage.
          if (iwriteout_dos .eq. 1) call dos (s)

          if (iseparate .eq. 1) close (s%logfile)
        end do
        close (1)
        ! end loop over all structures

! ===========================================================================
! ---------------------------------------------------------------------------
!               D E S T R O Y   A R R A Y S - F I N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Destroy datafile storage
        call destroy_Fdata_1C
        call destroy_Fdata_2C
        call destroy_Fdata_3c

! Destroy SYSTEM information.
        call destroy_species

        call cpu_time (time_end)

        write (ilogfile,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        write (*,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Working on structure - ', a25)

502     format (2x, '           ebs = ', f15.6)
503     format (2x, '     uii - uee = ', f15.6)
505     format (2x, '        uxcdcc = ', f15.6)
507     format (2x, '          ETOT = ', f15.6)
508     format (2x, '     Etot/atom = ', f15.6)
509     format (2x, ' Atomic Energy = ', f15.6)
510     format (2x, '     CohesiveE = ', f15.6)
511     format (2x, ' Cohesive Energy per atom  = ', f15.6)

! End Program
! ===========================================================================
        stop
        end program lightning
