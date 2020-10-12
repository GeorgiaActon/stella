program stella_adjoint

  use redistribute, only: scatter
  use job_manage, only: time_message, checkstop
  use run_parameters, only: nstep, fphi, fapar
  use stella_time, only: update_time, code_time, code_dt
  use dist_redistribute, only: kxkyz2vmu
  use time_advance, only: advance_stella
  use stella_diagnostics, only: diagnose_stella, nsave
  use stella_save, only: stella_save_for_restart
  use dist_fn_arrays, only: gnew, gvmu, gold
  use dist_fn_arrays, only: lambda_new, lambda_old, lambda_vmu
  use dist_fn_arrays, only: g_k, g_tau, omega
  use adjoint_convergence_g, only: get_omega, check_omega_convergence, get_gk
  use adjoint_convergence_g, only: check_gk_convergence
  use kt_grids, only: nakx, naky
  use vpamu_grids, only: integrate_vmu
  use adjoint_p_derviatives, only: perturb_p
  use advancing_g_analytic, only: g_analytic
  use stella_time, only : code_time
  use advancing_g_analytic, only: get_g_tau, get_tau1
  use fields_arrays, only: phi_old, phi
  use lambda_subroutines, only: advance_adjoint
!  use adjoint_advance, only : advance_adjoint_s3
  use stella_diagnostics, only: omega_vs_time
  
  logical :: debug = .true.
  logical :: stop_stella = .false.
  logical :: mpi_initialized = .false.

  integer :: istep
  integer :: istatus
  real, dimension (2) :: time_init = 0.
  real, dimension (2) :: time_diagnostics = 0.
  real, dimension (2) :: time_total = 0.

   complex, dimension (:,:), allocatable :: sum_omega
!   complex, dimension (:,:,:,:,:), allocatable :: sum_gk
   complex, dimension(:,:,:,:), allocatable :: dG_dp
   complex, dimension(:,:), allocatable :: dQ_dp
   complex, dimension(:), allocatable :: sum_lag, sum_beta_q, lag
   real, dimension (:,:), allocatable :: keep_track
   real :: tau, tau1, tf
   real :: max_diff, max_diffk 
   integer :: do_average
   integer :: istep_initial
   logical :: check_tau
   integer :: check_convergence, check_convergencek
   real, dimension(:), allocatable :: keep_track_g
   logical :: adjoint
   logical :: exclude_source
   logical :: stop_convergence
   integer :: change_p
   logical :: last_call = .True.
   
   stop_convergence = .False.
   
   call init_stella 
   call init_convergence
   allocate(keep_track(0:nstep,2))
   allocate(keep_track_g(0:nstep))

   allocate(lag(14))
  !     call init_convergence (sum_omega, sum_phik)
!  sum_omega = 0
!  sum_gk = 0
  do_average = 0
  check_convergence = 0
  check_convergencek = 0
  istep_initial = 0
  
  if (debug) write (*,*) 'stella::diagnose_stella'
  call diagnose_stella (0)
  
  if (debug) write (*,*) 'stella::advance_stella'
  do istep = 1, nstep
     if (debug) write (*,*) 'istep = ', istep
     call advance_stella (istep)
     if (debug) write (*,*) '1'
     if (nsave > 0 .and. mod(istep,nsave)==0) then
        call scatter (kxkyz2vmu, gnew, gvmu)
        call stella_save_for_restart (gvmu, code_time, code_dt, istatus, fphi, fapar)
     end if
     call update_time
     call time_message(.false.,time_diagnostics,' diagnostics')
     call diagnose_stella (istep)
     call time_message(.false.,time_diagnostics,' diagnostics')
     if (debug) write (*,*) '2'
     !!Convergence test to get analytic form of g
     if (do_average == 0) then
        if (debug) write (*,*) 'get_omega'
        call get_omega (istep, max_diff)
        if (debug) write(*,*) 'get_gk'
        call get_gk (istep, max_diffk)
        if (abs(max_diff) < 1E-002 .and. abs(max_diffk) < 1E-002) then
           do_average = 1
           istep_initial = istep
        end if
     else
        write(*,*) 'first stage of convergence complete'
        write (*,*)
        call check_omega_convergence (istep, istep_initial,  check_convergence)
        call check_gk_convergence (istep, istep_initial, check_convergencek)
        if (check_convergence == 1 .and. check_convergencek == 1) then
           stop_convergence = .True.
           write(*,*) 'g_k =', g_k(1,1,1,1,1)
           write(*,*)
           if (stop_convergence) exit
        end if
     end if
     keep_track_g(nstep-istep) = abs(gnew(1,1,1,1,100))
     if (mod(istep,10)==0) call checkstop (stop_stella)
     if (stop_stella) exit
  end do

  deallocate(sum_omega)
!  deallocate(sum_gk)
  
  if (debug) write (*,*) 'stella::finish_stella'
  call finish_stella (adjoint = .False.)

  tau = code_time
  code_time = 0.
  adjoint = .True.

  change_p = 0
  
  call init_stella (change_p, adjoint)

  exclude_source = .False.

  lambda_old = 0.
  
  call get_tf (tf)
  if (debug) write(*,*) 'tf =', tf, 'tau =', tau
  call get_tau1 (tau, tau1, tf, check_tau)

  if (check_tau) then
     call get_g_tau (tau)
     if (debug) write(*,*) 'g_tau =', g_tau(1,1,1,1,100)
!     call allocate_lam_beta
!     if (debug) write(*,*) 'allocate_lam_beta'
     ! call adjoint_initial_condition
     ! if (debug) write(*,*) 'adjoint_initial_condition'
     call diagnose_stella (0, adjoint)
     lambda_old = 1E-004
     lambda_new = 1E-004 
     write(*,*) 'mark2, lambda_new =', lambda_new(1,1,1,1,100)
     ! call perturb_p(dG_dp, dQ_dp)
     ! if (debug) write(*,*) 'peturb p'
     ! call init_stella
     ! call lagrangian_term1 (dG_dp, sum_lag)
     ! if (debug) write(*,*) 'lagrangian term 1- stage 1'
     ! call lagrangian_term2 (dQ_dp, sum_beta_q)
     ! if (debug) write(*,*) 'lagrangian term 2- stage 1'
     keep_track(0,1) = abs(lambda_old(1,1,1,1,100))
     keep_track(0,2) = 0
     call g_analytic (tf)
     keep_track_g(0) = abs(gnew (1,1,1,1,100))
     write(*,*) 'lambda_old =', lambda_old(1,1,1,1,100), 'lambda_new =', lambda_new(1,1,1,1,100), 'gnew =', gnew(1,1,1,1,100)
     write(*,*)
!     keep_track(1,2) = abs(gnew(1,1,1,1,100))
!     call scatter (kxkyz2vmu, lambda_old, gvmu)
!     write(*,*) 'step1 scatter'
     do istep = 1, nstep
        write(*,*) 'mark, lambda_new =', lambda_new(1,1,1,1,100)
        if (code_time .LE. tf - tau1) then
!           phi = 0
!           phi_old = 0
!           gvmu = lambda_vmu
           if (debug) write (*,*) 'istep = ', istep
!           call advance_stella (istep)
!           if (debug) write (*,*) 'stella::advance_stella'
           call g_analytic (tf)
           if (debug) write(*,*) 'g analytic'
           write(*,*) 'starting adjoint'
           !           call advance_adjoint (tf, tau, istep, exclude_source)
           call advance_adjoint (tf, tau, istep)
           write(*,*) 'lambda_new =', lambda_new(1,1,1,1,100), 'gnew=',gnew(1,1,1,1,100)
           write(*,*)
!           call adjoint (tf, tau1, istep)
           if (debug) write(*,*) 'adjoint advance'
           if (nsave > 0 .and. mod(istep,nsave)==0) then
              write(*,*) 'scatter adjoint'
              call scatter (kxkyz2vmu, lambda_new, gvmu)
              write(*,*) 'save_for restart adjoint'
              call stella_save_for_restart (gvmu, code_time, code_dt, istatus, fphi, fapar)
           end if
           write(*,*) 'assign'
           keep_track(istep,1) = abs(lambda_new(1,1,1,1,100))
           write(*,*) 'lambda_new written'
!           keep_track(istep,2) = abs(gnew(1,1,1,1,100))
           keep_track(istep,2) = istep
           write(*,*) 'gnew written', gnew(1,1,1,1,100), 'abs=', abs(gnew(1,1,1,1,100))
           keep_track_g(istep) = abs(gnew(1,1,1,1,100))
           write(*,*) 'istep written'
           ! call perturb_p(dG_dp, dQ_dp)
           ! call lagrangian_term1 (dG_dp, sum_lag)
           ! call lagrangian_term2 (dQ_dp, sum_beta_q)
           write(*,*) 'update time'
           call update_time
           write(*,*) 'code time =', code_time
           write(*,*)
           call time_message(.false.,time_diagnostics,' diagnostics')
!           gnew = lambda_new
           call diagnose_stella (istep, adjoint)
!           lambda_new = gnew
           !           lambda_vmu = gvmu
        else 
           exclude_source = .True.
           phi = 0
           phi_old = 0
           if (debug) write (*,*) 'istep = ', istep
           call g_analytic (tf)
           if (debug) write(*,*) 'g analytic'
           write(*,*) 'starting adjoint'
           call advance_adjoint (tf, tau, istep, exclude_source)
           write(*,*) 'lambda_new =', lambda_new(1,1,1,1,100), 'gnew=',gnew(1,1,1,1,100)
           write(*,*)
           if (debug) write(*,*) 'adjoint advance'
           if (nsave > 0 .and. mod(istep,nsave)==0) then
              write(*,*) 'scatter adjoint'
              call scatter (kxkyz2vmu, lambda_new, gvmu)
              write(*,*) 'save_for restart adjoint'
              call stella_save_for_restart (gvmu, code_time, code_dt, istatus, fphi, fapar)
           end if
           write(*,*) 'assign'
           keep_track(istep,1) = abs(lambda_new(1,1,1,1,100))
           write(*,*) 'lambda_new written'
           keep_track(istep,2) = istep
           write(*,*) 'gnew written', gnew(1,1,1,1,100), 'abs=', abs(gnew(1,1,1,1,100))
!           keep_track_g(istep) = abs(gnew(1,1,1,1,100))
           write(*,*) 'istep written'
           write(*,*) 'update time'
           call update_time
           write(*,*) 'code time =', code_time
           write(*,*)
           call time_message(.false.,time_diagnostics,' diagnostics')
           call diagnose_stella (istep, adjoint = .True.)
        end if
     end do
      
      if (debug) write (*,*) 'stella::finish_stella'
      call finish_stella (last_call, adjoint)
      
   else
      write (*,*) 'final time not long enough'
      write (*,*)
   end if

!   lag = sum_lag + sum_beta_q

   open(1, file = 'data1.dat', status='new')  
   do istep = 0, nstep  
      write(1,*) keep_track(istep,1), keep_track(istep,2), keep_track_g(istep)
   end do  
   close(1) 

!   call deallocate_adjoint
   
contains
  
  subroutine deallocate_adjoint
    use dist_fn_arrays, only : lambda_old, beta, g_tau, g_k
    deallocate(lag)
    deallocate(sum_lag)
    deallocate(sum_beta_q)
    deallocate(dG_dp)
    deallocate(dQ_dp)
    deallocate(g_tau)
    deallocate(g_k)
    deallocate(beta)
    deallocate(lambda_old)
  end subroutine deallocate_adjoint

  
  subroutine init_convergence 

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_diagnostics, only: navg
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: g_tau, g_k, omega
    
    allocate(g_tau(naky,nakx,-nzgrid:nzgrid, ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(g_k(naky,nakx,-nzgrid:nzgrid, ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(omega(naky,nakx))
  end subroutine init_convergence


  ! subroutine allocate_lam_beta

  !   use zgrid, only: nzgrid, ntubes
  !   use stella_layouts, only: vmu_lo
  !   use kt_grids, only: naky, nakx
  !   use dist_fn_arrays, only: lambda_old, beta, lambda_new, lambda_vmu
  !   use stella_layouts, only: kxkyz_lo
  !   use vpamu_grids, only: nvpa, nmu

  !   allocate(lambda_vmu(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
  !   allocate(lambda_old(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !   allocate(lambda_new(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !   allocate(beta(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc))
  ! end subroutine allocate_lam_beta

  subroutine adjoint_initial_condition

    use dist_fn_arrays, only: lambda_old, lambda_vmu, lambda_new 
    use lambda_subroutines, only: get_beta
    
    lambda_vmu = 0
    lambda_old = 0
    call get_beta (lambda_old)
    
  end subroutine adjoint_initial_condition


  ! subroutine adjoint (tf, tau, istep)

  !   use lambda_subroutines, only: advance_adjoint
  !   use dist_fn_arrays, only: lambda_old, beta, gnew

  !   real, intent (in) :: tf, tau
  !   integer, intent (in) :: istep
    
  !   call advance_adjoint (tf, tau, istep)
  !   write (*,*) 'lambda_new =',  lambda_new (1,1,1,1,100), 'gnew =', gnew(1,1,1,1,100), 'beta =', beta(1)
  !   write(*,*)
  ! end subroutine adjoint

  subroutine lagrangian_term1 (dG_dp, sum_lag)

    use dist_fn_arrays, only: lambda_old
    use stella_layouts, only: kxkyz_lo
    use vpamu_grids, only: nvpa, nmu
    use dist_redistribute, only: kxkyz2vmu
    use redistribute, only: gather, scatter
    use kt_grids, only: naky, nakx
    use species, only: nspec
    use stella_time, only : code_time
    use redistribute, only: gather, scatter
    use advancing_g_analytic, only: get_chi
    
    complex, dimension(:,:,:), allocatable :: lambda_v
    complex, dimension(:,:,:,:), allocatable, intent(in) :: dG_dp
    complex, dimension(:), allocatable, intent (in out) :: sum_lag
    complex, dimension(:), allocatable :: sum_lag1
    complex, dimension(:,:,:,:,:), allocatable :: chi
    complex, dimension(:,:,:), allocatable :: chi_v
    complex, dimension(:), allocatable :: lag
    
    integer :: ikxkyz, q, imu, change_p

    allocate(lambda_v(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(lag(14))
    allocate(sum_lag1(14))
    if (.not. allocated(sum_lag)) allocate(sum_lag(14))
    allocate(chi_v(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))


    if (code_time .LE. tf -  tau) then
       call scatter (kxkyz2vmu, lambda_old, lambda_v)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          do imu = 1, nmu
             do q = 1, nvpa
                do change_p = 1, 14
                   lag(change_p)= lambda_v(q,imu,ikxkyz)*dG_dp(ikxkyz,q,imu,change_p)
                end do
             end do
          end do
       end do
         
    else if (code_time .GT. tf - tau) then
       call get_chi (chi, tf, tau)
       call scatter (kxkyz2vmu, chi, chi_v)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          do imu = 1, nmu
             do q = 1, nvpa
                do change_p = 1, 14
                   lag(change_p)= -chi_v(q,imu,ikxkyz)*dG_dp(ikxkyz,q,imu,change_p)
                end do
             end do
          end do
       end do

    end if

    sum_lag1 = sum_lag + lag
    sum_lag = sum_lag1
	
    deallocate(sum_lag1)
    deallocate(lag)

  end subroutine lagrangian_term1

  
  subroutine lagrangian_term2 (dQ_dp, sum_beta_q)
    
    use dist_fn_arrays,only: beta
    use stella_layouts, only: kxkyz_lo

    complex, dimension(:,:), allocatable, intent (in) :: dQ_dp
    complex, dimension(:), allocatable :: beta_dQ_dp
    complex, dimension(:), allocatable, intent(in out) :: sum_beta_q
    complex, dimension(:), allocatable :: sum_beta_q1
   
    integer :: ikxkyz, change_p

    if (.not. allocated(sum_beta_q)) allocate(sum_beta_q(14))
    allocate(sum_beta_q1(14))
    allocate(beta_dQ_dp(14))
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       do change_p = 1,14
          beta_dQ_dp(change_p) = beta(ikxkyz)*dQ_dp(ikxkyz,change_p)
       end do
    end do

    sum_beta_q1= sum_beta_q + beta_dQ_dp
    sum_beta_q = sum_beta_q1

    deallocate(sum_beta_q1)
    deallocate(beta_dQ_dp)
    
  end subroutine lagrangian_term2

  subroutine get_tf (tf)
    use init_g, only: tstart
    use run_parameters, only: nstep
    use stella_time, only: code_dt

    real, intent (out) :: tf
    
    tf = tstart + code_dt*nstep
    
  end subroutine get_tf

    
  subroutine init_stella (change_p, adjoint)

    use mp, only: init_mp, broadcast
    use mp, only: proc0
    use file_utils, only: init_file_utils
    use file_utils, only: run_name
    use job_manage, only: checktime, time_message
    use physics_parameters, only: init_physics_parameters
    use physics_flags, only: init_physics_flags
    use physics_flags, only: nonlinear, full_flux_surface, include_parallel_nonlinearity
    use run_parameters, only: init_run_parameters
    use run_parameters, only: avail_cpu_time, nstep
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use species, only: init_species, read_species_knobs
    use species, only: nspec
    use zgrid, only: init_zgrid
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: init_geometry
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac
    use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
    use response_matrix, only: init_response_matrix
    use init_g, only: ginit, init_init_g
    use fields, only: init_fields, advance_fields
    use stella_time, only: init_tstart
    use init_g, only: tstart
    use stella_diagnostics, only: init_stella_diagnostics
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gnew
    use dist_fn, only: init_gxyz, init_dist_fn
    use time_advance, only: init_time_advance
    use extended_zgrid, only: init_extended_zgrid
    use kt_grids, only: init_kt_grids, read_kt_grids_parameters
    use kt_grids, only: naky, nakx, ny, nx, nalpha
    use vpamu_grids, only: init_vpamu_grids, read_vpamu_grids_parameters
    use vpamu_grids, only: nvgrid, nmu
    use stella_transforms, only: init_transforms

    implicit none
    logical :: exit, list, restarted
    character (500), target :: cbuff
    integer, optional, intent (in) :: change_p
    logical, optional, intent (in) :: adjoint

    ! initialize mpi message passing                                                                                
    if (.not.mpi_initialized) call init_mp
    mpi_initialized = .true.

    debug = debug .and. proc0

    if (debug) write (*,*) 'stella::init_stella::check_time'
    ! initialize timer                                                                                              
    call checktime(avail_cpu_time,exit)

    if (proc0) then
       if (debug) write (*,*) 'stella::init_stella::write_start_message'
       ! write message to screen with useful info                                                                   
       ! regarding start of simulation                                                                              
       call write_start_message
       if (debug) write (*,*) 'stella::init_stella::init_file_utils'
       ! initialize file i/o                                                                                        
       call init_file_utils (list)
       call time_message(.false.,time_total,' Total')
       call time_message(.false.,time_init,' Initialization')
       cbuff = trim(run_name)
    end if

    call broadcast (cbuff)
    if (.not. proc0) run_name => cbuff

    if (debug) write(6,*) "stella::init_stella::init_physics_flags"
    call init_physics_flags
    if (debug) write(6,*) "stella::init_stella::init_physics_parameters"
    call init_physics_parameters
    if (debug) write(6,*) "stella::init_stella::init_zgrid"
    call init_zgrid
    if (debug) write (6,*) "stella::init_stella::read_species_knobs"
    call read_species_knobs
    if (debug) write (6,*) "stella::init_stella::read_kt_grids_parameters"
    call read_kt_grids_parameters
    if (debug) write (6,*) "stella::init_stella::read_vpamu_grids_parameters"
    call read_vpamu_grids_parameters
    if (debug) write (6,*) "stella::init_stella::init_dist_fn_layouts"
    call init_dist_fn_layouts (nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
    if (nonlinear .or. full_flux_surface .or. include_parallel_nonlinearity) then
       if (debug) write (*,*) "stella::init_stella::init_transforms"
       call init_transforms
    end if
    if (debug) write(6,*) "stella::init_stella::init_geometry"
    if (present (change_p)) then
       call init_geometry (change_p)
    else
       call init_geometry
    end if
!    call init_geometry
    if (debug) write (6,*) 'stella::init_stella::init_species'
    call init_species
    if (debug) write(6,*) "stella::init_stella::init_init_g"
    call init_init_g
    if (debug) write(6,*) "stella::init_stella::init_run_parameters"
    call init_run_parameters
    if (debug) write (6,*) 'stella::init_stella::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (6,*) 'stella::init_stella::init_kt_grids'
    call init_kt_grids (geo_surf, twist_and_shift_geo_fac)
    if (debug) write (6,*) 'stella::init_stella::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write(6,*) "stella::init_stella::init_dist_fn"
    if (present(adjoint)) then 
       call init_dist_fn (adjoint)
    else
       call init_dist_fn
    end if
    if (debug) write (6,*) 'stella::init_stella::init_extended_zgrid'
    call init_extended_zgrid
    if (debug) write (6,*) 'stella::init_stella::init_fields'
    call init_fields
    if (debug) write (6,*) 'stella::init_stella::init_time_advance'
    call init_time_advance
    if (debug) write(6,*) "stella::init_stella::ginit"
    call ginit (restarted)
    if (debug) write(6,*) "stella::init_stella::init_gxyz"
    if (present(adjoint)) then 
       call init_gxyz (adjoint)
    else
       call init_gxyz
    end if
    if (debug) write(6,*) "stella::init_stella::init_response_matrix"
    if (stream_implicit .or. driftkinetic_implicit) call init_response_matrix

    if (.not.restarted) then
       if (debug) write (6,*) 'stella::init_stella::get_fields'
       ! get initial field from initial distribution function
       if (present(adjoint)) then
          call advance_fields (lambda_new, phi, apar, dist='gbar')
       else
          call advance_fields (gnew, phi, apar, dist='gbar')
       end if
    end if

    if (debug) write (6,*) 'stella::init_stella::init_stella_diagnostics'
    call init_stella_diagnostics (nstep)
    if (debug) write (6,*) 'stella::init_stella::init_tstart'
    call init_tstart (tstart)

    if (proc0) call time_message(.false.,time_init,' Initialization')

  end subroutine init_stella

  subroutine write_start_message

    use mp, only: proc0, nproc

    implicit none

    if (proc0) then
       if (nproc==1) then
          write (*,*) "Running on ", nproc, " processor"
       else
          write (*,*) "Running on ", nproc, " processors"
       end if
       write (*,*)
    end if

  end subroutine write_start_message

  subroutine finish_stella (last_call, adjoint)

    use mp, only: finish_mp
    use mp, only: proc0
    use file_utils, only: finish_file_utils
    use job_manage, only: time_message
    use physics_parameters, only: finish_physics_parameters
    use physics_flags, only: finish_physics_flags
    use run_parameters, only: finish_run_parameters
    use zgrid, only: finish_zgrid
    use species, only: finish_species
    use time_advance, only: time_gke, time_parallel_nl
    use time_advance, only: finish_time_advance
    use parallel_streaming, only: time_parallel_streaming
    use mirror_terms, only: time_mirror
    use dissipation, only: time_collisions
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use fields, only: time_field_solve
    use stella_diagnostics, only: finish_stella_diagnostics
    use response_matrix, only: finish_response_matrix
    use stella_geometry, only: finish_geometry
    use extended_zgrid, only: finish_extended_zgrid
    use vpamu_grids, only: finish_vpamu_grids
    use kt_grids, only: finish_kt_grids

    implicit none

    logical, intent (in), optional :: last_call
    logical, optional, intent (in) :: adjoint

    if (debug) write (*,*) 'stella::finish_stella::finish_stella_diagnostics'
    call finish_stella_diagnostics
    if (debug) write (*,*) 'stella::finish_stella::finish_response_matrix'
    call finish_response_matrix
    if (debug) write (*,*) 'stella::finish_stella::finish_fields'
    call finish_fields
    if (debug) write (*,*) 'stella::finish_stella::finish_time_advance'
    call finish_time_advance
    if (debug) write (*,*) 'stella::finish_stella::finish_extended_zgrid'
    call finish_extended_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_dist_fn'
    if (present(adjoint)) then 
       call finish_dist_fn (adjoint)
    else
       call finish_dist_fn
    end if
    if (debug) write (*,*) 'stella::finish_stella::finish_vpamu_grids'
    call finish_vpamu_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_kt_grids'
    call finish_kt_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_run_parameters'
    call finish_run_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_species'
    call finish_species
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_flags'
    call finish_physics_flags
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_parameters'
    call finish_physics_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_geometry'
    call finish_geometry
    if (debug) write (*,*) 'stella::finish_stella::finish_zgrid'
    call finish_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_file_utils'
    if (proc0) then
       call finish_file_utils
       call time_message(.false.,time_total,' Total')
       write (*,*)
       write (*,fmt=101) 'initialization:', time_init(1)/60., 'min'
       write (*,fmt=101) 'diagnostics:', time_diagnostics(1)/60., 'min'
       write (*,fmt=101) 'fields:', time_field_solve(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_field_solve(1,2)/60., 'min'
       write (*,fmt=101) 'mirror:', time_mirror(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_mirror(1,2)/60., 'min'
       write (*,fmt=101) 'stream:', time_parallel_streaming(1)/60., 'min'
       write (*,fmt=101) 'dgdx:', time_gke(1,5)/60., 'min'
       write (*,fmt=101) 'dgdy:', time_gke(1,4)/60., 'min'
       write (*,fmt=101) 'wstar:', time_gke(1,6)/60., 'min'
       write (*,fmt=101) 'collisions:', time_collisions(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_collisions(1,2)/60., 'min'
       write (*,fmt=101) 'ExB nonlin:', time_gke(1,7)/60., 'min'
       write (*,fmt=101) 'parallel nonlin:', time_parallel_nl(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_parallel_nl(1,2)/60., 'min'
       write (*,fmt=101) 'total implicit: ', time_gke(1,9)/60., 'min'
       write (*,fmt=101) 'total explicit: ', time_gke(1,8)/60., 'min'
       write (*,fmt=101) 'total:', time_total(1)/60., 'min'
       write (*,*)
    end if
101 format (a17,0pf8.2,a4)

    if (debug) write (*,*) 'stella::finish_stella::finish_mp'
    ! finish (clean up) mpi message passing                                                                         
    if (present(last_call)) then
       call finish_mp
       mpi_initialized = .false.
    end if

  end subroutine finish_stella

end program stella_adjoint

