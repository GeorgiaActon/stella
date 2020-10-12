module adjoint_convergence

  public :: get_omega
  public :: check_omega_convergence
  public :: check_phik_convergence
  public :: get_phik
  public :: get_avg_phi

  private

  logical :: write_omega

contains

  !!First convergence test
  subroutine get_omega (istep, max_diff)

    use kt_grids, only: nakx, naky
    use stella_time, only: code_time
    use stella_diagnostics, only: omega_vs_time, navg

    complex, dimension(:,:), allocatable :: diff
    integer, intent (in) :: istep
    integer :: ikx, iky
    real, intent (out) :: max_diff

    allocate(diff(naky,nakx))

     do ikx = 1, nakx
       do iky = 1, naky
          if (mod(istep,navg)+1 > 1) then 
             diff(iky, ikx) = omega_vs_time(mod(istep,navg)+1,iky,ikx)- omega_vs_time(mod(istep,navg),iky,ikx)
          else 
             diff(iky,ikx) = omega_vs_time(1,iky,ikx)- omega_vs_time(50,iky,ikx)
          end if
       end do
    end do

    max_diff = maxval(abs(diff))

    deallocate(diff)

  end subroutine get_omega


  subroutine get_phik (istep, max_diffk)

    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid
    use stella_diagnostics, only: omega_vs_time, navg
    use fields_arrays, only: phi, phik
    use stella_time, only: code_time

    implicit none 

    complex, dimension(:,:,:), allocatable :: diff
    integer, intent (in) :: istep
    integer :: ikx, iky, p 
    real, intent (out) :: max_diffk
    complex :: i

    i = (0,1)

    allocate(diff(naky,nakx,-nzgrid:nzgrid))

    do iky = 1, naky
       do ikx = 1, nakx
          do p = -nzgrid, nzgrid
             phik(iky,ikx,p,mod(istep,navg)+1) = phi(iky,ikx,p,1)/exp(-i*omega_vs_time(mod(istep,navg)+1,iky,ikx)*code_time)
             if (mod(istep,navg)+1 > 1) then
                diff(iky, ikx, p) = phik(iky,ikx,p,mod(istep,navg)+1) - phik(iky,ikx,p,mod(istep,navg))
             else
                diff(iky, ikx, p) = phik(iky,ikx,p,1) -phik(iky,ikx,p,50)
             end if
          end do
       end do
    end do

    max_diffk = maxval(abs(diff))

    deallocate(diff)

  end subroutine get_phik


  !!Second convergence test
  
  subroutine check_omega_convergence (istep, istep_initial, sum_omega, check_convergence, avg_omega2)

    use stella_diagnostics, only: omega_vs_time, navg
    use kt_grids, only: nakx, naky

    integer, intent (in) :: istep, istep_initial
    integer :: ikx, iky, p
    complex, dimension (:,:), allocatable :: avg_omega
    complex, dimension (:,:), allocatable, intent (out) :: avg_omega2
    complex, dimension (:,:), intent (in out) :: sum_omega
    complex, dimension (:,:), allocatable :: sum_omega2
    complex, dimension(:,:), allocatable :: diff_omega
    real :: max_diff_omega 
    integer, intent (out) :: check_convergence

    allocate(avg_omega(naky,nakx))
    allocate(diff_omega(naky,nakx))
    allocate(sum_omega2(naky,nakx))

    avg_omega = 0 

    do ikx = 1, nakx
       do iky = 1, naky
          sum_omega(iky,ikx) = sum_omega(iky,ikx) + omega_vs_time(mod(istep,navg)+1,iky,ikx)
       end do
    end do


    if (istep-istep_initial > navg + 20) then 
    
       avg_omega = sum_omega/(istep- istep_initial)

       do iky = 1, naky
          do ikx = 1, nakx
             sum_omega2(iky, ikx) = sum(omega_vs_time(:,iky,ikx))
          end do
       end do
       avg_omega2 = sum_omega2/navg

       diff_omega = avg_omega - avg_omega2
       
       max_diff_omega = maxval(abs(diff_omega))

       if (abs(max_diff_omega) < 1E-004) then 
          check_convergence = 1 
       else 
          check_convergence = 0
       end if

     end if

    deallocate(diff_omega, avg_omega, sum_omega2)

  end subroutine check_omega_convergence

  
  subroutine check_phik_convergence (istep, istep_initial, sum_phik, check_convergencek, avg_phik2)
    use fields_arrays, only: phi, phik
    use stella_diagnostics, only: omega_vs_time, navg
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid
    use stella_time, only: code_time

    complex, dimension(:,:,-nzgrid:), intent (in out) :: sum_phik
    complex, dimension(:,:,:), allocatable :: sum_phik2
    complex, dimension (:,:,:), allocatable :: avg_phik
    complex, dimension (:,:,:), allocatable, intent (out) :: avg_phik2
    complex, dimension (:,:,:), allocatable :: diff_phik
    integer, intent (in) :: istep, istep_initial
    integer :: ikx, iky, p
    real :: max_diff_phik
    integer, intent (out) :: check_convergencek 
    complex :: i

    i = (0,1)

    allocate(avg_phik(naky,nakx,-nzgrid:nzgrid))
    allocate(diff_phik(naky, nakx, -nzgrid:nzgrid))
    allocate(avg_phik2(naky,nakx,-nzgrid:nzgrid))
    allocate(sum_phik2(naky,nakx,-nzgrid:nzgrid))

    avg_phik = 0

    do ikx = 1, nakx
       do iky = 1, naky
          do p = -nzgrid, nzgrid
             phik(iky,ikx,p,1) = phi(iky,ikx,p,1)/exp(-i*omega_vs_time(mod(istep,navg)+1,iky,ikx)*code_time)
             sum_phik(iky,ikx,p) = sum_phik(iky,ikx,p) +  phik(iky,ikx,p,1)
          end do
       end do
    end do

    if (istep-istep_initial > navg+ 20) then

       avg_phik = sum_phik/(istep - istep_initial)  

       do iky = 1, naky
          do ikx = 1, nakx
             do p = -nzgrid, nzgrid
                sum_phik2 (iky,ikx,p) = sum(phik(iky,ikx,p,:))
             end do
          end do
       end do

       avg_phik2 = sum_phik2/ navg
       diff_phik = avg_phik - avg_phik2

       max_diff_phik = maxval(abs(diff_phik))

       write(*,*) avg_phik(1,1,1), avg_phik2(1,1,1), max_diff_phik
       
       if(max_diff_phik < 1E-004) then 
          check_convergencek = 1
       else
          check_convergencek = 0
       end if

    end if

    deallocate(diff_phik, avg_phik, sum_phik2)

  end subroutine check_phik_convergence


  subroutine get_avg_phi (avg_phik, avg_omega)
    use fields_arrays, only: phi_steady
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid

    complex, dimension (:,:,-nzgrid:), intent (in) :: avg_phik
    complex, dimension (:,:), intent (in) :: avg_omega

    integer :: ikx, iky, p

    do iky = 1, naky
       do ikx =  1, nakx
          do p = -nzgrid, nzgrid
             phi_steady(iky,ikx,p,1) = avg_phik(iky,ikx,p)
          end do
       end do
    end do

  end subroutine get_avg_phi




end module adjoint_convergence
