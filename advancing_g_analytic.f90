module advancing_g_analytic

  public :: g_analytic
  public :: s_term
  public :: get_chi
  public :: get_tau1
  public :: get_g_tau
  
  private
  
contains

  !!get tau
  subroutine get_tau1 (tau, tau1,  tf, check_tau)
    real, intent (in) :: tau
    real, intent (out) :: tau1
    logical, intent (out) :: check_tau
    real, intent (in) :: tf
    
    check_tau = .FALSE.
    tau1 = 0.4*tf
    if (tau<tau1) check_tau = .TRUE.
    
  end subroutine get_tau1


  !! Get analytic form of g
  subroutine g_analytic  (tf)

    use dist_fn_arrays, only : gnew, g_k, omega
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid,ntubes
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_time
    
    real, intent (in) :: tf
    real :: code_time1
    integer :: ivmu, ikx, iky, p, ia
    complex :: i
    
    ia = 1
    i = (0,1)
    code_time1 = tf - code_time
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ikx = 1, nakx
          do iky = 1, naky
             do p = -nzgrid, nzgrid
                gnew(iky,ikx,p,ia,ivmu) = g_k(iky,ikx,p,ia,ivmu) * exp(i*omega(iky,ikx)*code_time1)
             end do
          end do
       end do
    end do
    write(*,*) 'gnew original=',gnew(1,1,1,1,100), 'tf-code_time', tf - code_time
  end subroutine g_analytic

  subroutine get_g_tau (tau)
    
    use dist_fn_arrays, only: g_tau, g_k, omega
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, naky
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    integer :: ivmu, ikx, iky, p, ia
    complex :: i
    real, intent (in) :: tau
    i = (0,1)
    ia = 1    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ikx = 1, nakx
          do iky = 1, naky
             do p = -nzgrid, nzgrid
                g_tau(iky,ikx,p,ia,ivmu) = g_k(iky,ikx,p,ia,ivmu) * exp(i*omega(iky,ikx)*tau)
             end do
          end do
       end do
    end do

  end subroutine get_g_tau
  
  subroutine get_weights (xi_vpa,xi_mu,xi_z)
    use vpamu_grids, only: nvpa, vpa, vpa_max
    use vpamu_grids, only: nmu, mu, mu_max
    use zgrid, only: nzgrid, zed

    real, dimension (:), allocatable, intent (out) :: xi_vpa
    real, dimension (:), allocatable, intent (out) :: xi_mu
    real, dimension (:), allocatable, intent (out) :: xi_z
    integer :: p, q, imu
    
    allocate (xi_vpa(nvpa))
    allocate(xi_mu(nmu))
    allocate(xi_z(-nzgrid:nzgrid))

    do q = 1, nvpa
       xi_vpa(q) = exp(-1/abs(vpa(q)*vpa(q)-vpa_max*vpa_max))
    end do

    do imu = 1, nmu
       xi_mu(imu) = exp(-1/abs(mu(imu)*mu(imu) - mu_max*mu_max))
    end do
    
    do p = -nzgrid, nzgrid
       xi_z(p) = exp(-1/abs(zed(p)*zed(p) - zed(nzgrid)*zed(nzgrid)))
    end do

    write(*,*) 'xi_z =', xi_z(nzgrid/2), 'xi_mu =', xi_mu(nmu/2), 'xi_vpa', xi_vpa(nvpa/2)

  end subroutine get_weights

  
  subroutine get_chi (chi, tf, tau)

    use dist_fn_arrays, only : omega, gnew
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid,ntubes
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
    use dist_fn_arrays, only: g_tau

    real, dimension (:), allocatable :: xi_vpa
    real, dimension (:), allocatable :: xi_mu
    real, dimension (:), allocatable :: xi_z
    complex, dimension (:,:,:,:,:), allocatable, intent (out) :: chi
    complex, dimension (:,:,:,:,:), allocatable :: g_tau2, g2
    real, intent (in) :: tf, tau
    complex :: mod_g2, mod_g_tau2
    integer :: ivmu, imu, q, is, ikx, iky, p
    complex :: i
    
    i = (0,1)
    
    allocate (chi(naky,nakx,-nzgrid:nzgrid, ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate ( g2(naky,nakx,-nzgrid:nzgrid, ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate ( g_tau2(naky,nakx,-nzgrid:nzgrid, ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    call get_weights (xi_vpa,xi_mu,xi_z)
    call get_g_tau (tau)
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do ikx = 1, nakx
          do iky = 1, naky
             do p = -nzgrid, nzgrid
                g2(iky, ikx, p, 1, ivmu) = gnew(iky,ikx,p,1,ivmu)*gnew(iky,ikx,p,1,ivmu)
                g_tau2(iky,ikx,p,1,ivmu) = g_tau(iky,ikx,p,1,ivmu)*g_tau(iky,ikx,p,1,ivmu)
             end do
          end do
       end do
    end do

    mod_g2 = sum(g2)
    mod_g_tau2 = sum(g_tau2)
                

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do ikx = 1, nakx
          do iky = 1, naky
             do p = -nzgrid, nzgrid
                chi(iky,ikx,p,1,ivmu) = - 2*i /(tf-tau) *mod_g2&
                     * g_tau(iky,ikx,p,1,ivmu)/(mod_g_tau2**2) &
                     *xi_vpa(q) * xi_z(p) * xi_mu(imu)/(omega(iky,ikx))&
                     *(exp(i*omega(iky,ikx)*tf) - exp(i*omega(iky,ikx)*tau))
             end do
          end do
       end do
    end do
    
  end subroutine get_chi
  
  
  !! Get source term for lambda equation
  subroutine s_term (sterm, tf, tau)

    use dist_fn_arrays, only : gnew, g_tau
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid,ntubes
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    real, dimension (:), allocatable :: xi_vpa
    real, dimension (:), allocatable :: xi_mu
    real, dimension (:), allocatable :: xi_z
    complex, dimension (:,:,:,:,:), allocatable :: conj
    integer :: ivmu, ia, imu, q,is, ikx, iky, p
    real, intent (in) :: tf
    real, intent (in) :: tau
    complex, dimension(:,:,:,:,:), allocatable, intent (out) :: sterm
    
    allocate(sterm(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(conj(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    
    call get_weights (xi_vpa,xi_mu,xi_z)
    call get_g_tau (tau)
    
    ia = 1

    conj = conjg(g_tau)
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do ikx = 1, nakx
          do iky = 1, naky
             do p = -nzgrid, nzgrid
                sterm(iky,ikx,p,ia,ivmu) = (2/(tf-tau))*xi_vpa(q)*xi_z(p)*xi_mu(imu)*gnew(iky,ikx,p,ia,ivmu)/(g_tau(iky,ikx,p,ia,ivmu)*conj(iky,ikx,p,ia,ivmu))
             end do
          end do
       end do
    end do
    write(*,*) 'sterm original =', sterm(1,1,1,1,100), 'gnew=', gnew(1,1,1,1,100)
    ! do p = -nzgrid, nzgrid
    !    do q = 1, nvpa
    !       do imu = 1, nmu
    !          write(*,*) 'sterm original', sterm(p,q,imu)
    !       end do
    !    end do
    ! end do
  end subroutine s_term

  ! !! Get gamma lagrange multiplier
  ! subroutine get_gamma (avg_omega2,avg_gk2,integral)

  !   use dist_fn_arrays, only : gold, g_tau
  !   use kt_grids, only: nakx, naky
  !   use zgrid, only: nzgrid,ntubes
  !   use stella_layouts, only: vmu_lo
  !   use stella_time, only: code_dt, code_time
    
  !   complex, dimension (:,:,:,:,:), allocatable :: to_integrate
  !   complex, dimension (:,:,:,:,:), allocatable, intent (in out) :: integral
  !   complex, dimension (:,:), allocatable, intent (in) :: avg_omega2
  !   complex, dimension (:,:,:), allocatable, intent (in) :: avg_gk2
  !   complex, dimension (:,:,:,:,:), allocatable :: g_tf
  !   complex :: i
  !   integer :: p,q,is,imu,ia,iy,ikx,iky,ikxkyz,ivmu
  !   i=(0,1)
  !   ia=1
    
  !   allocate(to_integrate(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !   if(.not. allocated(integral)) allocate(integral(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !   allocate(gam(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !   allocate(g_tf(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    
  !   if (code_time < tau) then
  !      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
  !         do ikx = 1,nakx
  !            do iky =1,naky
  !               do p = -nzgrid,nzgrid               
  !                  to_integrate(iky,ikx,p,ia, ivmu) =&
  !                       2/(tf-tau)*gold(iky,ikx,p,ia,ivmu)**2/(g_tau(iky,ikx,p,ia,ivmu,1)**3)
  !                  integral(iky,ikx,p,ia,ivmu) = integral(iky,ikx,p,ia,ivmu) + &
  !                       to_integrate(iky,ikx,p,ia,ivmu)*code_dt
  !               end do
  !            end do
  !         end do
  !      end do
  !   else
  !      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
  !         do ikx = 1,nakx
  !            do iky =1,naky
  !               do p = -nzgrid,nzgrid
  !                  to_add(iky,ikx,p,ia,ivmu) = 2/(tf-tau)*1/(2*i*avg_omega2(iky,ikx))*(g_tf(iky,ikx,p,ia,ivmu)**2-g_tau(iky,ikx,p,ia,ivmu,1)**2)/(g_tau(iky,ikx,p,ia,ivmu,1)**3)
  !                  integral(iky,ikx,p,ia,ivmu) = integral(iky,ikx,p,ia,ivmu) + to_add(iky,ikx,p,ia,ivmu)
  !               end do
  !            end do
  !         end do
  !      end do
  !   end if

    
  ! end subroutine get_gamma


end module advancing_g_analytic
