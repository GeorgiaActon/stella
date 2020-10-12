module adjoint_p_derviatives


implicit none

public :: perturb_p

private


real, dimension (:,:), allocatable :: j0_initial
real, dimension (:,:), allocatable :: j1_initial
complex, dimension (:,:,:), allocatable :: dgdvpa0
real, dimension(:,:,:,:,:), allocatable :: darg0dz
complex, dimension(:,:,:,:,:), allocatable :: wdriftx_g0, wdrifty_g0
complex, dimension(:,:,:,:,:), allocatable :: wdriftx_phi0, wdrifty_phi0
complex, dimension(:,:,:,:,:), allocatable :: wstarv0
real, dimension (:), allocatable :: gradpar0
real, dimension (:,:), allocatable :: bmag0
complex, dimension(:,:,:), allocatable :: gnew0
real, dimension(:,:,:,:,:), allocatable :: arg0
real, dimension (:,:), allocatable :: dbdzed0
complex, dimension (:,:,:), allocatable :: dgdz0_v
complex, dimension (:), allocatable :: mu0
complex, dimension(:,:,:), allocatable :: phi0
complex, dimension (:,:,:,:), allocatable :: dphidz0
real, dimension (:), allocatable :: vpa0
real, dimension (:,:,:), allocatable :: dJ0_dz_exp_unpert
real, dimension (:,:,:), allocatable :: J0_exp_unpert
real, dimension (:,:,:), allocatable :: gamtot0
complex, dimension (:,:,:,:,:), allocatable :: dgdz0

contains 

  subroutine unperturbed_p

    use stella_geometry, only: dbdzed
    use vpamu_grids, only: vpa, vperp2
    use dist_fn_arrays, only: gnew
    use zgrid, only:nzgrid, delzed, ntubes
    use fields_arrays, only: phi
    use redistribute, only: gather, scatter
    use stella_layouts, only: vmu_lo, imu_idx, iv_idx
    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use gyro_averages, only: aj0v, aj1v
    use dist_fn_arrays, only: wstar
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use mirror_terms, only: get_dgdvpa_explicit
    use species, only: spec, nspec
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, naky, nalpha
    use dist_redistribute, only: kxkyz2vmu
    use stella_geometry, only: bmag, gradpar
    use dist_fn_arrays, only: kperp2
    use init_and_finish, only: allocate_stella, deallocate_stella
    use stella_geometry, only: get_dzed, get_dzed_complex
    use fields, only: gamtot
    use vpamu_grids, only:  maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nvpa, nmu
 
    complex, dimension(:,:,:), allocatable :: gnew_v
    
    integer :: p,q,is,imu,ia,ikx,iky,ikxkyz,ivmu

    allocate(gnew_v(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(j0_initial(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(j1_initial(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(dgdz0_v(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(dgdz0(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(dgdvpa0(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(dphidz0(naky,nakx,-nzgrid:nzgrid, ntubes))
    allocate(dbdzed0(nalpha,-nzgrid:nzgrid))
    allocate(darg0dz(nspec,-nzgrid:nzgrid,nmu,nakx,naky))
    allocate(wdriftx_g0(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdrifty_g0(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdriftx_phi0(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdrifty_phi0(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wstarv0(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(arg0(nspec,-nzgrid:nzgrid,nmu,nakx,naky))
    allocate(gnew0(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(gradpar0(-nzgrid:nzgrid))
    allocate(bmag0(nalpha,-nzgrid:nzgrid))
    allocate(phi0(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(mu0(nmu))
    allocate(vpa0(nvpa))
    allocate(dJ0_dz_exp_unpert(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(J0_exp_unpert(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate(gamtot0(naky,nakx,-nzgrid:nzgrid))
    

    call allocate_stella
    
    ia = 1

    gradpar0 = gradpar
    bmag0 = bmag

    gamtot0 = gamtot

    !! reallocate drift velocities
    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do p = -nzgrid,nzgrid
          wdriftx_g0(ia,p,q,imu,is)  = wdriftx_g(ia,p,ivmu)
          wdrifty_g0(ia,p,q,imu,is) = wdrifty_g(ia,p,ivmu)
          wdriftx_phi0(ia,p,q,imu,is) = wdriftx_phi(ia,p,ivmu)
          wdrifty_phi0(ia,p,q,imu,is) = wdrifty_phi(ia,p,ivmu)
          wstarv0(ia,p,q,imu,is) = wstar(ia,p,ivmu)
       end do
    end do

    
    !! reallocate bessel functions
    j0_initial = aj0v
    j1_initial = aj1v
    
    !! get dg/dz
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
       do ikx = 1, nakx
          do iky = 1, naky
             call get_dzed_complex (nzgrid, delzed, gnew(iky,ikx,:,ia,ivmu), dgdz0(iky,ikx,:,ia,ivmu))
          end do
       end do
    end do

    call scatter (kxkyz2vmu, dgdz0, dgdz0_v)

    !! get db/dz
    dbdzed0 = dbdzed

    !! get dphi/dz
    do ikx = 1, nakx
       do iky = 1, naky
          call get_dzed_complex (nzgrid, delzed, phi(iky,ikx,:,ia), dphidz0(iky,ikx,:,ia))
       end do
    end do    

    !! get dg/dvpa
    call scatter (kxkyz2vmu, gnew, gnew_v)
    gnew0 = gnew_v
    dgdvpa0 = gnew_v
    call get_dgdvpa_explicit (dgdvpa0)

    do iky = 1, naky
       do ikx = 1, nakx
          do is = 1, nspec
             do imu = 1, nmu
                do p = -nzgrid, nzgrid
                   arg0(is,p,imu,ikx,iky) = spec(is)%smz*sqrt(vperp2(ia,p,imu)*kperp2(iky,ikx,ia,p))/bmag(ia,p)
                   call get_dzed (nzgrid, delzed, arg0(is,:,imu,ikx,iky), darg0dz(is,:,imu,ikx,iky))
                end do
             end do
          end do
       end do
    end do

    !! parallel streaming terms
    vpa0=vpa 

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       p = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do q = 1, nvpa
             dJ0_dz_exp_unpert(q,imu,ikxkyz) = arg0(is,p,imu,ikx,iky)*aj1v(imu,ikxkyz)*maxwell_vpa(q)*maxwell_mu(1,p,imu)
             J0_exp_unpert(q,imu,ikxkyz) = aj0v(imu,ikxkyz)*maxwell_vpa(q)*maxwell_mu(1,p,imu)
          end do
       end do
    end do
    
    call deallocate_stella

    deallocate(gnew_v)

  end subroutine unperturbed_p

  subroutine perturb_p (dG_dp, dQ_dp)

    use zgrid, only: ntubes
    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use gyro_averages, only: aj0v
    use init_and_finish, only: allocate_stella, deallocate_stella 
    use species, only: spec, nspec
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, naky
    use fields, only: gamtot
    use vpamu_grids, only: integrate_vmu
    use millerlocal, only: del

    complex, dimension(:,:,:), allocatable :: dmirror_dp
    complex, dimension (:,:,:), allocatable :: prll_strm_term
    complex, dimension(:,:,:), allocatable :: star_term
    complex, dimension(:,:,:), allocatable :: vm_g_term, vm_phi_term
    
    complex, dimension(:,:,:,:), allocatable, intent (out) :: dG_dp
    !preal, dimension (:,:,:), allocatable :: gamtot0
    complex, dimension (:,:,:), allocatable :: to_avg1, to_avg2
    complex, dimension(:), allocatable :: to_sum1, to_sum2
    complex, dimension(:,:), allocatable, intent (out) :: dQ_dp
    complex, dimension(:,:), allocatable :: dJ0_dp
    complex, dimension(:,:,:), allocatable :: dgam_dp
    
    integer :: change_p
    integer :: ikx, iky, ikxkyz, p ,is, q , imu

    call deallocate_stella
    call unperturbed_p 
    
    do change_p = 1, 14
       
       write(*,*) del(change_p)
       write(*,*)
       
       call allocate_stella (change_p)

       if(.not. allocated(dG_dp)) allocate(dG_dp(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nvpa,nmu,14))
       if(.not. allocated(to_avg1)) allocate(to_avg1(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       if(.not. allocated(to_avg2)) allocate(to_avg2(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       if(.not. allocated(to_sum1)) allocate(to_sum1(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       if(.not. allocated(to_sum2)) allocate(to_sum2(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       if(.not. allocated(dQ_dp)) allocate(dQ_dp(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,14))
       allocate(dJ0_dp(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       allocate(dgam_dp(nakx, naky,-nzgrid:nzgrid))
       
       call get_mirror_term (change_p, dmirror_dp)
       call get_prll_strm_term (change_p, prll_strm_term)
       call get_omega_star_term (change_p, star_term)
       call get_drift_term (change_p, vm_g_term, vm_phi_term)
       
       
       !!full partial derivative
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          p = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             do q = 1, nvpa
                dG_dp(ikxkyz,q,imu,change_p)= prll_strm_term(ikxkyz,q,imu) + dmirror_dp(ikxkyz,q,imu)&
                     + vm_g_term(ikxkyz,q,imu) +vm_phi_term(ikxkyz,q,imu)+ star_term(ikxkyz,q,imu)
             end do
          end do
       end do


       !!Quasineutrality condition
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          do imu = 1, nmu
             do q = 1, nvpa
                dJ0_dp(imu,ikxkyz) = (aj0v(imu,ikxkyz)- j0_initial(imu,ikxkyz))/del(change_p)
                to_avg1(q,imu,ikxkyz) = dJ0_dp(imu,ikxkyz)*gnew0(q,imu,ikxkyz)
             end do
          end do
       end do

       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          p = iz_idx(kxkyz_lo,ikxkyz)
          call integrate_vmu (to_avg1(:,:,ikxkyz), p, to_sum1(ikxkyz))
       end do

       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          p = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             do q = 1, nvpa
                dgam_dp (ikx,iky,p) = (gamtot(iky,ikx,p) - gamtot0(iky,ikx,p))/del(change_p)
                to_avg2(q,imu,ikxkyz) = spec(is)%zt*dgam_dp(ikx,iky,p)*phi0(q,imu,ikxkyz)
             end do
          end do
       end do
       
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          p = iz_idx(kxkyz_lo,ikxkyz)
          call integrate_vmu (to_avg2(:,:,ikxkyz), p, to_sum2(ikxkyz))
       end do

       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          is = is_idx(kxkyz_lo,ikxkyz)
          dQ_dp(ikxkyz,change_p) = spec(is)%z*spec(is)%dens*(to_sum1(ikxkyz)+to_sum2(ikxkyz))
       end do

       call deallocate_stella

       deallocate(to_sum1)
       deallocate(to_sum2)
       deallocate(to_avg1)
       deallocate(to_avg2)
       deallocate(dJ0_dp)
       deallocate(dgam_dp)
       
    end do

    call deallocate_unperturbed

    call allocate_stella
    
  end subroutine perturb_p

  !! Mirror Term
  subroutine get_mirror_term (change_p, dmirror_dp)
    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid, delzed
    use millerlocal, only: del
    use species, only: spec
    use stella_geometry, only: bmag, gradpar
    use stella_geometry, only: get_dzed
    use kt_grids, only: nalpha
    
    complex, dimension(:), allocatable :: pert_gradpar_dbdz
    complex, dimension(:,:,:), allocatable, intent (out) :: dmirror_dp
    complex, dimension(:), allocatable :: gradpar_dbdz
    real, dimension (:,:), allocatable ::  dbdzed1
    
    integer :: p, ikxkyz, is, imu, q, ia
    integer, intent (in) :: change_p
    
    if(.not. allocated(pert_gradpar_dbdz)) allocate(pert_gradpar_dbdz(-nzgrid:nzgrid))
    if(.not. allocated(dmirror_dp)) allocate(dmirror_dp(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(gradpar_dbdz)) allocate(gradpar_dbdz(-nzgrid:nzgrid))
    if(.not. allocated(dbdzed1)) allocate(dbdzed1(nalpha,-nzgrid:nzgrid))

    ia = 1
    
    call get_dzed (nzgrid, delzed, bmag(ia,:), dbdzed1(ia,:))
    
    do p = -nzgrid, nzgrid
       gradpar_dbdz(p) = gradpar0(p)*dbdzed0(1,p)
       pert_gradpar_dbdz(p)= gradpar(p)*dbdzed1(1,p)
    end do
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       p = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do q = 1, nvpa
             dmirror_dp(ikxkyz,q,imu) = -spec(is)%stm*mu0(imu)*((pert_gradpar_dbdz(p)&
                  -gradpar_dbdz(p))/del(change_p))*dgdvpa0(q,imu,ikxkyz)
          end do
       end do
    end do

    deallocate (pert_gradpar_dbdz)
    deallocate (gradpar_dbdz)
    deallocate (dbdzed1)
    
  end subroutine get_mirror_term
  

  !!Parallel streaming terms 
  subroutine get_prll_strm_term (change_p, prll_strm_term)

    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid
    use millerlocal, only: del
    use species, only: spec
    use stella_geometry, only: gradpar
    use vpamu_grids, only:  maxwell_vpa, maxwell_mu
    use gyro_averages, only: aj1v, aj0v
    
    complex, dimension (:,:,:), allocatable :: vpa_term2_unpert_g, vpa_term2_unpert_phi
    complex, dimension (:,:,:), allocatable :: vpa_term2_unpert_dphi_dz
    real, dimension (:,:,:), allocatable :: dJ0_dz_exp_pert,J0_exp_pert
    complex, dimension (:,:,:), allocatable, intent (out) :: prll_strm_term
    complex, dimension (:,:,:), allocatable :: ps_term1, ps_term2, ps_term3
    complex, dimension(:,:), allocatable :: dgradpar_dp
    integer :: p, ikxkyz, is, imu, q, ikx, iky, ia
    integer, intent (in) :: change_p

    if(.not. allocated(vpa_term2_unpert_g)) allocate(vpa_term2_unpert_g(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc, nvpa,nmu))
    if(.not. allocated(vpa_term2_unpert_phi)) allocate(vpa_term2_unpert_phi(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc, nvpa,nmu))
    if(.not. allocated(vpa_term2_unpert_dphi_dz)) allocate(vpa_term2_unpert_dphi_dz(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc, nvpa,nmu))
    if(.not. allocated(dJ0_dz_exp_pert)) allocate(dJ0_dz_exp_pert(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if(.not. allocated(J0_exp_pert)) allocate(J0_exp_pert(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if(.not. allocated(ps_term1)) allocate(ps_term1(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nvpa,nmu))
    if(.not. allocated(ps_term2)) allocate(ps_term2(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nvpa,nmu))
    if(.not. allocated(ps_term3)) allocate(ps_term3(kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nvpa,nmu))
    if(.not. allocated(prll_strm_term)) allocate(prll_strm_term(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(dgradpar_dp)) allocate(dgradpar_dp(-nzgrid:nzgrid, nvpa))
    
    ia = 1
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       p = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do q = 1, nvpa
             vpa_term2_unpert_g(ikxkyz,q,imu) = spec(is)%stm*dgdz0_v(q,imu,ikxkyz)
             vpa_term2_unpert_phi(ikxkyz,q,imu) = spec(is)%stm*phi0(q,imu,ikxkyz)
             vpa_term2_unpert_dphi_dz(ikxkyz,q,imu) = spec(is)%stm*dphidz0(iky,ikx,p,ia)
             
             dgradpar_dp(p,q)= (gradpar(p) - gradpar0(p))/del(change_p)
             
             dJ0_dz_exp_pert (q,imu,ikxkyz)  = -aj1v(imu,ikxkyz)*maxwell_vpa(q)*maxwell_mu(1,p,imu)
             J0_exp_pert (q,imu,ikxkyz) = aj0v(imu,ikxkyz)*maxwell_vpa(q)*maxwell_mu(1,p,imu)
             
             ps_term1( ikxkyz,q,imu) = vpa0(q)*dgradpar_dp (p,q)*vpa_term2_unpert_g(ikxkyz,q,imu)
             ps_term2(ikxkyz,q,imu) = phi0(q,imu,ikxkyz)*vpa0(q) *(gradpar(p)*dJ0_dz_exp_pert(q,imu,ikxkyz) - gradpar0(p)*dJ0_dz_exp_unpert(q,imu,ikxkyz))/del(change_p)
             ps_term3(ikxkyz,q,imu) = dphidz0(iky,ikx,p,ia)*vpa0(q)*(gradpar(p)*J0_exp_pert(q,imu,ikxkyz)- gradpar0(p) * J0_exp_unpert(q,imu,ikxkyz))/del(change_p)
             
             prll_strm_term(ikxkyz,q,imu) = ps_term1(ikxkyz,q,imu) + ps_term2(ikxkyz,q,imu) + ps_term3(ikxkyz,q,imu)
          end do
       end do
    end do

    deallocate (vpa_term2_unpert_g)
    deallocate (vpa_term2_unpert_phi)
    deallocate (vpa_term2_unpert_dphi_dz)
    deallocate (dJ0_dz_exp_pert)
    deallocate (J0_exp_pert)
    deallocate (ps_term1)
    deallocate (ps_term2)
    deallocate (ps_term3)
    deallocate(dgradpar_dp)
    
  end subroutine get_prll_strm_term
  
  !! get omega* term
  subroutine get_omega_star_term (change_p, star_term)

    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid, ntubes
    use millerlocal, only: del
    use species, only: nspec
    use gyro_averages, only: aj0v
    use dist_fn_arrays, only: wstar
    use stella_layouts, only: vmu_lo, imu_idx, iv_idx
    
    complex, dimension(:,:,:,:,:), allocatable :: dwstar_dp
    complex, dimension(:,:,:), allocatable, intent (out) :: star_term
    complex, dimension(:,:,:,:,:), allocatable :: wstarv1
    complex :: i
    integer :: p, ikxkyz, is, imu, q, ikx, iky, ivmu, ia
    integer, intent (in) :: change_p
    
    ia = 1
    i = (0,1)
    
    if(.not. allocated(dwstar_dp)) allocate(dwstar_dp(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    if(.not. allocated(star_term)) allocate(star_term(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(wstarv1)) allocate(wstarv1(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do p = -nzgrid, nzgrid
          wstarv1(ia,p,q,imu,is)  = wstar(ia,p,ivmu)
       end do
    end do
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       p = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do q = 1, nvpa
             dwstar_dp(ia,p,q,imu,is) = (aj0v(imu,ikxkyz)*wstarv1(ia,p,q,imu,is) - j0_initial(imu,ikxkyz)*wstarv0(ia,p,q,imu,is))/del(change_p)
             star_term(ikxkyz,q,imu) = i*dwstar_dp(ia,p,q,imu,is)*phi0(q,imu,ikxkyz)
          end do
       end do
    end do

    deallocate(dwstar_dp)
    deallocate(wstarv1)
    
  end subroutine get_omega_star_term


  !! Get drift terms
  subroutine get_drift_term (change_p, vm_g_term, vm_phi_term)

    use stella_layouts, only: kxkyz_lo, iz_idx, is_idx, iky_idx, ikx_idx
    use vpamu_grids, only: nvpa, nmu
    use zgrid, only: nzgrid, ntubes
    use millerlocal, only: del
    use species, only: spec, nspec
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use vpamu_grids, only:  maxwell_vpa, maxwell_mu
    use gyro_averages, only: aj0v
    use kt_grids, only: akx, aky
    use stella_layouts, only: vmu_lo, imu_idx, iv_idx
    
    complex, dimension(:,:,:), allocatable :: dvmx_dp_g, dvmy_dp_g, dvmx_dp_phi, dvmy_dp_phi
    complex, dimension(:,:,:), allocatable, intent (out) :: vm_g_term, vm_phi_term
    real, dimension (:,:,:), allocatable :: J0_exp_pert
    complex, dimension(:,:,:,:,:), allocatable :: wdriftx_g1, wdrifty_g1, wdriftx_phi1, wdrifty_phi1
    integer :: p, ikxkyz, is, imu, q, ikx, iky, ia, ivmu
    integer, intent (in) :: change_p

    allocate(wdriftx_g1(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdrifty_g1(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdriftx_phi1(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    allocate(wdrifty_phi1(ntubes,-nzgrid:nzgrid,nvpa,nmu,nspec))
    if(.not. allocated(dvmx_dp_g)) allocate(dvmx_dp_g(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(dvmy_dp_g)) allocate(dvmy_dp_g(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(dvmx_dp_phi)) allocate(dvmx_dp_phi(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(dvmy_dp_phi)) allocate(dvmy_dp_phi(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(vm_g_term)) allocate(vm_g_term(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(vm_phi_term)) allocate(vm_phi_term(kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc,nvpa,nmu))
    if(.not. allocated(J0_exp_pert)) allocate(J0_exp_pert(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

    ia = 1

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       q = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do p = -nzgrid,nzgrid
          wdriftx_g1(ia,p,q,imu,is)  = wdriftx_g(ia,p,ivmu)
          wdrifty_g1(ia,p,q,imu,is) = wdrifty_g(ia,p,ivmu)
          wdriftx_phi1(ia,p,q,imu,is) = wdriftx_phi(ia,p,ivmu)
          wdrifty_phi1(ia,p,q,imu,is) = wdrifty_phi(ia,p,ivmu)
       end do
    end do
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       p = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do q = 1, nvpa
             J0_exp_pert (q,imu,ikxkyz) = aj0v(imu,ikxkyz)*maxwell_vpa(q)*maxwell_mu(1,p,imu)
             
             dvmx_dp_g(ikxkyz,q,imu) = (wdriftx_g1(ia,p,q,imu,is)&
                  - wdriftx_g0(ia,p,q,imu,is))/del(change_p)*akx(ikx)*gnew0(q,imu,ikxkyz)
             dvmy_dp_g(ikxkyz,q,imu) = (wdrifty_g1(ia,p,q,imu,is)&
                  - wdrifty_g0(ia,p,q,imu,is))/del(change_p)*aky(iky)*gnew0(q,imu,ikxkyz)
             dvmx_dp_phi(ikxkyz,q,imu) = spec(is)%zt*(wdriftx_phi1(ia,p,q,imu,is)*J0_exp_pert(q,imu,ikxkyz)- wdriftx_phi0(ia,p,q,imu,is)*J0_exp_unpert(q,imu,ikxkyz))/del(change_p)*akx(ikx)*phi0(q,imu,ikxkyz)
             dvmy_dp_phi(ikxkyz,q,imu) = spec(is)%zt*(wdrifty_phi1(ia,p,q,imu,is)*J0_exp_pert(q,imu,ikxkyz)- wdrifty_phi0(ia,p,q,imu,is)*J0_exp_unpert(q,imu,ikxkyz))/del(change_p)*akx(ikx)*phi0(q,imu,ikxkyz)
             
             vm_g_term(ikxkyz,q,imu) = dvmx_dp_g(ikxkyz,q,imu) + dvmy_dp_g(ikxkyz,q,imu)
             vm_phi_term(ikxkyz,q,imu) = dvmx_dp_phi(ikxkyz,q,imu) + dvmy_dp_phi(ikxkyz,q,imu)
          end do
       end do
    end do

    deallocate (dvmx_dp_g)
    deallocate (dvmy_dp_g)
    deallocate (dvmx_dp_phi)
    deallocate (dvmy_dp_phi)
    deallocate (J0_exp_pert)
    deallocate (wdriftx_g1)
    deallocate (wdrifty_g1)
    deallocate (wdriftx_phi1)
    deallocate (wdrifty_phi1)
    
  end subroutine get_drift_term



  subroutine deallocate_unperturbed
    deallocate(j0_initial)
    deallocate(j1_initial)
    deallocate(dgdvpa0)
    deallocate(darg0dz)
    deallocate(wdriftx_g0)
    deallocate(wdrifty_g0)
    deallocate(wdriftx_phi0)
    deallocate(wdrifty_phi0)
    deallocate(wstarv0)
    deallocate(gradpar0)
    deallocate(bmag0)
    deallocate(gnew0)
    deallocate(arg0)
    deallocate(dbdzed0)
    deallocate(dgdz0_v)
    deallocate(dgdz0)
    deallocate(mu0)
    deallocate(phi0)
    deallocate(dphidz0)
    deallocate(vpa0)
    deallocate(dJ0_dz_exp_unpert)
    deallocate(J0_exp_unpert)
    deallocate(gamtot0)  
  end subroutine deallocate_unperturbed

end module adjoint_p_derviatives
