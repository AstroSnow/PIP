!! For Cylindrical coordinates only
!! Ref: CANS md_diskflare, Hayashi+1996 ApJL
subroutine disk_flare
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,dxc,dyc,dzc,n_fraction,gra,flag_grav,scl_height,margin, &
       flag_cyl,debug_parameter,s_order,mpi_siz,mpi_pos
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision :: sseps = 0.2d0, npol = 3.d0, &
       eth = 2.d-3, emg = 3.d-5, tec0 = 1.d0, roc0 = 1.d-3
  double precision :: gpot(ix,jx,kx),gra0(ix,jx,kx)
  double precision :: ro_disk(ix,jx,kx),pr_disk(ix,jx,kx),vz_disk(ix,jx,kx)
  double precision :: ro_corona(ix,jx,kx),pr_corona(ix,jx,kx)
  double precision :: ro_tot(ix,jx,kx),p_tot(ix,jx,kx)
  double precision :: ss,sr,coef,b0,aa
  double precision :: psi0,te,tmp
  integer :: i,j,k,ibnd,ierr,flag_err_NUG
  double precision :: x_uni0,x_uni1,dx_uni,x_fact,dx_max
  double precision :: y_uni0,y_uni1,dy_uni,y_fact,dy_max  

  
  if(flag_cyl.ne.1) then
     if(my_rank.eq.0) then
        print*,'Cylindrical coordinates should be used: STOP'
!        call MPI_ABORT(mpi_comm_world,ierr)
        stop
     endif
  endif
  if(debug_parameter.ge.0.5d0) then
     if(my_rank.eq.0) then
        print*,'debug_parameter should be in the range of [0,0.5]'
!        call MPI_ABORT(mpi_comm_world,ierr)
        stop
     endif
  endif
  aa = debug_parameter
  
  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !! Note: If you use uniform and non-uniform meshes together,
  !!       then you have to call set_coordinate first.
  !!       Call set_coordinate_NUG afer that.
  !!set lower and upper coordinate
  start(1)=0.0d0 ;end(1)=15.0d0
  start(2)=0.0d0 ;end(2)=15.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  ! !---------------------------------------

  !Set coordinate (Non-uniform grid)--------------------------


  !! X-direction
  x_uni0 = 0.d0
  x_uni1 = 2.d0
  dx_uni = 1.d-2
  x_fact = 1.015d0
!  x_fact = 1.0d0
  dx_max = 10.d0*dx_uni
  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1), &
       ix,x,dx,dxc,s_order,x_uni0,x_uni1,dx_uni,x_fact,dx_max,mpi_pos(1),margin(1), &
       0,flag_err_NUG)
  if(flag_err_NUG.ne.0) then
     print*,'Grid # in x-direction is too small. STOP'
     stop
     !     call MPI_ABORT(mpi_comm_world,ierr)
  endif
  
  !! Y-direction
  y_uni0 = 0.d0
  y_uni1 = 1.5d0
  dy_uni = 5.d-3
  y_fact = 1.015d0
!  y_fact = 1.0d0
  dy_max = 10.d0*dy_uni
  call set_coordinate_NUG(mpi_siz(2)*(jx-2*margin(2))+2*margin(2), &
       jx,y,dy,dyc,s_order,y_uni0,y_uni1,dy_uni,y_fact,dy_max,mpi_pos(2),margin(2), &
       0,flag_err_NUG)
  if(flag_err_NUG.ne.0) then
     print*,'Grid # in y-direction is too small. STOP'
     stop
     !     call MPI_ABORT(mpi_comm_world,ierr)
  endif


  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=5  !! symmetric
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10 !! free
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=6  !! symmetric (at the disk mid-plane)
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10 !! free
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  
  if(flag_grav.eq.0) then
     print*,'=== flag_grav should be 1 ==='
!     call MPI_ABORT(mpi_comm_world,ierr)
     stop
  endif
  gpot(:,:,:) = 0.d0
  do j=1,jx ; do i=1,ix
     !     ss = x(i)
     ss = sqrt(x(i)**2+y(j)**2)
     if(ss.gt.sseps) then
        gpot(i,j,:) = -1.d0/ss
     else
        if(ss.gt.0.5d0*sseps) then
           gpot(i,j,:) = - (1.d0/sseps - (ss-sseps)/sseps**2)
        else
           gpot(i,j,:) = - 1.5d0/sseps
        endif
     endif
  enddo;enddo

  gra(:,:,:,:)=0.0d0  !! initialization
  do i=2,ix-1
     gra0(i,:,:) = - (gpot(i+1,:,:)-gpot(i-1,:,:))/(2.d0*dx(i))
  enddo
  gra(:,:,:,1) = gra0
  do j=2,jx-1
     gra0(:,j,:) = - (gpot(:,j+1,:)-gpot(:,j-1,:))/(2.d0*dy(j))
  enddo
  gra(:,:,:,2) = gra0

  !! BC for gra
  if(minval(x).lt.0.d0) then
     ibnd = 1 + margin(1)
     do i=1,margin(1)
        gra(ibnd-i,:,:,1) = -gra(ibnd-1+i,:,:,1)
     enddo
  endif
  if(minval(y).lt.0.d0) then
     ibnd = 1 + margin(2)
     do j=1,margin(2)
        gra(:,ibnd-j,:,2) = -gra(:,ibnd-1+j,:,2)
     enddo
  endif

  b0 = sqrt(emg)

  !! Disk
  coef = 1.d0/(npol+1.d0)/eth
  vz_disk = 0.d0;ro_disk=0.d0;pr_disk=0.d0 ! initialization
  do k=1,kx;do j=1,jx;do i=1,ix
     tmp = coef*( -(1.d0+gpot(i,j,k)) &
             + 1.d0/(2.d0*(1.d0-aa))*(1.d0-x(i)**(2.d0*aa-2.d0) )) + 1.d0
     if(tmp.gt.0.d0) then
        ro_disk(i,j,k) = tmp**npol
        vz_disk(i,j,k) = x(i)**(aa-1.d0)
     endif
  enddo;enddo;enddo
  where(ro_disk.gt.0.d0) pr_disk = eth * ro_disk**(1.d0+1.d0/npol)
!  where(ro_disk.gt.0.d0) te_disk = pr_disk/ro_disk * gm
  do k=1,jx;do j=1,jx;do i=1,ix
     ss = sqrt(x(i)**2+y(j)**2)
     B_x(i,j,:) = b0*3.d0*x(i)*y(j)/ss**2.5
     B_y(i,j,:) = b0*(2.d0*y(j)**2-x(i)**2)/ss**2.5
     B_z(i,j,:) = 0.d0     
  enddo;enddo;enddo
  
  ! psi0 = -1.d0 + 1.d0/2.d0/(1.d0-aa) + (npol+1.d0)*eth
  ! ro_disk = 0.d0; pr_disk = 0.d0; vz_disk = 0.d0 !! initialization
  ! do j=1,jx ; do i=1,ix
  !    sr = x(i)
  !    ss = sqrt(x(i)**2+y(j)**2)
  !    if(sr.gt.sseps) then
  !       te = (psi0+1.d0/ss-1.d0/2.d0/(1.d0-aa)*sr**(2.d0*aa-2.d0))/(npol+1.d0)*gm
  !       if(te.gt.0.d0) then
  !          ro_disk(i,j,:) = (te/gm/eth)**npol
  !          pr_disk(i,j,:) = te * ro_disk(i,j,:) / gm
  !          vz_disk(i,j,:) = sr**(aa-1.d0)
  !       endif
  !    endif

  !    ! B_x(i,j,:) = 0.d0
  !    ! B_y(i,j,:) = 0.d0
  !    ! B_z(i,j,:) = 0.d0
  !    B_x(i,j,:) = b0*3.d0*x(i)*y(j)/ss**2.5
  !    B_y(i,j,:) = b0*(2.d0*y(j)**2-x(i)**2)/ss**2.5
  !    B_z(i,j,:) = 0.d0

  ! enddo;enddo

  !! Corona
  ro_corona(:,:,:) = roc0*exp(-gm/tec0*(gpot(:,:,:)+1.d0))
  pr_corona(:,:,:) = ro_corona * tec0/gm

  !! Total
  ! do k=1,kx;do j=1,jx;do i=1,ix
  !    if(ro_disk(i,j,k).gt.0.d0) then
  !       ro_tot(i,j,k) = ro_disk(i,j,k)
  !       p_tot(i,j,k)  = pr_disk(i,j,k)
  !    else
  !       ro_tot(i,j,k) = ro_corona(i,j,k)
  !       p_tot(i,j,k)  = pr_corona(i,j,k)
  !    end if
  ! enddo;enddo;enddo

  ro_tot = ro_disk + ro_corona
  p_tot  = pr_disk + pr_corona

  ! ro_tot = ro_corona
  ! p_tot = pr_corona  
  
  ro_h = f_p_n * ro_tot
  ro_m = f_p_p * ro_tot
  p_h = f_p_n * p_tot
  p_m = f_p_p * p_tot  

  vx_h = 0.d0; vx_m = 0.d0
  vy_h = 0.d0; vy_m = 0.d0
  vz_h = vz_disk ; vz_m = vz_disk
  !vz_h = 0.d0 ; vz_m = 0.d0

  
  !!!========================================================
  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=10.0d0
     !tend=2.0d0     
     dtout=0.2d0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine DISK_FLARE
