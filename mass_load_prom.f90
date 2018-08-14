subroutine mass_load_prom
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,flag_divb &
       ,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,mpi_pos,mpi_siz,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
  use boundary_rot,only:PIPbnd    ! here all the boundary routines are called
  use mpi_rot,only:my_mpi_barrier
  implicit none
!  include "mpi.h"
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision pr_val
  real*4 :: lenx,lenz,znull,zaxis
  integer :: num,nx,nz,i,j,opt,div_x,div_y,start_x(2),start_y(2)
  character*40 :: filename
  double precision, allocatable:: b_mag(:,:)
  real*8 :: bx(1012,1012),by(1012,1012),bz(1012,1012),ro(1012,1012),pr(1012,1012)
  double precision, allocatable:: p_tot(:),ro_tot(:),temp(:),y_1(:)
  double precision t0,t1,trw,trp,b_norm


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
  ! load the file
      num=2
      if(num.eq.0) stop 'Quit'

      write (filename,fmt='("relaxed_initial_",i1,".dat")') num
      if (my_rank.eq.1) print *,'Reading data from '//filename
      open(unit=1,file=filename,form='unformatted',status='old')

      nx=1012; nz=1012           !  read array dimensions                                                                      
      if (my_rank.eq.1) print 10,nx ,nz
      if (my_rank.eq.0) print 20,ix,jx
 10   format(" nx =",i4,", nz =",i4)
 20   format(" ix=",i4,", jx=",i4)
 
  ! read in the grid
                                                                                                               
      if (my_rank.eq.2) print *, 'Start read'

      read (1) ((pr(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.3) print *, 'loaded pr'

      read (1) ((ro(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.3) print *, 'loaded ro'

      read (1) ((bx(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.3) print *, 'loaded b_x'

      read (1) ((by(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.0) print *, 'loaded B_y'

      read (1) ((bz(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.1) print *, 'loaded B_z'

      close(1)
      if (my_rank.eq.2) print *,'OK'


     ! NEED TO SET UP B FOR DIFFERENT CORES!!
!    if (my_rank.eq.3) print *, ro(:,20)
    div_x=ix-2*margin(1)
    div_y=jx-2*margin(2)
    if (my_rank.eq.3) print *, div_x,div_y, my_rank,'L84'




    start_x(1)=1+mpi_pos(1)*div_x
    start_x(2)=(mpi_pos(1)+1)*div_x+2*margin(1)

    start_y(1)=1+mpi_pos(2)*div_y
    start_y(2)=(mpi_pos(2)+1)*div_y+2*margin(2)

   if (my_rank.eq.0) print *, ix,jx,kx,'L95'
   if (my_rank.eq.1) print *, start_x,start_y,my_rank,'L96'
    ! deal with Margin!!

   !-----------------------------------------------------
  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-2.5d0 ;end(1)=2.5d0
  start(2)=0.0d0 ;end(2)=5.0d0
  start(3)=0.0d0 ;end(3)=10.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=2
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=2
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=4
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=2
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !density of lower fluid is unity
  nz=1000
  nx=1000

  flag_grav=1
  gra(:,:,:,:)=0.0d0
  scl_height=1.d0
  gra(:,:,:,2)=-1.0d0/(gm*scl_height)

  if (mpi_pos(2).eq.0) then
  gra(:,1:margin(2),:,2)=1.0d0/(gm*scl_height)
  endif
  if (mpi_pos(2).eq.(mpi_siz(2)-1)) then
  gra(:,jx-margin(2)+1:jx,:,2)=1.0d0/(gm*scl_height)
  endif


    P_m=DBLE(spread(pr(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))
!spread(spread(p_tot(start_y(1):start_y(2)) &
!      *f_p_p,1,ix),3,kx)
!    ro_m= 1.d0*f_p 
     ro_m= DBLE(spread(ro(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))
!spread(spread(ro_tot(start_y(1):start_y(2)) & 
!      *f_p ,1,ix),3,kx)         

    P_h= pr_val*f_p_n 
!     P_h= spread(spread(p_tot(start_y(1):start_y(2)) & 
!      *f_p_n ,1,ix),3,kx)
    ro_h= 1.d0*f_n 
!     ro_h= spread(spread(ro_tot(start_y(1):start_y(2)) &
!      *f_n ,1,ix),3,kx)



!   b_norm=sqrt(p_tot(328)/2.d0/beta)/bz(499,328)
!sqrt(p_tot(164)/2.d0/beta)/bz(499,164)
   if (my_rank.eq.0) print *, b_norm

  vx_h=0.0d0;vy_h=0.0d0;vz_h=0.0d0
  vx_m=0.0d0;vy_m=0.0d0;vz_m=0.0d0
!  B_x=0.0d0;B_y=0.0d0;B_z=0.0d0


    B_x=DBLE(spread(bx(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))

    if (my_rank.eq.3) print *, 'done bx', my_rank

    B_y=DBLE(spread(by(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))

    if (my_rank.eq.2) print *, 'done by',my_rank

    B_z=DBLE(spread(bz(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))
    if (my_rank.eq.1) print *, 'done bz',my_rank

!    print *, B_z(200,1:12,1)

!   B_x=0.0d0;B_y=0.0d0;B_z=0.0d0  

  !!!=================
  ! Add mass

!  print *, 'start_ro',my_rank

  print *, ro(503,328),f_p,jx,kx
  call my_mpi_barrier
!  STOP
   
  ro_m=ro_m
!+ro(503,328)*10.d0*f_p &
!       *spread(0.5d0*(1.d0-dtanh((abs(spread(x,2,jx))-0.05d0)/0.05d0)) &
!     *0.5d0*(1.d0-dtanh((abs(spread(y-1.45d0,1,ix))-0.35d0 )/0.05d0)),3,kx)

!     *spread(spread(1.d0+0.01d0*cos(z*4.d0*pi),1,ix),2,jx)  &

!  ro_h=ro_h+100.d0*f_n*spread(0.5d0*(1.d0-dtanh((abs(spread(x,2,jx))-0.1d0 )/0.1d0)) &
!     *0.5d0*(1.d0-dtanh((abs(spread(y-0.9d0,1,ix))-0.4d0 )/0.05d0)),3,kx)*ro_tot(328)
  !!!========================================================
  print *, 'finished ro',my_rank


  call my_mpi_barrier
  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  call my_mpi_barrier
  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=40.0d0
     dtout=tend/40.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine mass_load_prom
