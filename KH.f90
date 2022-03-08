subroutine KH
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
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
  double precision Atwood,ro_l,ro_u,vx_l,vx_u,w_lay,b0,theta
  double precision A(jx,3),b(jx),P_y(jx)
  integer i,j,k
  integer, allocatable:: new(:), old(:)
  integer size, seed(2), gseed(2), hiseed(2), zseed(2)
  real harvest(ix*jx*kx)
    data seed /123456789, 987654321/
    data hiseed /-1, -1/
    data zseed /0, 0/
   call random_seed(SIZE=size)

   ALLOCATE (new(size))
   ALLOCATE (old(size))
   CALL RANDOM_SEED(GET=old(1:size))
   new = old*(my_rank+1)
   CALL RANDOM_SEED(PUT=new(1:size))
   call random_number(HARVEST=harvest)

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
  !!set lower and upper coordinate
  start(1)=-0.5d0 ;end(1)=0.5d0
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-8.0d0 ;end(3)=8.0d0
  call set_coordinate(start,end)
  !---------------------------------------

  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=1
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=1
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=2
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=2
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !density of lower fluid is unity
  ro_l=1.0d0
  ro_u=1.0e2!10.0d0
  vx_l=ro_u/(ro_u+ro_l)*dsqrt(1.d0/10.d0)
!  if(flag_mhd.eq.1) then
!  vx_l=vx_l*sqrt(1.d0+2.d0/(gm*beta))
!  endif
  vx_u=-vx_l*ro_l/ro_u
  w_lay=0.003d0

  scl_height=0.50
  if(scl_height.eq.0.0) scl_height=1.0

!  ro_h=f_n*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  ro_m=f_p*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  vx_h=spread(spread(vx_l+(vx_u-vx_l)*(tanh(y/w_lay)+1.0)* &
!       0.5d0,1,ix),3,kx)

  ro_h=f_n*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0) &
       *0.5d0,1,ix),3,kx)
  ro_m=f_p*ro_h/f_n
  vx_h=abs(ro_l*vx_l)*spread(spread(tanh(y/w_lay),1,ix),3,kx)/ro_h
  vx_m=vx_h

  theta=2.d0*pi*0.d0/360.d0
  if(flag_mhd.eq.1) then
     b0=sqrt(2.0d0/(gm*beta))
     !B0=0.d0 !Hydrodynamic case for testing
     B_z=B0*cos(theta)
     B_x=B0*sin(theta)
     B_y=0.0d0
  endif


  P_h(:,:,:)=1.0/gm*f_p_n
  P_m(:,:,:)=1.0/gm*f_p_p

!  P_h(:,:,:)=2.5d0*f_p_n
!  P_m(:,:,:)=2.5d0*f_p_p

  vy_h=0.0d0;vz_h=0.0d0
  vy_m=0.0d0;vz_m=0.0d0

!  print *, 'done', my_rank

  do k=1,kx
     do j=1,jx
        do i=1,ix
        vy_h(i,j,k)=0.001d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!*dcos(x(i)*2.d0*pi)*exp(-2.d0*pi*abs(y(j))) &
!                 *(-0.5d0*dtanh((dabs(y(j))-1.5d0)/0.15d0)+0.5d0 )
!*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        vy_m(i,j,k)=0.001d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!dcos(x(i)*2.d0*pi)*exp(-2.d0*pi*abs(y(j))) &
!(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        enddo
     enddo
  enddo

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=20.0d0
     dtout=tend/40.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine KH



