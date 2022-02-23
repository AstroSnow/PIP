subroutine RT
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
  double precision Atwood,ro_l,ro_u,w_lay,kw
  double precision A(jx,3),b(jx),P_y(jx)
  double precision, allocatable:: p_tot(:),ro_tot(:),temp(:),b0(:)
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
  start(1)=-0.001d0 ;end(1)=0.001d0
  start(2)=-0.0065d0 ;end(2)=0.0035d0
  start(3)=-0.001d0 ;end(3)=0.001d0
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
  Atwood=0.5
  ro_l=1.0d0
!  ro_u=ro_l+16.0d0
  ro_u=1.5d0
!(1.0+Atwood)/(1.0-Atwood)
  w_lay=0.00003d0
  kw = 4.d0*pi
  
  if(flag_grav.eq.0) then
     print*,'=== flag_grav should be 1 ==='
     stop
  endif
  gra(:,:,:,:)=0.0d0
  scl_height=0.50
  if(scl_height.eq.0.0) scl_height=1.0  
  gra(:,:,:,2)=-0.1d0
 
  ro_h=f_n*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0)*0.5d0,1,ix),3,kx)
  ro_m=f_p*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0)*0.5d0,1,ix),3,kx)

   allocate(p_tot(jx))     

    do  j=1,jx
    p_tot(j)=2.0d0-ro_h(10,j,1)*y(j)*0.1d0
! Note grav is set to 0,1 explicitly here
    enddo


  P_h=f_n*spread(spread(p_tot,1,ix),3,kx)
  P_m=f_p*spread(spread(p_tot,1,ix),3,kx)


  vx_h=0.0d0;vy_h=0.0d0;vz_h=0.0d0
  vx_m=0.0d0;vy_m=0.0d0;vz_m=0.0d0


  !  vy_h=0.0001d0*spread(spread(cos(4*pi*x),2,jx),3,kx)* &
!  vy_h=0.0001d0*spread(spread(cos(kw*x),2,jx),3,kx)* &  
!       spread(spread(exp(-abs(y/w_lay)),1,ix),3,kx)

  do k=1,kx
     do j=1,jx
        do i=1,ix
        vy_h(i,j,k)=0.02d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        vy_m(i,j,k)=0.02d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
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
     tend=4.0d0
     dtout=tend/4.0/5.d0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine RT
