subroutine blob
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin
  use parameters,only:pi
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
  double precision Atwood,ro_l,ro_u,w_lay
  double precision A(jx,3),b(jx),P_y(jx)
  ! jack's vars 
  double precision debug
  double precision :: blobrad,blobht,pos
  double precision, allocatable:: p_tot(:),ro_tot(:),temp(:)
  double precision t0,t1,trw,trp
  double precision pr_val,ro_jump
  double precision xe,xs,ye,ys,ze,zs,dx0,dy0,dz0
  integer orix,oriy,oriz
  integer i,j,k,l,lmin,sx
  ! end my vars
  
  !print *, "welcome to blob.f90, please wipe your feet."
  
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
  ! stolen straight from RT.f90 (IN TERMS OF SCALE HEIGHTS)
  !i wanna make the space bigger
  start(1)=-50.0d0 ;end(1)=50.0d0
  start(2)=0.0d0 ;end(2)=100.0d0 !!NN
  start(3)=-50.0d0 ;end(3)=50.0d0
  call set_coordinate(start,end)
  !---------------------------------------

  
  ! ==================== v from old code - rewrite?
  ! set grid space
  ! dx0 not needed as dx is the replacement which is set in model.f90; these are now arrays, as the gird space may not be uniform, but currently all the same (probs dont change this)
  ! maybe replace all dx0 with dx(1) dx(l)
  dx0=(end(1)-start(1))/(ix-2*margin(1)+1)
  dy0=(end(2)-start(2))/(jx-2*margin(2)+1)
  dz0=(end(3)-start(3))/(kx-2*margin(3)+1)
!! AH: xe and xs had not been defined. I rewrote them as start and end to match the lines of code above. without this the gize was not being defined!


  ! set origin of grid space
  ! "at which grid point is xs associated with?"
  orix=margin(1)+1
  oriy=margin(2)+1
  oriz=margin(3)+1
  ! ==================== ^ from old code - rewrite?
   
  
  !!default boundary condition----------------------
  ! need to make sure these are right; no idea which is which or what I want tbh
  ! want reflection at bottom (and top?) and free on all other sides.
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=2
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=2
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------
  
  
  
  !!!=================
  !!!=========SHIT GETS MESSY BELOW HERE!====================
  !!!=================
  
  
  
  !!!========================================================
  !density of lower fluid is unity
  !think I maybe want to change this bit a bit or smthn? 
  ro_jump=9.0d0
  blobrad=10.0d0
  blobht=170.0d0
!! AH: I have changed these to fit the smaller box I imposed
  

      
      
  !settin up gravity? idk what scl_height is
!! AH: scl_height is the pressure scale height (see pdf I made for you some time ago). This should be set to 1 for almost all simulations (but some test problems use a different value to compare with famous solutions).
  flag_grav=1
  gra(:,:,:,:)=0.0d0
!  scl_height=0.50    !!!!! AH: This is not needed
! maybe remove all instances of scl_height b/c it should automatically be 1 from the temp
  if(scl_height.eq.0.0) scl_height=1.0 
  gra(:,oriy:(jx-margin(2)-1),:,2)=-1.0/(gm*scl_height)
  ! above: gravity for the main middle bit, except it goes too high
  !below: bottom boundary, gravity switches
  gra(:,1:margin(2),:,2)=1.0/(gm*scl_height)
  !now do the same for the top boundary
  gra(:,(jx-margin(2)):jx,:,2)=1.0/(gm*scl_height)

! AH: for the sym boundary, you have to have symmetric gravity as well 

!*spread(spread( &
!       0.5d0*(tanh((y+0.7)/0.15)+1)*(1.0-tanh((y-0.7)/0.15)),1,ix),3,kx)
!!!!!AH: This is not needed (localises gravity in space). This works nicely for the RT (as it reduces boundary problems), but for your problem you want to   
  
  
  !!!====================
  ! can i just use jx instead of sx since y is always going to be "up?"

  !! AH: Yes you can 
  !!!====================
  sx = jx
  
  !???????????? I have a feeling this may not be correct anymore
  allocate(p_tot(sx),ro_tot(sx),temp(sx))
       
  !set up temp; DOES THIS WORK? i dont get this new lingo being used with spread and y and etc
  t0 = 1.0d0
  t1 = 100.0d0
!!!! AH: 1000 is too high (unphysically?) and will make the simualtion harder
!150 maybs? coz taky says so
!  trp = 20.0d0
  trp=5.0d0 !!NN
  trw = 2.0d0
!! AH: This is the width of the transistion region, so though you need a few grid points to resolve it, using 5 pressure scale heights is too much.
! maybe 1? ...2?
  temp(:)=t0 + (t1 - t0) * 0.5d0*(tanh((y(:)-trp)/trw)+1.d0)

  !print *, temp
  !set up pressure? from the temp
  !oriy might need reassigning or smth?
! AH: You have to first assign pr_val as pr_val=1/gm

  pr_val=1.d0/gm

  do  l=1,sx
    if (l .eq. 1) then
    !andy set it up like this just because it worked (its a hack, not a beautiful solution)
      p_tot(l)=pr_val
    else
      lmin=l-1
      p_tot(l)=p_tot(lmin) - dy0 * p_tot(lmin)/temp(lmin)
    endif  
!     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,dy0,p_tot(l)
!    print *,y(l),p_tot(l)
  enddo
!  stop

  !set up the density from the pressure?
  do l=1,sx
    ro_tot(l)=p_tot(l) *gm/temp(l)
  enddo
  
  
  !now make the ARRAYS!
  ! is this ok? ??? guessing based on the kind of thing in RT.f90
!  do l=1,jx   
    P_m=spread(spread(p_tot*f_p_p,1,ix),3,kx)
!!! AH: I rewrote this to use spread. The function is spread(VARIABLE,DIRECTION,SIZE)
    ro_m=spread(spread(ro_tot*f_p ,1,ix),3,kx)         
    P_h=spread(spread(p_tot*f_p_n ,1,ix),3,kx)
    ro_h=spread(spread(ro_tot*f_n ,1,ix),3,kx)
!print *, 'hi'
!if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *, p_m(1,:,1)
  
!stop

  
  
    ! convinced the next line is total crap - do i need x(i) y(j) and shit? no?
!! AH: calculating 'pos' gives the geometric postion from the centre of the blob. So this gives a simple way of writing a localised blob. Why do you think this is wrong? 

! see below from andy 18/9/14 @ 1am; ± may be wrong but yeah
!  ro_m(i,j,k)=ro_m(i,j,k)-ro_jump*0.5d0*(tanh((pos-blob_rad)/blob_tr)-1.0d0)
! but first define blob_tr (0.1d0 ?)

      do k=1,kx
        do j=1,jx
          do i=1,ix   
            pos = sqrt( x(i)**2 + (y(j) - blobht)**2 + z(k)**2 )
            if (pos .le. blobrad) then
            ro_m(i,j,k)=ro_m(i,j,k)+ro_jump*ro_m(i,j,k)
            ro_h(i,j,k)=ro_h(i,j,k)+ro_jump*ro_h(i,j,k)
          endif
        enddo
      enddo
    enddo
  
  
  ! can i / would it be better to do like:
  !ro_h=f_n*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0)*0.5d0,1,ix),3,kx)
  !ro_m=f_p*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0)*0.5d0,1,ix),3,kx)
  
!! AH as this ahs no x and z dependence it would create a massive sheet of huge extent in the vertical dependence (also only has a start y range). That is a huge amount of mass you are going to make fall.
  
  !====== MAGNETIC FIELD ======
  do k=1,kx
    do j=1,jx
      do i=1,ix     
        B_x(i,j,k)=0.0d0
!sqrt(2.d0*P_tot(280)/beta)
!! AH: again I have set this to 0. Obviously we want a bit of magnetic field action, but for debugging purposes KISS (Keep It Simple Stupid)
        B_z(i,j,k)=0
      enddo
    enddo
  enddo
  
    
  !!!=================
  !!!=================BELOW HERE IS COOL AGAIN I THINK MAYBE============
  !!!=================
  
  
  !start stationary
  vx_h=0.0d0;vy_h=0.0d0;vz_h=0.0d0
  vx_m=0.0d0;vy_m=0.0d0;vz_m=0.0d0


  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------


  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=40.0d0
     dtout=tend/40.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------


end subroutine blob
