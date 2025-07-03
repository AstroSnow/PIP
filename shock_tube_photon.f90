subroutine shock_tube_photon
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin,T0,n0,Nexcite,&
       f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref,colrat,gm_ion,gm_rec,expinttab,&
       flag_rad,radrat,rad_temp,gm_rec_rad,gm_ion_rad,flag_IR,n_levels,n0,T0,nexcite0
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
  use PIP_rot, only:get_col_ion_coeff,expintread,get_radrat_fixed,set_NLTE_equilibrium
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx)
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision Atwood,ro_l,ro_u,vx_l,vx_u,w_lay,b0,theta
  double precision A(jx,3),b(jx),P_y(jx)
  double precision ::mask(ix,jx,kx)  
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr
  integer i,j,k,ii
  integer, allocatable:: new(:), old(:)
  integer size, seed(2), gseed(2), hiseed(2), zseed(2)
  real harvest(ix*jx*kx)
  double precision:: T0down,T0up,n0up,n0down
  double precision:: nnup,nndown,pnup,pndown,ppup,ppdown,b0up,b0down,ptot
  double precision:: f_nup,f_pup,f_ndown,f_pdown, radius, r0,v0
  double precision:: f_p_nup,f_p_pup,f_p_ndown,f_p_pdown
  double precision::Nexciteup(n_levels+1),Nexcitedown(n_levels+1),Eion(6)
  double precision,parameter::kbhat=1.38064852,mehat=9.10938356,hhat=6.62607004
  double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!random number generator
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set the reference values for temperature [K] and electron number density [m^-3]

!small jump
!T0down=6555.84250200805d0
!T0up=5500.d0
!n0up=7.5e16
!n0down=n0up*10.d0
!bigger jump
!T0down=7319.689479843136d0
!T0up=5500.d0
!n0up=7.0d16
!n0down=n0up*30.d0
!Tube values
T0down=6000.d0
T0up=T0down!6000.d0
n0up=1.0d16
n0down=n0up*1.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Check that the initial conditions are consistent with settings
if (abs(T0up-T0) .gt. 1.0e-6) then 
    print*,'T_norm in settings neq T0up. Normalisation wont work'
    print*,'T0= ',T0
    print*,'T0up= ',T0up
    T0up=T0
    T0down=T0
endif
!Check that the initial conditions are consistent with settings
if (abs(n0up-n0) .gt. 1.0e-6) then 
    print*,'n0 in settings neq n0up. Normalisation wont work'
    print*,'n0= ',n0
    print*,'n0up= ',n0up
    stop
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Flag to make sure I have the right ionisation routine
if(flag_IR .ne. 4) then
	print*,'set flag_IR=4 for this routine'
	!stop
endif
  
if (my_rank.eq.0) print*,'Calculating LTE excitation state'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Allocate arrays 
!This is only used for the n_level stuff
!Should be within a flag or deleted for MHD to save memory
allocate(nexcite0(n_levels+1))
allocate(Nexcite(ix,jx,kx,n_levels+1)) !Allocate the fractional array
allocate(Colrat(ix,jx,kx,n_levels+1,n_levels+1))
call expintread
if (flag_rad .ge. 2) allocate(radrat(ix,jx,kx,n_levels+1,n_levels+1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ionisation energies
Eion=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
Eion=Eion/13.6*2.18e-18 !Convert to joules (to be dimensionally correct)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the LTE level populations using Saha-Boltzman
Nexciteup(n_levels+1)=n0up
do i=1,n_levels
Nexciteup(i)=(2.d0/n0up/2.d0*(2.d0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.d0/2.d0)*exp(-Eion(i)/kboltz/T0up))
enddo
Nexciteup(1:n_levels)=n0up/Nexciteup(1:n_levels)
Nexciteup=Nexciteup/n0up
Nexciteup=Nexciteup/sum(Nexciteup(:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call set_NLTE_equilibrium(T0up,n0up,nexciteup,1.0d-3,10000)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Nexcitedown(n_levels+1)=n0down
do i=1,n_levels
    Nexcitedown(i)=(2.d0/n0down/2.d0*(2.d0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.d0/2.d0)*exp(-Eion(i)/kboltz/T0down))
enddo
Nexcitedown(1:n_levels)=n0down/Nexcitedown(1:n_levels)
Nexcitedown=Nexcitedown/n0down
Nexcitedown=Nexcitedown/sum(Nexcitedown(:))


!Calculate the xion fraction and pressure ratio for the LTE state
print*,'overwritting LTE state'
Nexciteup(1:n_levels)=Nexciteup(1:n_levels)*100.0
Nexcitedown(1:n_levels)=Nexcitedown(1:n_levels)*100.0


f_nup=sum(Nexciteup(1:n_levels))/(sum(Nexciteup(1:n_levels))+Nexciteup(n_levels+1)) !neutral fraction
f_pup=1.d0-f_nup   !Ion fraction
nnup=n0up/f_pup-n0up
pnup=nnup*T0up*3.d0/5.d0!f_nup/(f_nup+2.0d0*f_pup)
ppup=n0up*T0up*6.d0/5.d0
f_p_nup=pnup/n0up!(pnup+ppup)
f_p_pup=ppup/n0up!(pnup+ppup)

print*,nexciteup
print*,nnup/n0up,n0up/n0up
!stop

f_p_ini=n0up!f_pup
f_p_p_ini=ppup!f_p_pup
n0fac=1.d0!f_pup

f_ndown=sum(Nexcitedown(1:n_levels))/(sum(Nexcitedown(1:n_levels))+Nexcitedown(n_levels+1))
f_pdown=1.d0-f_ndown!Nexcite(1,1,1,6)
nndown=n0down/f_pdown-n0down
pndown=nndown*T0down*3.d0/5.d0!f_nup/(f_nup+2.0d0*f_pup)
ppdown=n0down*T0down*6.d0/5.d0
f_p_ndown=pndown/n0up!(pnup+ppup)
f_p_pdown=ppdown/n0up!(pnup+ppup)

f_p=f_pup
f_n=f_nup
f_p_p=f_p_pup
f_p_n=f_p_nup
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=0.0d0 ;end(1)=10000.0d0
  start(2)=0.0d0 ;end(2)=2.0d0
  start(3)=-10.0d0 ;end(3)=10.0d0
  call set_coordinate(start,end)
  !---------------------------------------

!Boundary conditions. 
!flag_bnd(1) is left
!flag_bnd(2) is right
!flag_bnd(3) is top
!flag_bnd(4) is bottom
!Values:  1 (periodic), 2(symmetric).......
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=3
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=10
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

  !!!========================================================
  !density of lower fluid is unity
  ro_l=nnup+n0up !total (plasma + neutral) density
  ro_u=nndown+n0down



  vx_l=ro_u/(ro_u+ro_l)*dsqrt(1.d0/10.d0) !Velocity (unused here)
!  if(flag_mhd.eq.1) then
!  vx_l=vx_l*sqrt(1.d0+2.d0/(gm*beta))
!  endif
  vx_u=-vx_l*ro_l/ro_u !Shear velocity
  w_lay=0.003d0
  w_lay=0.1d0

!Put the up and down fractions as densities
  Nexciteup(1:n_levels)=Nexciteup(1:n_levels)*(nnup+n0up)
  Nexcitedown(1:n_levels)=Nexcitedown(1:n_levels)*(nndown+n0down)
  Nexciteup(n_levels+1)=Nexciteup(n_levels+1)*(nnup+n0up)
  Nexcitedown(n_levels+1)=Nexcitedown(n_levels+1)*(nndown+n0down)

  theta=2.d0*pi*0.d0/360.d0

!Calculate the initial magentic field strength based on the plasma-beta given in settings file  
  if (flag_pip.eq.0) then
     ptot=ppdown
     b0=sqrt(2.0d0/(ptot/ppup/5.0*3.0*beta))
  else
     ptot=ppdown+pndown
     b0=sqrt(2.0d0/(ptot/ppup/5.0*3.0*beta))
     print*,'check total pressure for two fluid case'
  endif
  b0down=b0


print*,ppup,ppdown
!print*,b0up,b0down

!Set velocity arrays to zero
  vy_h=0.0d0;vz_h=0.0d0
  vy_m=0.0d0;vz_m=0.0d0


  r0=0.5 !Radius
  v0=0.01 !Reference velocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!Set Shock jumps
  mask=spread(spread(pi*x,2,jx),3,kx)
  phi_p=0.0d0
  theta_p=pi/2.0d0
  tmp=0.0

!print*,'Setting B0=0'
!  B0=0.0d0

  v_l=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,B0,0.0d0/)
  v_r=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,-B0,0.0d0/)

print*,v_l(7),v_r(7)
print*,'B0',B0
!stop

  where(mask<0)
     ro_h=f_n*v_l(1)
     p_h=f_p_n*v_l(2)
     ro_m=f_p*v_l(1)
     p_m=f_p_p*v_l(2)
     vx_h=v_l(3)
     vy_h=v_l(4)
     vz_h=v_l(5)
     vx_m=v_l(3)
     vy_m=v_l(4)
     vz_m=v_l(5)
     b_para=v_l(6)
     b_perp=v_l(7)
  elsewhere
     ro_h=f_n*v_r(1)
     p_h=f_p_n*v_r(2)
     ro_m=f_p*v_r(1)
     p_m=f_p_p*v_r(2)
     vx_h=v_r(3)
     vy_h=v_r(4)
     vz_h=v_r(5)
     vx_m=v_r(3)
     vy_m=v_r(4)
     vz_m=v_r(5)
     b_para=v_r(6)
     b_perp=v_r(7)
  end where
!  wtr=0.001d0
  wtr=dx(1)*10.d0

  do k=1,kx;do j=1,jx;do i=1,ix
     ro_h(i,j,k)=nnup!f_n*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     ro_m(i,j,k)=n0up!f_p*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     p_h(i,j,k)=pnup!f_p_n*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     p_m(i,j,k)=ppup!f_p_p*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vx_h(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vx_m(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vy_h(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vy_m(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vz_h(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vz_m(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     b_para(i,j,k)=(v_l(6)+(v_r(6)-v_l(6))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     b_perp(i,j,k)=(v_l(7)+(v_r(7)-v_l(7))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     Nexcite(i,j,k,1:n_levels+1)=Nexciteup(1:n_levels+1)
  enddo;enddo;enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  b_theta=b_perp*sin(tmp)
  b_phi  =b_perp*cos(tmp)
  b_x=b_para*sin(theta_p)*cos(phi_p)+ &
       b_theta*cos(theta_p)*cos(phi_p)-b_phi*sin(phi_p)
  b_y=b_para*sin(theta_p)*sin(phi_p)+ &
       b_theta*cos(theta_p)*sin(phi_p)+b_phi*cos(phi_p)
  b_z=b_para*cos(theta_p)-b_theta*sin(theta_p) 

!Normalisation to P_m=1/gamma, ro_m=1, such that the temperature of the plasma (T=gamma P/ ro) is 1
P_m=P_m/5.0*3.0/ppup
P_h=P_h/5.0*3.0/ppup
ro_m=ro_m/n0up!*T0up
ro_h=ro_h/n0up!*T0up
Nexcite=Nexcite/n0up

!Set the fraction for the ionisation
!f_p_ini=f_pup
!f_p_p_ini=f_p_pup

!Make sure that the level populations are properly defined
if (flag_IR .eq. 4) then
Nexcite(:,:,:,n_levels+1)=ro_m
do i=1,n_levels
    Nexcite(:,:,:,i)=Nexcite(:,:,:,i)*ro_h/sum(Nexcite(:,:,:,1:n_levels),dim=4)
enddo
!set an initial reference value for normalisation
do i=1,n_levels+1
    nexcite0(i)=nexciteup(i)/n0up
enddo
endif

if ((flag_MHD .eq. 1) .and. (flag_PIP .ne. 1)) then
    print*,'MHD model'
    ro_m=ro_m!+ro_h
    print*,'Density range = ',maxval(ro_m),minval(ro_m)
    P_m=P_m!+P_h
    print*,'Pressure range = ',maxval(P_m),minval(P_m)
    flag_IR=0 !make sure IR is turned off
    flag_rad=0 !set radiation to zero, possibly not wanted
endif

!print*,ro_m(20,1,1),ro_h(20,1,1)
print*,p_m(20,1,1),p_h(20,1,1)
!print*,b_y(1:20,1,1)
!stop
!print*,p_m(20,1,1)/ro_m(20,1,1)*5.0/6.0*T0
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

end subroutine shock_tube_photon



