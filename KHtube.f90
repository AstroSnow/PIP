subroutine KHtube
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin,T0,n0,Nexcite,&
       f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref,colrat,gm_ion,gm_rec,expinttab,&
       flag_rad,radrat,rad_temp,gm_rec_rad,gm_ion_rad,flag_IR,n_levels
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
  use PIP_rot, only:get_col_ion_coeff,expintread,get_radrat_fixed
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

!Set the reference values
!Test values
!T0down=7319.689479843136d0
!T0up=5500.d0
!n0up=7.5e16
!n0down=n0up*10.d0
!small jump
!T0down=6555.84250200805d0
!T0up=5500.d0
!n0up=7.5e16
!n0down=n0up*10.d0
!bigger jump
T0down=7319.689479843136d0
T0up=5500.d0
n0up=7.0e16
n0down=n0up*30.d0
!Tube values
!T0down=100000.d0
!T0up=6000.d0
!n0up=7.5e16
!n0down=n0up*1.d0

if(flag_IR .ne. 4) then
	print*,'set flag_IR=4 for this routine'
	!stop
endif
  
if (my_rank.eq.0) print*,'Calculating LTE excitation state'

!Allocate arrays
allocate(Nexcite(ix,jx,kx,n_levels+1)) !Allocate the fractional array
allocate(Colrat(ix,jx,kx,n_levels+1,n_levels+1))
call expintread
if (flag_rad .ge. 2) allocate(radrat(ix,jx,kx,n_levels+1,n_levels+1))

Eion=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
Eion=Eion/13.6*2.18e-18 !Convert to joules (to be dimensionally correct)

Nexciteup(n_levels+1)=n0up
do i=1,n_levels
Nexciteup(i)=(2.d0/n0up/2.d0*(2.d0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.d0/2.d0)*exp(-Eion(i)/kboltz/T0up))
enddo
Nexciteup(1:n_levels)=n0up/Nexciteup(1:n_levels)
Nexciteup=Nexciteup/n0up
Nexciteup=Nexciteup/sum(Nexciteup(:))

Nexcitedown(n_levels+1)=n0down
do i=1,n_levels
    Nexcitedown(i)=(2.d0/n0down/2.d0*(2.d0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.d0/2.d0)*exp(-Eion(i)/kboltz/T0down))
enddo
Nexcitedown(1:n_levels)=n0down/Nexcitedown(1:n_levels)
Nexcitedown=Nexcitedown/n0down
Nexcitedown=Nexcitedown/sum(Nexcitedown(:))

f_nup=sum(Nexciteup(1:n_levels))/(sum(Nexciteup(1:n_levels))+Nexciteup(n_levels+1))
f_pup=1.d0-f_nup!Nexcite(1,1,1,6)
nnup=n0up/f_pup-n0up
pnup=nnup*T0up*3.d0/5.d0!f_nup/(f_nup+2.0d0*f_pup)
ppup=n0up*T0up*6.d0/5.d0
f_p_nup=pnup/n0up!(pnup+ppup)
f_p_pup=ppup/n0up!(pnup+ppup)

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

  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-2.0d0 ;end(1)=2.0d0
  start(2)=0.0d0 ;end(2)=2.0d0
  start(3)=-10.0d0 ;end(3)=10.0d0
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
  ro_l=nnup+n0up
  ro_u=nndown+n0down

  vx_l=ro_u/(ro_u+ro_l)*dsqrt(1.d0/10.d0)
!  if(flag_mhd.eq.1) then
!  vx_l=vx_l*sqrt(1.d0+2.d0/(gm*beta))
!  endif
  vx_u=-vx_l*ro_l/ro_u
  w_lay=0.003d0
  w_lay=0.1d0

!Put the up and down fractions as densities
  Nexciteup(1:n_levels)=Nexciteup(1:n_levels)*(nnup+n0up)
  Nexcitedown(1:n_levels)=Nexcitedown(1:n_levels)*(nndown+n0down)
  Nexciteup(n_levels+1)=Nexciteup(n_levels+1)*(nnup+n0up)
  Nexcitedown(n_levels+1)=Nexcitedown(n_levels+1)*(nndown+n0down)

  theta=2.d0*pi*0.d0/360.d0
  
  if (flag_pip.eq.0) then
     ptot=ppdown
     b0=sqrt(2.0d0/(ptot/ppup/5.0*3.0*beta))
  else
     ptot=ppdown+pndown
     b0=sqrt(2.0d0/(ptot/ppup/5.0*3.0*beta))
     print*,'check total pressure for two fluid case'
  endif
!  if(flag_mhd.eq.1) then
!     !b0=sqrt(2.0d0/(gm*beta))
!     b0=sqrt(2.0d0/(ptot/ppup/5.0*3.0*beta))
!     !B_z=B0*cos(theta)
!     !B_x=B0*sin(theta)
!     !B_y=0.0d0
!  endif
  b0down=b0
  !b0up=dsqrt(2.0*(0.5*b0down**2+ppdown/5.0*3.0/ppup - 3.0/5.0))! dsqrt(ppdown/ppup)

print*,ppup,ppdown
print*,b0up,b0down
!stop

!  P_h(:,:,:)=1.0/gm*f_p_n
!  P_m(:,:,:)=1.0/gm*f_p_p

!  P_h(:,:,:)=2.5d0*f_p_n
!  P_m(:,:,:)=2.5d0*f_p_p

  vy_h=0.0d0;vz_h=0.0d0
  vy_m=0.0d0;vz_m=0.0d0

!  print *, 'done', my_rank

  r0=0.5
  v0=0.01
  do k=1,kx
     do j=1,jx
        do i=1,ix
        
        radius=dsqrt(y(j)**2+x(i)**2)
        
        b_x(i,j,k)=0.d0
        b_y(i,j,k)=0.d0
        vy_m(i,j,k)=0.d0
        vz_m(i,j,k)=0.d0
        vy_h(i,j,k)=0.d0
        vz_h(i,j,k)=0.d0
        if (radius .GT. 2.0*r0) then
          ro_m(i,j,k)=n0down
          ro_h(i,j,k)=nndown
          P_m(i,j,k)=ppdown
          P_h(i,j,k)=pndown
          b_z(i,j,k)=b0down
          vx_m(i,j,k)=v0
          vx_h(i,j,k)=v0
          do ii=1,n_levels
                  Nexcite(i,j,k,ii)=Nexcitedown(ii)
          enddo
          Nexcite(i,j,k,n_levels+1)=Nexcitedown(n_levels+1)

        else
          ro_m(i,j,k)=n0up+(n0down-n0up)*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          ro_h(i,j,k)=nnup+(nndown-nnup)*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          P_m(i,j,k) =ppup+(ppdown-ppup)*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          P_h(i,j,k) =pnup+(pndown-pnup)*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          !print*,ro_h(i,j,k),ro_m(i,j,k),P_h(i,j,k),P_m(i,j,k)
          !b_z(i,j,k) =b0up+(b0down-b0up)*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          if (flag_pip .eq. 0) then
              b_z(i,j,k)=dsqrt(2.0*(0.5*b0down**2+ppdown/5.0*3.0/ppup - 3.0/5.0*P_m(i,j,k)/ppup))
          else
              b_z(i,j,k)=dsqrt(2.0*(0.5*b0down**2+(ppdown+pndown)/5.0*3.0/ppup - 3.0/5.0*(P_m(i,j,k)+P_h(i,j,k))/ppup))
          endif

          !vx_m(i,j,k)=v0*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          !vx_h(i,j,k)=v0*(tanh((radius-r0)/w_lay)+1.0)*0.5d0

          do ii=1,n_levels
  	          Nexcite(i,j,k,ii)=Nexciteup(ii)+(Nexcitedown(ii)-Nexciteup(ii))*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
          enddo
          Nexcite(i,j,k,n_levels+1)=Nexciteup(n_levels+1)+(Nexcitedown(n_levels+1)&
          -Nexciteup(n_levels+1))*(tanh((radius-r0)/w_lay)+1.0)*0.5d0
        endif

        vx_m(i,j,k)=0.5d0*(dsin(3.14d0*(z(k))/10.d0))*v0*(1.0-tanh((radius-r0)/w_lay))
        vx_h(i,j,k)=0.5d0*(dsin(3.14d0*(z(k))/10.d0))*v0*(1.0-tanh((radius-r0)/w_lay))
        
!vy_h(i,j,k)=0.01d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!*dcos(x(i)*2.d0*pi)*exp(-2.d0*pi*abs(y(j))) &
!                 *(-0.5d0*dtanh((dabs(y(j))-1.5d0)/0.15d0)+0.5d0 )
!*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        !vy_m(i,j,k)=0.01d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!dcos(x(i)*2.d0*pi)*exp(-2.d0*pi*abs(y(j))) &
!(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        enddo
     enddo
  enddo

!Normalisation
P_m=P_m/5.0*3.0/ppup
P_h=P_h/5.0*3.0/ppup
ro_m=ro_m/n0up!*T0up
ro_h=ro_h/n0up!*T0up
Nexcite=Nexcite/n0up


!Make sure that the level populations are properly defined
if (flag_IR .eq. 4) then
Nexcite(:,:,:,n_levels+1)=ro_m
do i=1,n_levels
    Nexcite(:,:,:,i)=Nexcite(:,:,:,i)*ro_h/sum(Nexcite(:,:,:,1:n_levels),dim=4)
enddo
endif

if ((flag_MHD .eq. 1) .and. (flag_PIP .ne. 1)) then
    print*,'MHD model'
    ro_m=ro_m!+ro_h
    print*,'Density range = ',maxval(ro_m),minval(ro_m)
    P_m=P_m!+P_h
    print*,'Pressure range = ',maxval(P_m),minval(P_m)
    flag_IR=0
    flag_rad=0
endif

!print*,nnup+(nndown-nnup)*(tanh(y(1)/w_lay)+1.0)*0.5d0
!print*,sum(Nexcite(1,1,1,1:5)),ro_h(1,1,1)
!print*,Nexcite(1,1,1,6),ro_m(1,1,1)

!print*,ro_m+ro_h
!print*,(P_m+P_h)/(ro_m+ro_h)*5.0/3.0!/(f_p_p_ini/f_p_ini*5.0/6.0)
!print*,nndown+n0down,nnup+n0up
!print*,f_p_pup,f_p_pdown/(pndown+ppdown)
!print*,(P_m+P_h)/(ro_m+ro_h)*5.0/6.0!/(f_p_p_ini/f_p_ini*5.0/6.0)
!print*,(P_m)/(ro_m)*5.0/6.0*(5.0/3.0*ppup/n0up)
!print*,ppup/n0up*5.0/6.0,ppdown/n0down*5.0/6.0
!print*,f_p_pup/f_pup,f_p_pdown/f_pdown

!print*,'NEED TO DO SOME NORMALISATION'

!stop

!call get_col_ion_coeff(T0+0.d0*U_m(:,:,:,1),n0+0.d0*U_m(:,:,:,1),Gm_ion,Gm_rec)
!call get_radrat_fixed(rad_temp,T0+0.d0*U_m(:,:,:,1),T0+0.d0*U_m(:,:,:,1),n0+0.d0*U_m(:,:,:,1),Gm_ion_rad,Gm_rec_rad)
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

end subroutine KHtube



