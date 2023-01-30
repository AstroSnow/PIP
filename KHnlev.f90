subroutine KHnlev
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin,T0,n0,Nexcite,&
       f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref,colrat,gm_ion,gm_rec,expinttab,&
       flag_rad,radrat,rad_temp,gm_rec_rad,gm_ion_rad,flag_IR
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
  integer i,j,k
  integer, allocatable:: new(:), old(:)
  integer size, seed(2), gseed(2), hiseed(2), zseed(2)
  real harvest(ix*jx*kx)
  double precision:: T0down,T0up,n0up,n0down
  double precision:: nnup,nndown,pnup,pndown,ppup,ppdown
  double precision:: f_nup,f_pup,f_ndown,f_pdown
  double precision:: f_p_nup,f_p_pup,f_p_ndown,f_p_pdown
  double precision::Nexciteup(6),Nexcitedown(6),Eion(6)
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
T0down=6555.84250200805d0
T0up=5500.d0
!T0down=12284.274572709573
n0up=7.5e16
n0down=n0up*10.d0
!n0down=n0up*50.d0

if(flag_IR .ne. 4) then
	print*,'set flag_IR=4 for this routine'
	stop
endif
  
if (my_rank.eq.0) print*,'Calculating LTE excitation state'

!Allocate arrays
allocate(Nexcite(ix,jx,kx,6)) !Allocate the fractional array
allocate(Colrat(ix,jx,kx,6,6))
call expintread
if (flag_rad .ge. 2) allocate(radrat(ix,jx,kx,6,6))

Eion=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
Eion=Eion/13.6*2.18e-18 !Convert to joules (to be dimensionally correct)

Nexciteup(6)=n0up
Nexciteup(1)=(2.0/n0up/2.d0*(2.0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(1)/kboltz/T0up))
Nexciteup(2)=(2.0/n0up/8.d0*(2.0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(2)/kboltz/T0up))
Nexciteup(3)=(2.0/n0up/18.d0*(2.0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(3)/kboltz/T0up))
Nexciteup(4)=(2.0/n0up/32.d0*(2.0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(4)/kboltz/T0up))
Nexciteup(5)=(2.0/n0up/50.d0*(2.0*pi*mehat*kbhat*T0up/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(5)/kboltz/T0up))
Nexciteup(1:5)=n0up/Nexciteup(1:5)
Nexciteup=Nexciteup/n0up
Nexciteup=Nexciteup/sum(Nexciteup(:))

Nexcitedown(6)=n0down
Nexcitedown(1)=(2.0/n0down/2.d0*(2.0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(1)/kboltz/T0down))
Nexcitedown(2)=(2.0/n0down/8.d0*(2.0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(2)/kboltz/T0down))
Nexcitedown(3)=(2.0/n0down/18.d0*(2.0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(3)/kboltz/T0down))
Nexcitedown(4)=(2.0/n0down/32.d0*(2.0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(4)/kboltz/T0down))
Nexcitedown(5)=(2.0/n0down/50.d0*(2.0*pi*mehat*kbhat*T0down/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(5)/kboltz/T0down))
Nexcitedown(1:5)=n0down/Nexcitedown(1:5)
Nexcitedown=Nexcitedown/n0down
Nexcitedown=Nexcitedown/sum(Nexcitedown(:))

f_nup=sum(Nexciteup(1:5))/sum(Nexciteup(1:6))
f_pup=1.d0-f_nup!Nexcite(1,1,1,6)
nnup=n0up/f_pup-n0up
pnup=nnup*T0up*3.d0/5.d0!f_nup/(f_nup+2.0d0*f_pup)
ppup=n0up*T0up*6.d0/5.d0
f_p_nup=pnup/n0up!(pnup+ppup)
f_p_pup=ppup/n0up!(pnup+ppup)

f_p_ini=n0up!f_pup
f_p_p_ini=ppup!f_p_pup
n0fac=1.d0!f_pup

f_ndown=sum(Nexcitedown(1:5))/sum(Nexcitedown(1:6))
f_pdown=1.d0-f_ndown!Nexcite(1,1,1,6)
nndown=n0down/f_pdown-n0down
pndown=nndown*T0down*3.d0/5.d0!f_nup/(f_nup+2.0d0*f_pup)
ppdown=n0down*T0down*6.d0/5.d0
f_p_ndown=pndown/n0up!(pnup+ppup)
f_p_pdown=ppdown/n0up!(pnup+ppup)

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
  ro_l=nnup+n0up
  ro_u=nndown+n0down
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

!  ro_h=f_n*spread(spread(ro_l+(ro_u-ro_l)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  ro_m=f_p*ro_h/f_n
!  ro_h=spread(spread(f_nup+(f_ndown-f_nup)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  ro_m=spread(spread(f_pup+(f_pdown-f_pup)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  P_m=spread(spread(f_p_pup+(f_p_pdown-f_p_pup)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
!  P_h=spread(spread(f_p_nup+(f_p_ndown-f_p_nup)*(tanh(y/w_lay)+1.0) &
!       *0.5d0,1,ix),3,kx)
  ro_h=spread(spread(nnup+(nndown-nnup)*(tanh(y/w_lay)+1.0) &
       *0.5d0,1,ix),3,kx)
  ro_m=spread(spread(n0up+(n0down-n0up)*(tanh(y/w_lay)+1.0) &
       *0.5d0,1,ix),3,kx)
  P_m=spread(spread(ppup+(ppdown-ppup)*(tanh(y/w_lay)+1.0) &
       *0.5d0,1,ix),3,kx)
  P_h=spread(spread(pnup+(pndown-pnup)*(tanh(y/w_lay)+1.0) &
       *0.5d0,1,ix),3,kx)
  do i=1,5
  	Nexcite(:,:,:,i)=ro_h*spread(spread(Nexciteup(i)+(Nexcitedown(i)-Nexciteup(i))*(tanh(y/w_lay)+1.0) &
  	     *0.5d0,1,ix),3,kx)
  enddo
 Nexcite(:,:,:,6)=ro_m*spread(spread(Nexciteup(6)+(Nexcitedown(6)-Nexciteup(6))*(tanh(y/w_lay)+1.0) &
             *0.5d0,1,ix),3,kx)
!  ro_m=f_p*ro_h/f_n
  vx_h=abs((nnup+n0up)*vx_l)*spread(spread(tanh(y/w_lay),1,ix),3,kx)/(nnup+n0up)
  vx_m=vx_h

print*,sum(Nexcite(1,1,1,1:5)),ro_h(1,1,1)
print*,Nexcite(1,1,1,6),ro_m(1,1,1)

stop

  theta=2.d0*pi*0.d0/360.d0
  if(flag_mhd.eq.1) then
     b0=sqrt(2.0d0/(gm*beta))
     B_z=B0*cos(theta)
     B_x=B0*sin(theta)
     B_y=0.0d0
  endif


!  P_h(:,:,:)=1.0/gm*f_p_n
!  P_m(:,:,:)=1.0/gm*f_p_p

!  P_h(:,:,:)=2.5d0*f_p_n
!  P_m(:,:,:)=2.5d0*f_p_p

  vy_h=0.0d0;vz_h=0.0d0
  vy_m=0.0d0;vz_m=0.0d0

!  print *, 'done', my_rank

  do k=1,kx
     do j=1,jx
        do i=1,ix
        vy_h(i,j,k)=0.01d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
!*dcos(x(i)*2.d0*pi)*exp(-2.d0*pi*abs(y(j))) &
!                 *(-0.5d0*dtanh((dabs(y(j))-1.5d0)/0.15d0)+0.5d0 )
!*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
        vy_m(i,j,k)=0.01d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0)
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

end subroutine KHnlev



