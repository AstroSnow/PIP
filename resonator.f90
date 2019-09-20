subroutine resonator
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,T0,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option,gra,flag_grav,total_prc,&
       margin,dxc,s_order,mpi_siz,mpi_pos
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use procedures
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision :: temperature(1:ix+1),energy(1:ix+1),rho(1:ix+1)
  double precision ::mask(ix,jx,kx)
  double precision ::v_para(ix,jx,kx),v_perp(ix,jx,kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),B0,P_tot,P_top,P_base, P_grad,rho_ref,lpint,lnint,x0
  double precision lambda_p0,lambda_p1,lambda_h0,lambda_h1,h, P_ref,lambda_pa(1:ix),lambda_na(1:ix)
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr,rng, dg,h0,hn,ttemp,htemp,L0,Rgas,kb_const,mh_const, gra0, grah1,grah2,tmax
  double precision ro0,P0,gran, ro_h_temp, ro_m_temp, p_h_temp, p_m_temp, gra_temp,tn
  double precision, DIMENSION(:),ALLOCATABLE :: t_read,h_read,ln_read,lp_read,gra_read
  double precision, DIMENSION(:),ALLOCATABLE :: ro_h_read,ro_m_read,p_h_read,p_m_read
  integer i,j,k,nxt,ti,xphot,dummy_pass,flag_err


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(debug_direction.eq.1) then
! Loaded VALC temperature profile, add subphotospheric temperature profile, and solve equilibrium
! Gravity rolled to zero on upper boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     print*,'=== flag_PIP should be 1 ==='
     stop
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-200000.0d0 ;end(1)=16000000.0d0          !400.0d0/f_p
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  rng=end(1)-start(1)
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

  if(flag_grav.eq.0) then
     print*,'=== flag_grav should be 1 ==='
     stop
  endif
  gra(:,:,:,:)=0.0d0 
!  gra(:,:,:,1)=-274.0d0 !Normalised further down
  gra0=-274.0d0
  !!!========================================================
  !parameters---------------------
  B0=0.0d0
  rho_ref=1.0d0
  P_ref=1.0d0
  wtr=0.1d0
  Rgas= 8.3145
  kb_const=1.38e-23
  mh_const=1.67e-27

! Boundary Conditions
  v_l=(/1.0d0,P_base,0.0d0,0.0d0,0.0d0,B0*1.0d0,0.0d0,0.0d0/)
  v_r=(/1.0d0,P_top,0.0d0,0.0d0,0.0d0,B0*1.0d0,0.0d0,0.0d0/)
!     Density, pressure, velocity * 3, magnetic field *3

!read temperature profile
  open(unit=101,file="Maltby_prof.dat",status='old')
!  open(unit=101,file="maltby_16384.dat",status='old')
!  open(unit=101,file="VALC_25000.dat",status='old')
  read(101,*) nxt,h0,hn
  L0=(hn-h0)/end(1) !IS THIS THE RIGHT CHOICE?
!  L0=1.0e4
!  gra(:,:,:,1)=gra(:,:,:,1)/L0 !Gravity normalisation
!  gra(:,:,:,1)=1.d0/gm !Gravity normalisation?
  ALLOCATE(t_read(1:nxt))
  ALLOCATE(h_read(1:nxt))
  ALLOCATE(ln_read(1:nxt))
  ALLOCATE(lp_read(1:nxt))
  ALLOCATE(gra_read(1:nxt))
  do i=1,nxt
	read(101,*) htemp,ttemp
	h_read(i)=1.0e3*htemp
	t_read(i)=ttemp
!	print*,h_read(i),t_read(i)
	gra_read(i)=gra0 !ROLE DOWN GRAVITY NEAR TOP
  enddo
  close(101)
  print*,nxt,h0,hn,maxval(h_read)

!!!!!!!!!!!!!!!!!!!!!!!
!Set normalisation parameters based on totals TO DO
!T0=1.0e4 !This is read from the parameters file
!find t(x)=t0
x0=interpol(t_read,h_read,T0)
print*,'T=T0 at height ', x0
!L0=!pressure scale height at T0
!!!!!!!!!!!!!!!!!!!!!!!

grah1=x(ix-1)-x(ix-1)/3.d0
grah2=x(ix-3)
tmax=interpol(h_read,t_read,grah1) !set the max temperature
print*,grah1,grah2
!Calculate lambda for the read values, roll gravity to zero at the top
  do i=1,nxt
	if (h_read(i) .GT. grah2) then
		gra_read(i)=0.0d0
		t_read(i)=tmax
	elseif (h_read(i) .GT. grah1) then
		gra_read(i)=gra_read(i)*(1.0d0+COS(pi*(h_read(i)-grah1)/(grah2-grah1)))/2.0d0
		t_read(i)=tmax
	endif
!	lp_read(i)=gra(1,1,1,1)/2.0d0/t_read(i)/Rgas !more constants?
!	ln_read(i)=gra(1,1,1,1)/t_read(i)/Rgas
	lp_read(i)=gra_read(i)*mh_const/2.0d0/t_read(i)/kb_const !more constants?
	ln_read(i)=gra_read(i)*mh_const/t_read(i)/kb_const
!	print*,ln_read(i)
!	print*,h_read(i),t_read(i),gra_read(i)
  enddo
print*,gra_read(1),t_read(4),Rgas,ln_read(4)

!Interpolate the temperature to the simulation grid
  temperature(:)=t_read(1)
  ti=1
  do i=1,ix
	if (x(i) .GE. h_read(nxt-1)) then
		temperature(i)= t_read(nxt)
		gra(i,:,:,1)=gra_read(nxt)
	elseif (x(i) .LT. 0.0d0) then
		gra(i,:,:,1)=gra0
		temperature(i)=t_read(1) + x(i) * gra(i,1,1,1) / ((1.d0/(gm-1.d0))+1.0d0) * mh_const/kb_const
!		temperature(i)=t_read(1)
	else
		temperature(i)=interpol(h_read,t_read,x(i))
		gra(i,:,:,1)=interpol(h_read,gra_read,x(i))
!		print*,x(i),gra(i,1,1,1)
	endif
  enddo
  temperature(ix+1)=temperature(ix) !Does this exist?

!Set parameters at the photosphere
  xphot=1
  do while (x(xphot) .LT. 0.d0)
	xphot=xphot+1
  enddo 
print*,'xphot = ',xphot
  ro_h(xphot,:,:)=f_n*rho_ref 
  ro_m(xphot,:,:)=f_p*rho_ref 
  p_h(xphot,:,:)=ro_h(xphot,:,:)*temperature(xphot)*kb_const/mh_const
  p_m(xphot,:,:)=ro_m(xphot,:,:)*temperature(xphot)*2.0d0*kb_const/mh_const

  h=x(4)-x(3)
!  print*,h
  lnint=0.0d0
  lpint=0.0d0
!vertical integration for pressure, ideal gas law for density, use mid points
  do i=xphot+1,ix-1 !Integrate up from the photosphere
	lnint=lnint+integrate(h_read,ln_read,x(i),h,3,1)
	lpint=lpint+integrate(h_read,lp_read,x(i),h,3,1)
	p_h(i,:,:)=p_h(xphot,:,:)*exp(lnint)
	p_m(i,:,:)=p_m(xphot,:,:)*exp(lpint)
	ro_h(i,:,:)=p_h(i,:,:)*mh_const/temperature(i)/kb_const
	ro_m(i,:,:)=p_m(i,:,:)*mh_const/temperature(i)/2.0d0/kb_const
  enddo
print*,'Integrated up'
  lnint=0.0d0
  lpint=0.0d0
  do i=xphot-1,2,-1 !Integrate down from the photosphere
	lnint=lnint-0.5*h*(gra(i,1,1,1)*mh_const/temperature(i)/kb_const+gra(i-1,1,1,1)*mh_const/temperature(i-1)/kb_const)
	lpint=lpint-0.5*h*(gra(i,1,1,1)*mh_const/temperature(i)/kb_const/2.d0+gra(i-1,1,1,1)*mh_const/temperature(i-1)/kb_const/2.d0)
	p_h(i,:,:)=p_h(xphot,:,:)*exp(lnint)
	p_m(i,:,:)=p_m(xphot,:,:)*exp(lpint)
	ro_h(i,:,:)=p_h(i,:,:)*mh_const/temperature(i)/kb_const
	ro_m(i,:,:)=p_m(i,:,:)*mh_const/temperature(i)/2.0d0/kb_const
  enddo
print*,'Integrated down'
!ro_m(2,:,:)=ro_m(3,:,:)
ro_m(1,:,:)=ro_m(2,:,:)
!p_m(2,:,:)=p_m(3,:,:)
p_m(1,:,:)=p_m(2,:,:)
!ro_h(2,:,:)=ro_h(3,:,:)
ro_h(1,:,:)=ro_h(2,:,:)
!p_h(2,:,:)=p_h(3,:,:)
p_h(1,:,:)=p_h(2,:,:)
!p_h(ix-1,:,:)=p_h(ix-2,:,:)
p_h(ix,:,:)=p_h(ix-1,:,:)
!p_m(ix-1,:,:)=p_m(ix-2,:,:)
p_m(ix,:,:)=p_m(ix-1,:,:)
!ro_h(ix-1,:,:)=ro_h(ix-2,:,:)
ro_h(ix,:,:)=ro_h(ix-1,:,:)
!ro_m(ix-1,:,:)=ro_m(ix-2,:,:)
ro_m(ix,:,:)=ro_m(ix-1,:,:)

!Need to set ix-1:ix values

print*,'ro_n/(ro_h+ro_m) at 0 = ',ro_h(3,1,1)/(ro_m(3,1,1)+ro_h(3,1,1))
print*,'ro_n/(ro_h+ro_m) at ix = ',ro_h(ix-2,1,1)/(ro_m(ix-2,1,1)+ro_h(ix-2,1,1))

!initialise other arrays
  do k=1,kx;do j=1,jx;do i=1,ix
     vx_h(i,j,k)=0.0d0
     vx_m(i,j,k)=0.0d0
     vy_h(i,j,k)=0.0d0
     vy_m(i,j,k)=0.0d0
     vz_h(i,j,k)=0.0d0
     vz_m(i,j,k)=0.0d0
     b_x(i,j,k)=0.0d0
     b_y(i,j,k)=0.0d0
     b_z(i,j,k)=0.0d0
  enddo;enddo;enddo

!!!!!!!!!!!!!!!!!!!!!!!
!Set normalisation parameters based on totals TO DO
!T0=1.0e4 !This is read from the parameters file
!find T(x)=T0
x0=interpol(temperature(xphot:ix),x(xphot:ix),T0)
L0=-kb_const*T0/gra0/mh_const !pressure scale height at T0, using neutrals
ro0=interpol(x,ro_h(:,1,1),x0)
P0=interpol(x,p_h(:,1,1),x0)*gm
gran=abs(interpol(x,gra(:,1,1,1),x0))*gm
print*,'T=',T0,' at x=', x0,',L0=',L0,', ro0=',ro0, 'P0=',P0
x=x/L0
start=start/L0
end=end/L0
ro_h=ro_h/ro0
ro_m=ro_m/ro0
p_h=p_h/P0
p_m=p_m/P0
gra=gra/gran
!!!!!!!!!!!!!!!!!!!!!!!
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write initial atmosphere
OPEN(UNIT=102, FILE="maltbyatmos_res.txt", ACTION="write", STATUS="replace")
write(102,*) ix, start(1), end(1), L0, T0 , ro0, P0, gran
do i=1,ix 
	write(102,*) x(i), ro_h(i,1,1), ro_m(i,1,1), p_h(i,1,1), p_m(i,1,1), gra(i,1,1,1)
enddo
close(102)

  !!set boundary condition----------------------
  flag_bnd(1)=20
  flag_bnd(2)=10
  flag_bnd(3)=10
  flag_bnd(4)=10
  flag_bnd(5)=10
  flag_bnd(6)=10
  !-------------------------------------------------

print*,'Atmosphere saved. Change debug_direction flag to continue'
stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elseif(debug_direction.eq.10) then
! Loaded saved stable profile with PIP (made using debug_direction=10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     print*,'=== flag_PIP should be 1 ==='
     stop
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

!  dummy_pass=1

!  if (my_rank .GT. 0) then
!	mpi_recv(dummy_pass,1,MPI_INT,my_rank-1,1,MPI_COMM_WORLD)
!  endif

  open(unit=102,file='valcatmos_res.txt',status='old')
!  open(unit=102,file='maltbyatmos_res.txt',status='old')
  read(102,*) nxt, h0, hn, L0, Tn , ro0, P0, gran
  ALLOCATE(h_read(1:nxt))
  ALLOCATE(ro_h_read(1:nxt))
  ALLOCATE(ro_m_read(1:nxt))
  ALLOCATE(p_h_read(1:nxt))
  ALLOCATE(p_m_read(1:nxt))
  ALLOCATE(gra_read(1:nxt))
  do i=1,nxt
!	write(102,*) x(i), ro_h(i,1,1), ro_m(i,1,1), p_h(i,1,1), p_m(i,1,1), gra(i,1,1,1)
	read(102,*) htemp, ro_h_temp, ro_m_temp, p_h_temp, p_m_temp, gra_temp
!	print*,htemp, ro_h_temp, ro_m_temp, p_h_temp, p_m_temp, gra_temp
	h_read(i)=htemp
	ro_h_read(i)=ro_h_temp
	ro_m_read(i)=ro_m_temp
	p_h_read(i)=p_h_temp
	p_m_read(i)=p_m_temp
	gra_read(i)=gra_temp
  enddo

  close(102)

!  if (my_rank .LT. total_prc) then
!	mpi_send(dummy_pass,1,MPI_INT,my_rank+1,1,MPI_COMM_WORLD)
!  endif

  print*,'data read on processor',my_rank

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=h0 ;end(1)=hn          !400.0d0/f_p
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  rng=end(1)-start(1)
  call set_coordinate(start,end)
  !---------------------------------------
  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
       ix,x,dx,dxc,s_order,0.0d0-start(1),20.d0,dx(1)*0.5d0,1.001d0,dx(1)*10.0d0, &
       mpi_pos(1),margin(1),0,flag_err)
!print*,dx
x=x+h0 !fix for having negative x after NUG
!  call set_coordinate_NUG_T(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
!       ix,x,dx,dxc,s_order,start(1),x(4),dx(1),10.0d0,dx(1)*4.0d0, &
!       mpi_pos(1),margin(1),0,flag_err,5.0d0*p_h_read(:)/ro_h_read(:)/3.0d0)
!Set domain values
print*,'initial height',x(1),h0,start(1)
print*,'synthesised atmosphere height',hn,'calculated height',maxval(x)

!stop

!Set domain values
  gra(:,:,:,:)=0.0d0 
  do i=1,ix
!	print*,'x',x(i),h_read(i)
	ro_h(i,:,:)=interpol(h_read,ro_h_read,x(i))
!	print*,'ro_h',ro_h(i,1,1),ro_h_read(i)
	ro_m(i,:,:)=interpol(h_read,ro_m_read,x(i))
!	print*,'ro_m',ro_m(i,1,1),ro_m_read(i)
	p_h(i,:,:)=interpol(h_read,p_h_read,x(i))
!	print*,'p_h',p_h(i,1,1),p_h_read(i)
	p_m(i,:,:)=interpol(h_read,p_m_read,x(i))
!	print*,'p_m',p_m(i,1,1),p_m_read(i)
	gra(i,:,:,1)=interpol(h_read,gra_read,x(i))
!	print*,'gra',gra(i,1,1,1),gra_read(i)
!stop
  enddo

if (x(1) .lt. h0) then
	ro_h(1,:,:)=ro_h(2,:,:)
	ro_m(1,:,:)=ro_m(2,:,:)
	p_h(1,:,:)=p_h(2,:,:)
	p_m(1,:,:)=p_m(2,:,:)
	gra(1,:,:,1)=gra(2,:,:,1)
endif

print*,'finished rank',my_rank

!initialise other arrays
  do k=1,kx;do j=1,jx;do i=1,ix
     vx_h(i,j,k)=0.0d0
     vx_m(i,j,k)=0.0d0
     vy_h(i,j,k)=0.0d0
     vy_m(i,j,k)=0.0d0
     vz_h(i,j,k)=0.0d0
     vz_m(i,j,k)=0.0d0
     b_x(i,j,k)=0.0d0
     b_y(i,j,k)=0.0d0
     b_z(i,j,k)=0.0d0
  enddo;enddo;enddo


  !!set boundary condition----------------------
  flag_bnd(1)=20
  flag_bnd(2)=10
  flag_bnd(3)=10
  flag_bnd(4)=10
  flag_bnd(5)=10
  flag_bnd(6)=10
  !-------------------------------------------------

!print*,x(1:4)
!print*,ro_h(1:4,1,1)
!print*,ro_m(1:4,1,1)
!print*,p_h(1:4,1,1)
!print*,p_m(1:4,1,1)
!print*,gra(1:4,1,1,1)

!stop

endif


  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=0.4d0
     dtout=tend/10.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine resonator
