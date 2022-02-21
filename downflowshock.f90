subroutine downflowshock
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option,gra,flag_grav,total_prc, &
       dxc,s_order,margin,mpi_siz, mpi_pos,T0
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
  double precision lambda_p0,lambda_p1,lambda_h0,lambda_h1,h, P_ref,lambda_pa(1:ix+1),lambda_na(1:ix+1)
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr,rng, dg,h0,hn,ttemp,htemp,L0,Rgas,kb_const,mh_const, gra0,grah1,grah2
  integer i,j,k,nxt,ti,xphot,flag_err
  double precision, DIMENSION(:),ALLOCATABLE :: t_read,h_read,ln_read,lp_read,gra_read
  double precision, DIMENSION(:),ALLOCATABLE :: ro_h_read,ro_m_read,p_h_read,p_m_read
  double precision ro0,P0,gran, ro_h_temp, ro_m_temp, p_h_temp, p_m_temp, gra_temp,tn,grah3,grah4,t_temp,h_temp
  double precision, DIMENSION(ix) :: T_model, h_model

!============================================================
! Set different possible initial profiles
!============================================================
if(debug_direction.eq.1) then
! Loaded temperature profile with PIP, gravity fixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     print*,'=== flag_PIP should be 1 ==='
     f_n=0.d0
     f_p=1.0d0
     f_p_n=0.d0
     f_p_p=1.d0
!     stop
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-2.0e7 ;end(1)=1.0e8
  start(2)=-1.0d6 ;end(2)=1.0d6
  start(3)=-1.0d0 ;end(3)=1.0d0
  rng=end(1)-start(1)
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=1
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=1
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=10
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  if(flag_grav.eq.0) then
     print*,'=== flag_grav should be 1 ==='
     stop
  endif
  gra(:,:,:,:)=0.0d0 
!  gra(:,:,:,1)=-274.0d0 !Normalised further down
  gra(:,:,:,1)=-274.0d0
  gra0=-274.0d0

print*,'setting gravity'
!Upper limits for the gravity reduction
grah1=9.0e7!x(ix-1)-x(ix-1)/3.d0
grah2=10.0e7!x(ix-3)
!Lower limits for the gravity reduction
grah3=-1.0e6
grah4=-2.0e6
  do i=1,ix
	if (x(i) .GT. grah2) then
		gra(i,:,:,1)=0.0d0
	elseif (x(i) .GT. grah1) then
		gra(i,:,:,1)=gra(i,:,:,1)*(1.0d0+COS(pi*(x(i)-grah1)/(grah2-grah1)))/2.0d0
	elseif (x(i) .LT. grah4) then
		gra(i,:,:,1)=0.0d0
	elseif (x(i) .LT. grah3) then
!		gra(i,:,:,1)=gra(i,:,:,1)*(1.0d0-COS(pi*(x(i)-grah1)/(grah2-grah1)))/2.0d0
		gra(i,:,:,1)=gra(i,:,:,1)*(1.0d0+cos(pi*(x(i)-grah3)/(grah4-grah3)))/2.0d0
	endif
  enddo

  !!!========================================================
  !parameters---------------------
  B0=0.0d0
  rho_ref=1.0d0
  P_ref=1.0d0
  wtr=0.1d0
  Rgas= 8.3145d0
  kb_const=1.38e-23
  mh_const=1.67e-27

! Boundary Conditions
  v_l=(/1.0d0,P_base,0.0d0,0.0d0,0.0d0,B0*1.0d0,0.0d0,0.0d0/)
  v_r=(/1.0d0,P_top,0.0d0,0.0d0,0.0d0,B0*1.0d0,0.0d0,0.0d0/)
!     Density, pressure, velocity * 3, magnetic field *3

!print*,'setting temperature profile'
!  open(unit=101,file="sunquakes_height_temp_file_dimensionless.txt",status='old')
  do i=1,ix
!	read(101,*) h_temp,t_temp
!	h_model(i)=h_temp*165264.d0
!	t_model(i)=t_temp*6583.d0
    t_model(i)=1.0e4+(1.0e6-1.0e4)*0.5d0*(tanh(x(i)/5.d0/(x(2)-x(1)))+1.d0)
!	print*,h_read(i),t_read(i)
!	gra_read(i)=gra0 !ROLE DOWN GRAVITY NEAR TOP
  enddo
!  close(101)

!Calculate temperature and lambda
print*,h_model,t_model
  do i=1,ix
!	temperature(i)=T0!**2.0d0*mh_const/gm/kb_const
!	temperature(i)=interpol(h_model,t_model,x(i))
    temperature(i)=t_model(i)
    if (temperature(i) .lt. T0) then 
        temperature(i) = T0
    endif
	print*,x(i),temperature(i)
	lambda_pa(i)=gra(i,1,1,1)*mh_const/2.0d0/temperature(i)/kb_const
	lambda_na(i)=gra(i,1,1,1)*mh_const/temperature(i)/kb_const

  enddo

  ro_h(1,:,:)=f_n*rho_ref 
  ro_m(1,:,:)=f_p*rho_ref 
  p_h(1,:,:)=ro_h(1,:,:)*temperature(1)*kb_const/mh_const
  p_m(1,:,:)=ro_m(1,:,:)*temperature(1)*2.0d0*kb_const/mh_const

  h=x(4)-x(3)
!  print*,h
  lnint=0.0d0
  lpint=0.0d0
print*,'vertical integration'
!vertical integration for pressure, ideal gas law for density, use mid points
  do i=2,ix-1
	lnint=lnint+integrate(x,lambda_na,x(i),h,1,1)
	lpint=lpint+integrate(x,lambda_pa,x(i),h,1,1)

	p_h(i,:,:)=p_h(1,:,:)*exp(lnint)
	p_m(i,:,:)=p_m(1,:,:)*exp(lpint)
	ro_h(i,:,:)=p_h(i,:,:)*mh_const/temperature(i)/kb_const
	ro_m(i,:,:)=p_m(i,:,:)*mh_const/temperature(i)/2.0d0/kb_const
	print*,x(i),p_m(i,1,1),p_m(i,1,1)/ro_m(i,1,1)/2.d0*mh_const/kb_const,lpint,gra(i,1,1,1)
  enddo

!  gra(1,:,:,1)=-gra(4,:,:,1) !Symmetric gravity. Helps a bit.
!  gra(2,:,:,1)=-gra(3,:,:,1) !Symmetric gravity. Helps a bit.
!  ro_m(2,:,:)=ro_m(3,:,:)
!  ro_m(1,:,:)=ro_m(2,:,:)
!  p_m(2,:,:)=p_m(3,:,:)
!  p_m(1,:,:)=p_m(2,:,:)
!  ro_h(2,:,:)=ro_h(3,:,:)
!  ro_h(1,:,:)=ro_h(2,:,:)
!  p_h(2,:,:)=p_h(3,:,:)
!  p_h(1,:,:)=p_h(2,:,:)
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

!!!!!!!!!!!!!!!!!!!!!!!
!Set normalisation parameters based on totals TO DO
!T0=1.0e4 !This is read from the parameters file
!find T(x)=T0
!x0=interpol(temperature,x,T0)
i=floor(ix/2.0d0)
!x0=x(i-1)+(temperature(i)-temperature(i-1))*(x(i)-x(i-1))/(temperature(i)-temperature(i-1)) 
x0=x(i)
print*,'normalisation height', x0
!print*,'temperature at normalised height',interpol(x,temperature,x0)
L0=2.0d0*kb_const*temperature(i)/abs(gra(i,1,1,1))/mh_const !pressure scale height
!L0=-2.0d0*kb_const*temperature(i)/gra0/mh_const !IF GRAVITY IS ZERO USE THIS!!!
ro0=ro_m(i,1,1)            !interpol(x,ro_h(:,1,1),x0)
P0=p_m(i,1,1)*gm/2.d0           !interpol(x,p_h(:,1,1),x0)*gm
gran=abs(gra(i,1,1,1))*gm/2.d0  !abs(interpol(x,gra(:,1,1,1),x0))*gm
!gran=abs(gra0)*gm  !IF GRAVITY IS ZERO USE THIS!!!
print*,'T=',temperature(i),' at x=', x0,',L0=',L0,', ro0=',ro0, 'P0=',P0, 'G0',gran
x=x/L0
start=start/L0
end=end/L0
ro_h=ro_h/ro0
ro_m=ro_m/ro0
p_h=p_h/P0
p_m=p_m/P0
gra=gra/gran
!dxc=dxc/L0
!dx=dx/L0
!!!!!!!!!!!!!!!!!!!!!!!

print*,'Photospheric temperature',T0*p_h(2,1,1)/ro_h(2,1,1)*5.0d0/3.0d0

!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write initial atmosphere
!OPEN(UNIT=102, FILE="downflowshock_PIP_hr.txt", ACTION="write", STATUS="replace")
OPEN(UNIT=102, FILE="df_PIP_xi0_05.txt", ACTION="write", STATUS="replace")
write(102,*) ix, start(1), end(1), L0, T0 , ro0, P0, gran
do i=1,ix 
	write(102,*) x(i), ro_h(i,1,1), ro_m(i,1,1), p_h(i,1,1), p_m(i,1,1), gra(i,1,1,1)
enddo
close(102)
stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 B0=0.0d0!/20.0d0


!initialise other arrays
  do k=1,kx;do j=1,jx;do i=1,ix
     vx_h(i,j,k)=0.0d0
     vx_m(i,j,k)=0.0d0
     vy_h(i,j,k)=0.0d0
     vy_m(i,j,k)=0.0d0
     vz_h(i,j,k)=0.0d0
     vz_m(i,j,k)=0.0d0
     b_x(i,j,k)=0.3d0*B0
!     b_y(i,j,k)=0.5d0*B0*(1.0d0+tanh(x(i)/dx(1)/10.d0))
     b_y(i,j,k)=B0!*tanh(x(i)/dx(1)/10.d0)
     b_z(i,j,k)=0.0d0
  enddo;enddo;enddo

!  b_x(1,:,:)=-B0
!  b_x(2,:,:)=-B0


  !!set boundary condition----------------------
  flag_bnd(1)=10  !3 for asymmetric B
  flag_bnd(2)=10
  flag_bnd(3)=10
  flag_bnd(4)=10
  flag_bnd(5)=10
  flag_bnd(6)=10
  !-------------------------------------------------
stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elseif(debug_direction.eq.10) then
! Loaded saved stable profile with PIP (made using debug_direction=10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
!     print*,'=== flag_PIP should be 1 ==='
     f_n=0.0d0
     f_p=1.0d0
     f_p_n=0.0d0
     f_p_p=1.0d0
!     stop
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

!  open(unit=102,file='const_MHD_g0.txt',status='old')
!  open(unit=102,file='const_cs8958_MHD_grar.txt',status='old')
!  open(unit=102,file='fal.txt',status='old')
!  open(unit=102,file='const_cs8958_PIP_grar_xin_09.txt',status='old')
!  open(unit=102,file='const_rev.txt',status='old')
!  open(unit=102,file='const_PIP_xin_0.99_xlong_c.txt',status='old')
!open(unit=102,file='downflowshock_PIP_hr.txt',status='old') !xi0=0.9, corona still ionised
open(unit=102,file='df_PIP_xi0_05.txt',status='old') !xi0=0.5
  read(102,*) nxt, h0, hn, L0, Tn , ro0, P0, gran
!h0=0.0d0
!print*,'Overwriting h0 as:',h0
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
!	ro_h(i,:,:)=ro_h_temp
!	ro_m(i,:,:)=ro_m_temp
!	p_h(i,:,:)=p_h_temp
!	p_m(i,:,:)=p_m_temp
!	gra(i,:,:,1)=gra_temp
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
!print*,x(1)
  !---------------------------------------
  !Set nonuniform grid
!  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
!       ix,x,dx,dxc,s_order,0.0d0-start(1),20.d0,dx(1)*1.0d0,1.001d0,dx(1)*20.0d0, &
!       mpi_pos(1),margin(1),0,flag_err)
!print*,dx
!x=x+h0 !fix for having negative x after NUG


!Set domain values
print*,'initial height',x(1),h0,start(1)
print*,'synthesised atmosphere height',hn,'calculated height',maxval(x)
!  gra(:,:,:,:)=0.0d0 
 do i=1,ix
!	print*,x(i),h_read(i)
	ro_h(i,:,:)=interpol(h_read,ro_h_read,x(i))
!	print*,'ro_h',ro_h_read(i)
	ro_m(i,:,:)=interpol(h_read,ro_m_read,x(i))
!	print*,'i',i,ro_m(i,1,1),ro_m_read(i)
	p_h(i,:,:)=interpol(h_read,p_h_read,x(i))
!	print*,'p_h',p_h_read(i)
	p_m(i,:,:)=interpol(h_read,p_m_read,x(i))
!	print*,'p_m',p_m_read(i)
	gra(i,:,:,1)=interpol(h_read,gra_read,x(i))
!	print*,'gra',gra_read(i)
  enddo

if (x(1) .LT. h0) then
	ro_h(1,:,:)=ro_h(2,:,:)
	ro_m(1,:,:)=ro_m(2,:,:)
	p_h(1,:,:)=p_h(2,:,:)
	p_m(1,:,:)=p_m(2,:,:)
	gra(1,:,:,1)=-gra(2,:,:,1)
endif

!if (x(2) .LT. 0) then
!	ro_h(2,:,:)=ro_h(3,:,:)
!	ro_m(2,:,:)=ro_m(3,:,:)
!	p_h(2,:,:)=p_h(3,:,:)
!	p_m(2,:,:)=p_m(3,:,:)
!	gra(2,:,:,1)=-gra(3,:,:,1)
!endif

print*,'finished reading rank',my_rank

 B0=1.0d0!1.0d0!/20.0d0

!Neutral density lump
 do i=1,ix
    if ((x(i) .ge. 1.0) .and. (x(i) .lt. 1.2)) then 
!    	ro_h(i,:,:)=ro_h(i,:,:)+10.0
!        ro_h(i,:,:)=ro_h(i,:,:)+100.0*(1.0-(dcos(x(i)*2.0*pi/0.2)+1.0)/2.0)
    endif
 enddo
!initialise other arrays
  do k=1,kx;do j=1,jx;do i=1,ix
     vx_h(i,j,k)=0.0d0
     vx_m(i,j,k)=0.0d0
     vy_h(i,j,k)=0.0d0
     vy_m(i,j,k)=0.0d0
     vz_h(i,j,k)=0.0d0
     vz_m(i,j,k)=0.0d0
     b_x(i,j,k)=1.0d0*B0
!     b_y(i,j,k)=0.5d0*B0*(1.0d0+tanh(x(i)/dx(1)/10.d0))
     b_y(i,j,k)=0.0d0*B0!*tanh(x(i)/dx(1)/10.d0)
     b_z(i,j,k)=0.0d0
  enddo;enddo;enddo

!  b_x(1,:,:)=-B0
!  b_x(2,:,:)=-B0


  !!set boundary condition----------------------
  flag_bnd(1)=10  !3 for asymmetric B
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

end subroutine downflowshock
