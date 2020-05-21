module PIP_rot
  use globalvar,only:ix,jx,kx,ac,xi_n,gm_rec,gm_ion,nvar_h,nvar_m,&
       flag_pip_imp,gm,n_fraction,t_ir,col,x,y,z,beta,T0, n0,my_rank,flag_IR_type,flag_col,arb_heat,nout
  use scheme_rot,only:get_Te_HD,get_Te_MHD,cq2pv_HD,cq2pv_MHD,get_vel_diff
  use parameters,only:T_r_p,deg1,deg2,pi
  implicit none
  integer,save::col_type,IR_type,xin_type,is_IR,IR_T_dependence
  double precision factor,factor2,mu_p,mu_n,T_ionization,factor3
  double precision :: rec_fac,ion_fac
  double precision :: ioneq,f_n,f_p, f_p_p,tfac 
contains
  subroutine initialize_collisional(flag_col)
    integer,intent(inout)::flag_col
    if (flag_col.eq.0) return
    allocate(ac(ix,jx,kx),xi_n(ix,jx,kx))
    if (flag_col.ge.2 .and. flag_col.le.7) then
      col_type=flag_col-1
    else
      col_type=0
      if(my_rank.eq.1) print*, 'WARNING: constant alpha used'
    endif
!mod((flag_col/10),10)
    flag_col=mod(flag_col,10)
    T_ionization=T_r_p
    factor=2.7d0*sqrt(T0*T_r_p)*T0*T_r_p/5.6e-16/n0
 !   factor=2.7d0*sqrt(T0*T_r_p)*T0*T_r_p/5.6e-16/n0
!    factor2=1.0d0/((T0/T_r_p)**deg1/factor+(T0/T_r_p)**deg2*exp(-T_r_p/T0))
!    factor2=1.0d0/(T_r_p/T0/factor+sqrt(T0/T_r_p)*exp(-T_r_p/T0))
    factor2=1.0d0/(T_ionization/factor+exp(-T_ionization)/sqrt(T_ionization))

    mu_p=0.5d0
    mu_n=1.0d0

  end subroutine initialize_collisional




  subroutine set_collisional(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision,parameter::r_x=0.2d0,r_y=1.0d0,r_z=5.0d0
    double precision,allocatable:: Te_m(:,:,:), Te_h(:,:,:), vd(:,:,:,:)
    integer i,j,k
    select case(col_type)
    case(0)
       ac(:,:,:)=col
    case(1)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)/sqrt(beta/2.d0*gm)
    case(2)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)
    case(3)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*sqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:))) &
              /sqrt(beta/2.d0*gm)
    case(4)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*sqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:)))
    case(5)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*sqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:))) &
              /sqrt(beta/2.d0*gm)*((beta/2.d0*gm) &
              /(Te_h(:,:,:)+Te_m(:,:,:))/2.d0+gm*pi/16.d0*sum(vd**2,dim=4))**0.125d0
              
    case(6)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*sqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*sqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:)))&
              /((Te_h(:,:,:)+Te_m(:,:,:))/2.d0+gm*pi/16.d0*sum(vd**2,dim=4))**0.125d0

    end select
  end subroutine set_collisional

  subroutine initialize_IR(flag_IR)
    integer,intent(inout)::flag_IR
    if (flag_IR.eq.0) return    
    allocate(Gm_rec(ix,jx,kx),Gm_ion(ix,jx,kx))
    IR_type=flag_IR
  end subroutine initialize_IR

  function rec_temperature(Te)
    double precision Te(ix,jx,kx)
    double precision rec_temperature(ix,jx,kx)
    if(IR_T_dependence.eq.0) then
       rec_temperature=T_ionization/te/factor*factor2       
    else if(IR_T_dependence.eq.1) then
       rec_temperature=n_fraction
    endif
  end function rec_temperature

  function ion_temperature(Te)
    double precision Te(ix,jx,kx)
    double precision ion_temperature(ix,jx,kx)
    if(IR_T_dependence.eq.0) then
       ion_temperature=sqrt(Te/T_ionization)*exp(-T_ionization/Te)*factor2  
    else if(IR_T_dependence.eq.1) then
       ion_temperature=1.0-n_fraction
    endif
  end function ion_temperature

  subroutine set_IR(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision Te_n(ix,jx,kx),Te_p(ix,jx,kx),Te_e(ix,jx,kx)
    double precision xi_n_tmp(ix,jx,kx)
    double precision Te_0
    select case(IR_type)
    case(1)
	!Formulation from Jeffery paper
	!WORK IN PROGRESS - Need to define normalisation quantities
	! 
	!Get species temperatures
	call get_Te_HD(U_h,Te_n)
	call get_Te_MHD(U_m,Te_p)
	factor=exp(-13.2d0/(T0/11605.0d0)) !exp(-E0/T0) in electron volts
	factor2=2.7*(13.2d0*11605.0d0)**-2.0d0*T0**(1.0d0/2.0d0)*n0
	factor3=5.6e-16*(13.2d0*11605.0d0)**-2.0d0*T0**(-1.0d0)*n0**2.0d0
	rec_fac=factor2/factor3 !is this right?
	ion_fac=factor3/t_ir
	Gm_rec=Te_p**(-0.5d0)*U_m(:,:,:,1)**2.d0*rec_fac*ion_fac
	Gm_ion=Te_p**0.5d0*U_m(:,:,:,1)*factor**(1.0d0/Te_p)*ion_fac
!	print*,'factor',factor
!	print*,'factor2',factor2
!	print*,'factor3',factor3
    case(2)
	!Formulation from Popescu+2019 paper
	!Empirical estimates for the rates
	!WORK IN PROGRESS
	call get_Te_HD(U_h,Te_n)
	call get_Te_MHD(U_m,Te_p)
	!Calculate electron temperature in eV
	Te_0=T0/1.1604e4
	rec_fac=2.6e-19*(n0*1.0e6)/sqrt(Te_0)*t_ir  !n0 converted to m^-3
!	ele_n=U(:,:,:,1)*rho0/mh_si
!	psi_ion=13.6d0
!	A_ion=2.91e-14
!	k_ion=0.39d0
!	x_ion=0.232d0
	factor=exp(-13.6d0/Te_0)
	factor2=2.91e-14*(n0*1.0e6)*(13.6d0/Te_0)**0.39d0
	ion_fac=factor2/t_ir
	Gm_rec=U_m(:,:,:,1)/sqrt(Te_p)*rec_fac
!	Gm_ion=factor**(-Te_p)*U_m(:,:,:,1)*Te_p**(1.0d0-0.39d0)/(Te_p*0.232d0+13.6d0/Te_0)*ion_fac
	Gm_ion=(n0*1.0e6*U_m(:,:,:,1))*2.91e-14*exp(-13.6d0/Te_0/Te_p)*(13.6d0/Te_0/Te_p)**0.39d0/(0.232d0+13.6d0/Te_0/Te_p)*t_ir
!	print*,Gm_rec(1,1,1)/Gm_ion(1,1,1),U_h(1,1,1,1)/U_m(1,1,1,1),(2.6e-19/sqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)* &
!		(13.6/Te_0)**0.39*exp(-13.6/Te_0))
!	print*,'Gm_rec',Gm_rec(1,1,1)/U_m(1,1,1,1)*t_ir,(n0*1.0e6)*(2.6e-19/sqrt(Te_0))
!	print*,'Gm_ion',Gm_ion(1,1,1)/U_m(1,1,1,1)*t_ir,(n0*1.0e6)*(2.91e-14/(0.232d0+13.6d0/Te_0)*(13.6d0/Te_0)**0.39d0*exp(-13.6d0/Te_0))
    case(3)
	!Formulation from Popescu+2019 paper
	!Empirical estimates for the rates
	!Quarentine attempt
	call get_Te_HD(U_h,Te_n)
	call get_Te_MHD(U_m,Te_p)
	!Calculate electron temperature in eV
	Te_0=T0/1.1604e4
	rec_fac=2.6e-19*(n0*1.0e6)/dsqrt(Te_0)  !n0 converted to m^-3
!	factor=dexp(-13.6d0/Te_0)
!	factor2=(13.6d0/Te_0)**0.39d0
!	factor3=1.d0/(0.232d0+13.6d0/Te_0)
!	ion_fac=factor*factor2*factor3

	!initial equilibrium fractions
	ioneq=(2.6e-19/dsqrt(Te_0))/(2.91e-14/(0.232d0+13.6d0/Te_0)*(13.6d0/Te_0)**0.39d0*dexp(-13.6d0/Te_0))
	f_n=ioneq/(ioneq+1.0d0)
	f_p=1.0d0-f_n
	f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
	tfac=beta/2.0d0*f_p_p*5.0d0/6.0d0/f_p

	Gm_rec=U_m(:,:,:,1)/dsqrt(Te_p)*t_ir/f_p*dsqrt(tfac)
	Gm_ion=2.91e-14*(n0*1.0e6)*U_m(:,:,:,1)*dexp(-13.6d0/Te_0/Te_p*tfac)*(13.6d0/Te_0/Te_p*tfac)**0.39d0
	Gm_ion=Gm_ion/(0.232d0+13.6d0/Te_0/Te_p*tfac)/rec_fac/f_p *t_ir
!print*,'set_ir:',maxval(gm_rec),maxval(gm_ion)
    if(nout .eq. 0) then
	allocate(arb_heat(ix,jx,kx))
	if(flag_col .eq. 2) then
		arb_heat=Gm_ion*U_h(:,:,:,1)*(13.6d0/gm/T0/8.6173e-5)
	elseif(flag_col .eq. 3) then
		arb_heat=Gm_ion*U_h(:,:,:,1)*(13.6d0*beta/T0/2.d0/8.6173e-5)
	else
		print*,'option not included!'
		stop
	endif
    endif
    end select
  end subroutine set_IR
  

  subroutine get_initial_xin(Pr_tot,Te_tot,N_tot,xi_n0)
    !return initial xi_n and total number density from
    !  total pressure and temperature
    ! assumption 
    !  -ionization equilibrium
    !    f(T)[f_T] is temperature dependence of ionization and recombination
    !    ionization rate    : nu_ion(T)*n_p*n_n
    !    recombination rate : nu_rec(T)*n_p**3
    !    f(T)=nu_ion(T)/nu_rec(T)
    double precision,intent(in)::Pr_tot(ix,jx,kx),Te_tot(ix,jx,kx)
    double precision,intent(out)::xi_n0(ix,jx,kx),N_tot(ix,jx,kx)
    double precision N_pr(ix,jx,kx),N_i(ix,jx,kx),N_n(ix,jx,kx),f_T(ix,jx,kx)
    f_T=ion_temperature(Te_tot)/rec_temperature(Te_tot)
    N_pr=Pr_tot/Te_tot*gm    
    N_i=-f_T+sqrt(f_T*(f_T+N_pr))
    N_n=N_i*N_i/f_T    
    N_tot=N_i+N_n
    xi_n0=N_n/N_tot        
  end subroutine get_initial_xin

  subroutine get_NT_from_PX(Pr_tot,xi_n0,N_tot,Te_tot)
    !return initial temperature and total number density from
    !  total pressure and neutral fraction
    ! assumption 
    !  -ionization equilibrium
    !    f(T)[f_T] is temperature dependence of ionization and recombination
    !    ionization rate    : nu_ion(T)*n_p*n_n
    !    recombination rate : nu_rec(T)*n_p**3
    !    f(T)=nu_ion(T)/nu_rec(T)
    !    *** f_T must independent of temperature
    double precision,intent(in)::Pr_tot(ix,jx,kx),xi_n0(ix,jx,kx)
    double precision,intent(out)::N_tot(ix,jx,kx)
    double precision,intent(inout)::Te_tot(ix,jx,kx)
    double precision N_pr(ix,jx,kx),N_i(ix,jx,kx),N_n(ix,jx,kx),f_T(ix,jx,kx)
    f_T=ion_temperature(Te_tot)/rec_temperature(Te_tot)
    N_i=xi_n0/(1.0d0-xi_n0)*f_T
    N_n=xi_n0/(1.0d0-xi_n0)*N_i
    N_tot=N_i+N_n
    Te_tot=gm*Pr_tot/(2.0d0*N_i+N_n)    
  end subroutine get_NT_from_PX


  subroutine initialize_xin(flag_amb,flag_col)
    integer,intent(inout)::flag_amb,flag_col
    if (flag_amb.eq.0) return
    if(flag_col.eq.0) flag_col=1
    xin_type=mod((flag_amb/10),10)
    flag_amb=mod(flag_amb,10)
  end subroutine initialize_xin
  
  subroutine set_xin(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    integer i,j,k
    double precision,parameter::r_x=10.0d0,r_y=10.0d0,r_z=5.0d0
    double precision tmp,alpha_1,alpha_0
    if(xin_type.eq.0) then
       xi_n=n_fraction
    elseif(xin_type.eq.1)then
       alpha_0=n_fraction/(1-n_fraction)
       alpha_1=1.0d0
       do k=1,kx;do j=1,jx;do i=1,ix
          !       tmp=alpha_0+(alpha_1-alpha_0)*min(1.0d0,y(j)/r_y)
          tmp=alpha_0+(alpha_1-alpha_0)*((tanh((x(i)-r_x)/2.0d0)+1.0d0)*0.5d0 &
               +(1.0d0-tanh((x(i)+r_x)/2.0d0))*0.5d0)
          xi_n(i,j,k)=tmp/(1.0d0+tmp)
       enddo;enddo;enddo       
    endif
  end subroutine set_xin

  !  subroutine source_PIP(U_h0,U_m0,U_h,U_m,dt_coll_i,dt_sub,S_h,S_m)
  subroutine source_PIP(U_h0,U_m0,U_h,U_m,dt_coll_i,S_h,S_m)  
    double precision,intent(inout):: S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
    double precision,intent(inout):: U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision,intent(inout):: U_h0(ix,jx,kx,nvar_m),U_m0(ix,jx,kx,nvar_m)
    double precision:: dS(ix,jx,kx,nvar_h)
    !    double precision, intent(inout):: dt_coll_i,dt_sub
    double precision, intent(inout):: dt_coll_i
    double precision temp(ix,jx,kx),te(ix,jx,kx),nte(ix,jx,kx)
    double precision pr(ix,jx,kx),npr(ix,jx,kx)
    double precision de(ix,jx,kx),nde(ix,jx,kx)
    double precision vx(ix,jx,kx),nvx(ix,jx,kx)
    double precision vy(ix,jx,kx),nvy(ix,jx,kx)
    double precision vz(ix,jx,kx),nvz(ix,jx,kx)
    double precision bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision kapper(ix,jx,kx),lambda(ix,jx,kx)
    double precision A(ix,jx,kx),B(ix,jx,kx),D(ix,jx,kx),ion_pot(ix,jx,kx)    
    dS=0.0
    call cq2pv_hd(nde,nvx,nvy,nvz,npr,U_h)
    call cq2pv_mhd(de,vx,vy,vz,pr,bx,by,bz,U_m)
    te=pr/de*gm*mu_p
    nte=npr/nde*gm*mu_n
    
    if(flag_pip_imp.eq.1) then
       temp=dt_coll_i*(u_h(:,:,:,1)+u_m(:,:,:,1))
       dS(:,:,:,2)=(u_m(:,:,:,1)*u_h(:,:,:,2)-u_h(:,:,:,1)*u_m(:,:,:,2)) &
            *(1.0d0-exp(-ac*temp))/temp
       dS(:,:,:,3)=(u_m(:,:,:,1)*u_h(:,:,:,3)-u_h(:,:,:,1)*u_m(:,:,:,3)) &
            *(1.0d0-exp(-ac*temp))/temp          
       dS(:,:,:,4)=(u_m(:,:,:,1)*u_h(:,:,:,4)-u_h(:,:,:,1)*u_m(:,:,:,4)) &
            *(1.0d0-exp(-ac*temp))/temp
       
       lambda=ac*(nde+de)
       kapper=3.0d0*(gm-1.0d0)*ac*(mu_p*nde+mu_n*de)
       
       B=-1.0d0/3.0d0*gm*kapper/((nde+de)*(kapper-lambda))* &
            (de*((nvx-vx)*vx+(nvy-vy)*vy+(nvz-vz)*nvz)+  &
            nde*((nvx-vx)*nvx+(nvy-vy)*nvy+(nvz-vz)*nvz))
       D=-gm/6.0d0*kapper*(de-nde)/((de+nde)*(kapper-2*lambda))* &
            ((nvx-vx)**2+(nvy-vy)**2+(nvz-vz)**2)
       dS(:,:,:,5)=-(3.0d0/gm)*ac*de*nde/kapper*( &
            +(exp(-kapper*dt_coll_i)-1.0d0)*(nte-te)  &
            -(B+D)*exp(-kapper*dt_coll_i) &
            + B*exp(-lambda*dt_coll_i)+D*exp(-2.0d0*lambda*dt_coll_i))/dt_coll_i

       U_h0(:,:,:,1:5)=U_h0(:,:,:,1:5)-dt_coll_i*ds(:,:,:,1:5)
       U_m0(:,:,:,1:5)=U_m0(:,:,:,1:5)+dt_coll_i*ds(:,:,:,1:5)
!       S_h(:,:,:,1:5)=S_h(:,:,:,1:5)-dS(:,:,:,1:5) 
!       S_m(:,:,:,1:5)=S_m(:,:,:,1:5)+dS(:,:,:,1:5) 
    else
       dS(:,:,:,2)=ac*(u_m(:,:,:,1)*u_h(:,:,:,2)-u_h(:,:,:,1)*u_m(:,:,:,2))
       dS(:,:,:,3)=ac*(u_m(:,:,:,1)*u_h(:,:,:,3)-u_h(:,:,:,1)*u_m(:,:,:,3))
       dS(:,:,:,4)=ac*(u_m(:,:,:,1)*u_h(:,:,:,4)-u_h(:,:,:,1)*u_m(:,:,:,4)) 
       dS(:,:,:,5)=ac*nde*de*(0.5d0*((nvx**2-vx**2)+(nvy**2-vy**2)+ &
            (nvz**2-vz**2)) + 3.0d0/gm/2.0d0*(nte-te))       

       S_h(:,:,:,1:5)=S_h(:,:,:,1:5)-dS(:,:,:,1:5) 
       S_m(:,:,:,1:5)=S_m(:,:,:,1:5)+dS(:,:,:,1:5)
    endif

    if(IR_type.ge.1) then
       ds(:,:,:,1)=Gm_rec*de-Gm_ion*nde
       ds(:,:,:,2)=Gm_rec*de*vx-Gm_ion*nde*nvx
       ds(:,:,:,3)=Gm_rec*de*vy-Gm_ion*nde*nvy
       ds(:,:,:,4)=Gm_rec*de*vz-Gm_ion*nde*nvz
	if(flag_IR_type .eq. 0) then !New formulation
	ds(:,:,:,5)=0.5d0*(Gm_rec*de*(vx*vx+vy*vy+vz*vz)- &
	    Gm_ion*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
	    Gm_ion*npr/(gm-1.d0) + 0.5d0*Gm_rec*pr/(gm-1.d0)
	S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)
	ion_pot=0.0d0
		if(flag_col .eq. 2) then
			ion_pot=Gm_ion*nde*(13.6d0/gm/T0/8.6173e-5)
		elseif(flag_col .eq. 3) then
			ion_pot=Gm_ion*nde*(13.6d0*beta/T0/2.d0/8.6173e-5)
		else
			print*,'option not included!'
			stop
		endif
!print*,maxval(Gm_ion),maxval(ion_pot)
!stop
	!print*,maxval(abs(ion_pot)),maxval(abs(arb_heat)),maxval(abs(ion_pot-arb_heat))
	S_m(:,:,:,5)=S_m(:,:,:,5)-ion_pot+arb_heat
!	print*,'New type'
	else if(flag_IR_type .eq. 1) then
	ds(:,:,:,5)=0.5d0*(Gm_rec*de*(vx*vx+vy*vy+vz*vz)- &
	    Gm_ion*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
	    Gm_ion*npr/(gm-1.d0) + 0.5d0*Gm_rec*pr/(gm-1.d0) !Khomenko
	S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)
!	print*,'Khomenko type'
	else if(flag_IR_type .eq. 2) then
	ds(:,:,:,5)=0.5d0*(Gm_rec*de*(vx*vx+vy*vy+vz*vz)- &
	    Gm_ion*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
		    Gm_ion*npr/(gm-1.d0) + Gm_rec*pr/(gm-1.d0) !Singh
	S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)
!	print*,'Singh type'
	endif

    endif    
    return
  end subroutine source_PIP

end module PIP_rot
