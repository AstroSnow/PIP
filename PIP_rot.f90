module PIP_rot
  use globalvar,only:ix,jx,kx,ac,xi_n,gm_rec,gm_ion,nvar_h,nvar_m,&
       flag_pip_imp,gm,n_fraction,t_ir,col,x,y,z,beta,T0, n0,my_rank,flag_IR_type,flag_col,arb_heat,nout,flag_restart,Colrat,&
        Nexcite,n0,f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref
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
    factor=2.7d0*dsqrt(T0*T_r_p)*T0*T_r_p/5.6e-16/n0
 !   factor=2.7d0*dsqrt(T0*T_r_p)*T0*T_r_p/5.6e-16/n0
!    factor2=1.0d0/((T0/T_r_p)**deg1/factor+(T0/T_r_p)**deg2*exp(-T_r_p/T0))
!    factor2=1.0d0/(T_r_p/T0/factor+dsqrt(T0/T_r_p)*exp(-T_r_p/T0))
    factor2=1.0d0/(T_ionization/factor+dexp(-T_ionization)/dsqrt(T_ionization))

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
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)/dsqrt(beta/2.d0*gm)
    case(2)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)
    case(3)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*dsqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:))) &
              /dsqrt(beta/2.d0*gm)
    case(4)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*dsqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:)))
    case(5)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*dsqrt(1.d0 + &
              9.d0*pi/64.d0*gm/2.d0*sum(vd**2.d0,dim=4)/(Te_h(:,:,:)+Te_m(:,:,:))) &
              /dsqrt(beta/2.d0*gm)*((beta/2.d0*gm) &
              /(Te_h(:,:,:)+Te_m(:,:,:))/2.d0+gm*pi/16.d0*sum(vd**2,dim=4))**0.125d0
              
    case(6)
    allocate(Te_h(ix,jx,kx), Te_m(ix,jx,kx), vd(ix,jx,kx,3))
    call get_Te_MHD(U_m,Te_m)
    call get_Te_HD(U_h,Te_h)
    call get_vel_diff(vd,U_h,U_m)
    ac(:,:,:)=col*dsqrt((Te_h(:,:,:)+Te_m(:,:,:))/2.d0)*dsqrt(1.d0 + &
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
       ion_temperature=dsqrt(Te/T_ionization)*dexp(-T_ionization/Te)*factor2  
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
	factor=dexp(-13.2d0/(T0/11605.0d0)) !exp(-E0/T0) in electron volts
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
	rec_fac=2.6e-19*(n0*1.0e6)/dsqrt(Te_0)*t_ir  !n0 converted to m^-3
!	ele_n=U(:,:,:,1)*rho0/mh_si
!	psi_ion=13.6d0
!	A_ion=2.91e-14
!	k_ion=0.39d0
!	x_ion=0.232d0
	factor=dexp(-13.6d0/Te_0)
	factor2=2.91e-14*(n0*1.0e6)*(13.6d0/Te_0)**0.39d0
	ion_fac=factor2/t_ir
	Gm_rec=U_m(:,:,:,1)/dsqrt(Te_p)*rec_fac
!	Gm_ion=factor**(-Te_p)*U_m(:,:,:,1)*Te_p**(1.0d0-0.39d0)/(Te_p*0.232d0+13.6d0/Te_0)*ion_fac
	Gm_ion=(n0*1.0e6*U_m(:,:,:,1))*2.91e-14*dexp(-13.6d0/Te_0/Te_p)*(13.6d0/Te_0/Te_p)**0.39d0/(0.232d0+13.6d0/Te_0/Te_p)*t_ir
!	print*,Gm_rec(1,1,1)/Gm_ion(1,1,1),U_h(1,1,1,1)/U_m(1,1,1,1),(2.6e-19/dsqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)* &
!		(13.6/Te_0)**0.39*exp(-13.6/Te_0))
!	print*,'Gm_rec',Gm_rec(1,1,1)/U_m(1,1,1,1)*t_ir,(n0*1.0e6)*(2.6e-19/dsqrt(Te_0))
!	print*,'Gm_ion',Gm_ion(1,1,1)/U_m(1,1,1,1)*t_ir,(n0*1.0e6)*(2.91e-14/(0.232d0+13.6d0/Te_0)*(13.6d0/Te_0)**0.39d0*exp(-13.6d0/Te_0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	!print*,f_p_p,'f_p_p'
	
	if(mod(flag_col,2) .eq. 1) then
		tfac=0.5d0*f_p_p/f_p
		!print*,tfac
	elseif(mod(flag_col,2) .eq. 0) then
		tfac=beta/2.0d0*f_p_p*5.0d0/6.0d0/f_p
	else
		print*,'option not included!'
		stop
	endif
	
	!print*,Te_p(1,1,1),'temperature'
	!print*,'mod(2,2)=',mod(2,2)	
	!print*,'mod(3,2)=',mod(3,2)
	!print*,'mod(flag_col,2)=',mod(flag_col,2)
	!stop

	Gm_rec=U_m(:,:,:,1)/dsqrt(Te_p)*t_ir/f_p*dsqrt(tfac)
	Gm_ion=2.91e-14*(n0*1.0e6)*U_m(:,:,:,1)*dexp(-13.6d0/Te_0/Te_p*tfac)*(13.6d0/Te_0/Te_p*tfac)**0.39d0
	Gm_ion=Gm_ion/(0.232d0+13.6d0/Te_0/Te_p*tfac)/rec_fac/f_p *t_ir
!print*,'set_ir:',maxval(gm_rec),maxval(gm_ion)
    if(nout .eq. 0 .and. flag_restart.eq.-1) then
	allocate(arb_heat(ix,jx,kx))
	if(mod(flag_col,2) .eq. 1) then
		arb_heat=Gm_ion*U_h(:,:,:,1)*(13.6d0/gm/T0/8.6173e-5)
	elseif(mod(flag_col,2) .eq. 0) then
		arb_heat=Gm_ion*U_h(:,:,:,1)*(13.6d0*beta/T0/2.d0/8.6173e-5)
	else
		print*,'option not included!'
		stop
	endif
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(4) !6 level formula using Leenaarts+2007 and Johnson 1972

    !Get the species temperature
	call get_Te_HD(U_h,Te_n)
	call get_Te_MHD(U_m,Te_p)
    

    !Normalise the temperature 
	if(mod(flag_col,2) .eq. 1) then
		tfac=0.5d0*f_p_p_ini/f_p_ini
		!print*,tfac
	elseif(mod(flag_col,2) .eq. 0) then
		tfac=beta/2.0d0*f_p_p_ini*5.0d0/6.0d0/f_p_ini
	else
		print*,'option not included!'
		stop
	endif

    if(nout .eq. 0 .and. flag_restart.eq.-1) then
	    allocate(Colrat(ix,jx,kx,6,6)) !Allocate the rate array
        call get_col_ion_coeff(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec)
        Gm_rec_ref=Gm_rec(1,1,1)  !Get the normalisation recombination rates
        !allocate(Nexcite(ix,jx,kx,6)) !Allocate the fractional array (allocated in IC)
    endif
!print*,Gm_rec_ref
!stop
    !Calculate the coefficients
!print*,Te_p(1,1,1),T0,tfac,Te_p(1,1,1)*T0/tfac,f_p_ini,f_p_p_ini
!print*,U_m(1,1,1,1)*n0/n0fac,U_m(1,1,1,1),n0,n0fac
!stop

!Get the dimensional ionisation and recombination rates
    call get_col_ion_coeff(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec) 

!Normalise the rates based on intial recombination rate
    Gm_ion=Gm_ion/Gm_rec_ref*t_ir
    Gm_rec=Gm_rec/Gm_rec_ref*t_ir

    !print*,Gm_rec(1,1,1),Gm_ion(1,1,1)
    end select
  end subroutine set_IR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_col_ion_coeff(Telec,Nelec,Gm_ion,Gm_rec)
    !Calculate the excitation and ionisation coefficients from Johnson 1972
    !Assume a 6 level hydrogen atom (1=ground, 2=1st excitation, ...., 6=ionised) 
    double precision,intent(in)::Telec(ix,jx,kx),Nelec(ix,jx,kx)
    double precision,intent(out)::Gm_ion(ix,jx,kx),Gm_rec(ix,jx,kx)
    !Universal constants
    double precision,parameter::melec=9.10938356e-31 !Electron mass [kg]
    double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
    double precision,parameter::a0bohr=5.29e-11 !Bohr's radius [m] 
    double precision,parameter::hplank=6.62607004e-34 !Planks constant [m2 kg s^-1]
    double precision,parameter::kbhat=1.38064852d0,mehat=9.10938356d0,hhat=6.62607004d0
    double precision::gweight(6) !Statistical weighting of excitation states
    double precision::Eion(6) !Energy to ionise
    double precision::garr(3)
    double precision::Colex(6,6),dneut(6)
    integer::k,j,i,ii,jj
    double precision::rn,bn,xrat,Enn,yhat,rnn,zhat,gauntfac,fnn
    double precision::Ann,Bnn,E0y,E1y,E2y,E0z,E1z,E2z
    double precision::yn,zn,ziyn,zizn,An0,Bn0
    double precision::dntot
    integer::nmaxloc
    !This is all consant so can be moved somewhere else?
!    gweight(1)=2.0*1.0**2 !ground state
!    gweight(2)=2.0*2.0**2 !1st level excitation
!    gweight(3)=2.0*3.0**2 !2nd level excitation
!    gweight(4)=2.0*4.0**2 !3rd level excitation
!    gweight(5)=2.0*5.0**2 !4th level excitation
!    gweight(6)=1.0 !Ionised hydrogen
    gweight(1)=2.d0 !ground state
    gweight(2)=8.d0 !1st level excitation
    gweight(3)=18.d0 !2nd level excitation
    gweight(4)=32.d0 !3rd level excitation
    gweight(5)=50.d0 !4th level excitation
    gweight(6)=1.d0 !Ionised hydrogen

    Eion=[13.6d0,3.4d0,1.51d0,0.85d0,0.54d0,0.0d0] !in eV
    Eion=Eion/13.6d0*2.18e-18 !Convert to joules (to be dimensionally correct)

    !set colex to zero initially
    colex(:,:)=0.d0
    !Loop over the grid
    do k=1,kx;do j=1,jx; do i=1,ix
        !loop over the Excitation states
            do ii=1,5
                if (ii .eq. 1) then 
                    rn=0.45d0 !Equation 31
                else
                    rn=1.94d0*dble(ii)**(-1.57d0) !Equation 32
                endif

                if (ii .eq. 1) then 
                    bn=-0.603d0 !Equation 25
                else 
                    bn=(1.0d0/dble(ii))*(4.0d0-18.63d0/dble(ii)+36.24d0/(dble(ii)**2)-&
                                    28.09d0/(dble(ii)**3)) !Equation 26
                endif

                do jj=ii+1,5
                    xrat=1.0d0-(dble(ii)/dble(jj))**2 !ratio of transition energy to ionisation energy of lower level
                    Enn=Eion(ii)-Eion(jj) !Difference in ionisation energies
                    yhat=Enn/kboltz/Telec(i,j,k) !Equation 37

                    rnn=rn*xrat !Equation 30

                    zhat=rnn+Enn/kboltz/Telec(i,j,k)   !Equation 38

                    !Gaunt factors from table 1 and equation 4
                    !Constants so can be moved for speed
                    if (ii .eq. 1) then 
                        gauntfac=1.1330d0-0.4059d0/xrat+0.07014d0/(xrat**2)
                    else if (ii .eq. 2) then 
                        gauntfac=1.0785d0-0.2319d0/xrat+0.02947d0/(xrat**2)
                    else
                        gauntfac=0.9935d0+0.2328d0/dble(ii)-0.1296d0/(dble(ii)**2)&
                            -(1.0d0/xrat)*(1.0d0/dble(ii))*(0.6282d0-0.5598d0/dble(ii)+0.5299d0/(dble(ii)**2))&
                            +(1.0d0/xrat)**2*(1.0d0/dble(ii)**2)*(0.3887d0-1.181d0/dble(ii)+1.470d0/(dble(ii)**2))
                    endif

                    fnn=32.0d0/3.0d0/dsqrt(3.d0)/pi*dble(ii)/dble(jj)**3/(xrat**3)*gauntfac !Equation 3

                    Ann=2.0d0*dble(ii)**2/xrat*fnn !Equation 11
                    Bnn=4.0d0*dble(ii)**4/(dble(jj)**3)/(xrat**2)*(1.0d0+4.0d0/3.0d0/xrat+bn/(xrat**2)) !Equation 23

                    !Exponential fits (THESE NEED WORK)
                    call ionexpfittest(yhat,1.d0,0.001d0,0.0001d0,E1y)   !Equation 8
                    call ionexpfittest(zhat,1.d0,0.001d0,0.0001d0,E1z)
                    call ionexpfittest(yhat,2.d0,0.001d0,0.0001d0,E2y)
                    call ionexpfittest(zhat,2.d0,0.001d0,0.0001d0,E2z)

                    !Rate coefficient for excitation Equation 36 using electron mass
                    Colex(ii,jj)=Nelec(i,j,k)*dsqrt(8.d0*kboltz*Telec(i,j,k)/pi/melec)*2.0d0*&
                        dble(ii)**2/xrat*pi*a0bohr**2*yhat**2*&
                              (Ann*((1.0d0/yhat+0.5d0)*E1y-(1.0d0/zhat+0.5d0)*E1z)+&
                              (Bnn-Ann*dlog(2.d0*dble(ii)**2/xrat))*(E2y/yhat-E2z/zhat))

                enddo

                !now do the ionisation states

                !Ionisation coefficient
                yn=Eion(ii)/kboltz/Telec(i,j,k) 
                zn=rn+Eion(ii)/kboltz/Telec(i,j,k)

                !Call the hokey exponential integral
                call ionexpfittest(yn,1.d0,0.001d0,0.0001d0,E1y)   !Equation 8
                call ionexpfittest(zn,1.d0,0.001d0,0.0001d0,E1z)
                call ionexpfittest(yn,2.d0,0.001d0,0.0001d0,E2y)
                call ionexpfittest(zn,2.d0,0.001d0,0.0001d0,E2z)
                call ionexpfittest(yn,0.d0,0.001d0,0.0001d0,E0y)
                call ionexpfittest(zn,0.d0,0.001d0,0.0001d0,E0z)

                ziyn=E0y-2.0d0*E1y+E2y !Equation 42
                zizn=E0z-2.0d0*E1z+E2z !Equation 42

                if (ii .eq. 1) then 
                    garr(1)= 1.1330d0
                    garr(2)=-0.4059d0
                    garr(3)= 0.07014d0
                else if (ii .eq. 2) then 
                    garr(1)= 1.0785d0
                    garr(2)=-0.2319d0
                    garr(3)= 0.02947d0
                else
                    garr(1)=0.9935d0+0.2328d0/dble(ii)-0.1296d0/(dble(ii)**2)
                    garr(2)=(-1.0d0/dble(ii))*(0.6282d0-0.5598d0/dble(ii)+0.5299d0/(dble(ii)**2))
                    garr(3)=(1.0d0/dble(ii))**2*(0.3887d0-1.181d0/dble(ii)+1.470d0/(dble(ii)**2))
                endif
                
                !Equation 20
                An0=32.0d0/3.0d0/dsqrt(3.d0)/pi*dble(ii)+garr(1)/3.0d0+garr(2)/4.0d0+garr(3)/5.0d0
                Bn0=2.0d0/3.0d0*dble(ii)**2*(5.0d0+bn) !Equation 24

                !Ionisation coefficients Equation 35 using electron mass
                Colex(ii,6)=Nelec(i,j,k)*dsqrt(8.d0*kboltz*Telec(i,j,k)/pi/melec)&
                    *2.0d0*dble(ii)**2*pi*a0bohr**2*yn**2*&
                    (An0*(E1y/yn-E1z/zn)+&
                    (Bn0-An0*dlog(2.d0*dble(ii)**2))*(ziyn-zizn))
            enddo

            !Rates
            do ii=1,5
                !Excitation rates
                do jj=1,5
                    if (ii .eq. jj) then 
                        colrat(i,j,k,ii,jj)=0.d0
                    else if (jj .gt. ii) then 
                        colrat(i,j,k,ii,jj)=gweight(ii)/gweight(jj)*Colex(ii,jj)*dsqrt(Telec(i,j,k)) 
                    else
                        colrat(i,j,k,ii,jj)=Colex(jj,ii)*dsqrt(Telec(i,j,k))*&
                            exp(-(Eion(jj)-Eion(ii))/kboltz/Telec(i,j,k))
                    endif
                enddo
                !Ionisation rates
                colrat(i,j,k,ii,6)=Colex(ii,6)*dsqrt(Telec(i,j,k))*exp(-(Eion(6)-Eion(ii))/kboltz/Telec(i,j,k))
            enddo

            do ii=1,5
                colrat(i,j,k,6,ii)=nelec(i,j,k)*colrat(i,j,k,ii,6)*gweight(ii)/gweight(6)&
                    /2.0d0*(2.0d0*pi*mehat*kbhat*Telec(i,j,k)/hhat/hhat*1.0e14)**(-3.0d0/2.0d0)*&
                    exp(Eion(ii)/kboltz/Telec(i,j,k))
            enddo
!print*,nelec(i,j,k),Telec(i,j,k)
!print*,colex
!print*,colrat(i,j,k,:,:)
!stop
            ! Calculate the change in each species
            dntot=0.d0
            do ii=1,6
                do jj=1,6 
                    dneut(ii)=Nexcite(i,j,k,jj)*colrat(i,j,k,jj,ii) - Nexcite(i,j,k,ii)*colrat(i,j,k,ii,jj)
                enddo

                    if (Nexcite(i,j,k,ii) .ne. maxval(Nexcite(i,j,k,:))) then 
                        dntot=dntot+dneut(ii)
                    endif
            enddo
            nmaxloc=maxloc(Nexcite(i,j,k,:),DIM=1)
            dneut(nmaxloc)=-dntot

            Gm_ion(i,j,k)=0.d0
            Gm_rec(i,j,k)=0.d0
            do ii=1,5
                Gm_ion(i,j,k)=Gm_ion(i,j,k)+Nexcite(i,j,k,ii)*colrat(i,j,k,ii,6)
                Gm_rec(i,j,k)=Gm_rec(i,j,k)+Nexcite(i,j,k,6)*colrat(i,j,k,6,ii)
            enddo
!print*,gm_ion(1,1,1),gm_rec(1,1,1)
!stop
            !divide rates by total neutral/plasma density for consistency
            Gm_ion(i,j,k)=Gm_ion(i,j,k)/(Nexcite(i,j,k,1)+Nexcite(i,j,k,2)+Nexcite(i,j,k,3)&
                          +Nexcite(i,j,k,4)+Nexcite(i,j,k,5))        
            Gm_rec(i,j,k)=Gm_rec(i,j,k)/Nexcite(i,j,k,6)
    enddo;enddo;enddo

  end subroutine get_col_ion_coeff  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ionexpfittest(zhat,istar,stepsize,tol,sol)
    !numerical integration of the exponential
    !This needs work
    double precision,intent(in)::zhat,istar,stepsize,tol
    double precision,intent(out)::sol
    double precision::a,b,fa,fb,sol0,dif,solt

    a=1.d0
    b=a+stepsize
    fa=exp(-zhat*a)*a**(-istar)
    fb=exp(-zhat*b)*b**(-istar)
    sol0=0.5d0*(b-a)*(fa+fb)
    sol=sol0

    dif=1.d0

    do while (dif .gt. tol)
        a=b
        b=a+stepsize
        fa=exp(-zhat*a)*a**(-istar)
        fb=exp(-zhat*b)*b**(-istar)
        solt=0.5d0*(b-a)*(fa+fb)
        sol=sol+solt
        dif=solt/sol0
    enddo

  end subroutine ionexpfittest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine hydrogen_excitation_update(dt,rom,roh)
! update the hydrogen excitation states
  double precision,intent(in)::dt,rom(ix,jx,kx),roh(ix,jx,kx)
  double precision::dneut(6)
  double precision::dntot
  integer::i,j,k,ii,jj,nmaxloc
!The new number of electrons/protons is:
    Nexcite(:,:,:,6)=rom
    do k=1,kx;do j=1,jx; do i=1,ix
            ! Calculate the change in each neutral species
            dntot=0.d0
            do ii=1,5
                do jj=1,6 
                    dneut(ii)=Nexcite(i,j,k,jj)*colrat(i,j,k,jj,ii) - Nexcite(i,j,k,ii)*colrat(i,j,k,ii,jj)
                enddo
            enddo
!print*,dneut
!print*,colrat(i,j,k,:,:)
!print*,Nexcite(i,j,k,:)
            Nexcite(i,j,k,1:5)=Nexcite(i,j,k,1:5)+dt*dneut(1:5)/n0
            where (Nexcite>=0.d0)
                Nexcite = Nexcite
            elsewhere
                Nexcite = 0.d0
            end where
!print*,Nexcite(i,j,k,:)
            dntot=sum(Nexcite(i,j,k,1:5))
            Nexcite(i,j,k,1:5)=roh(i,j,k)*Nexcite(i,j,k,1:5)/dntot
    enddo;enddo;enddo
!print*,Nexcite(1,1,1,6),rom(1,1,1),sum(Nexcite(1,1,1,1:5)),roh(1,1,1)
  end subroutine hydrogen_excitation_update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
    N_i=-f_T+dsqrt(f_T*(f_T+N_pr))
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
    double precision dneut(ix,jx,kx,6)
    integer i,j   
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
		if(mod(flag_col,2) .eq. 1) then
			ion_pot=Gm_ion*nde*(13.6d0/gm/T0/8.6173e-5)
		elseif(mod(flag_col,2) .eq. 0) then
			ion_pot=Gm_ion*nde*(13.6d0*beta/T0/2.d0/8.6173e-5)
		else
			print*,'option not included!'
			stop
		endif
!print*,maxval(Gm_ion),maxval(ion_pot),maxval(abs(arb_heat))
!print*,maxval(abs(ion_pot)),maxval(abs(arb_heat)),maxval(abs(ion_pot-arb_heat))
	S_m(:,:,:,5)=S_m(:,:,:,5)-ion_pot+arb_heat
!	print*,'New type'
	else if(flag_IR_type .eq. 1) then
	ds(:,:,:,5)=0.5d0*(Gm_rec*de*(vx*vx+vy*vy+vz*vz)- &
	    Gm_ion*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
	    Gm_ion*npr/(gm-1.d0) + 0.5d0*Gm_rec*pr/(gm-1.d0) !Khomenko
	S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)

	else if(flag_IR_type .eq. 2) then
	ds(:,:,:,5)=0.5d0*(Gm_rec*de*(vx*vx+vy*vy+vz*vz)- &
	    Gm_ion*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
		    Gm_ion*npr/(gm-1.d0) + Gm_rec*pr/(gm-1.d0) !Singh
	S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)

    else if(flag_IR_type .eq. 3) then
    !6 level hydrogen atom
    
    !Calculate the delta densities conserving total by omiting largest value  
        do i=1,6
            do j=1,6 
                dneut(:,:,:,i)=Nexcite(:,:,:,j)*colrat(:,:,:,j,i)-Nexcite(:,:,:,i)*colrat(:,:,:,i,j)
            enddo
        enddo

    !The source has to be calculated elsewhere though I think.
	endif

    endif    
    return
  end subroutine source_PIP

end module PIP_rot
