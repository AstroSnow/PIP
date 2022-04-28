module PIP_rot
  use globalvar,only:ix,jx,kx,ac,xi_n,gm_rec,gm_ion,nvar_h,nvar_m,&
       flag_pip_imp,gm,n_fraction,t_ir,col,x,y,z,beta,T0, n0,my_rank,flag_IR_type,flag_col,arb_heat,nout,flag_restart,Colrat,&
        Nexcite,n0,f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref,expinttab,&
        rad_temp,flag_rad,gm_ion_rad,gm_rec_rad,radrat
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
    if (flag_rad .eq. 1) allocate(Gm_rec_rad(ix,jx,kx),Gm_ion_rad(ix,jx,kx))
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

!        allocate(expinttab(4,10000)) !table for the exponential integral table
!        call expintread
 
!	    allocate(Colrat(ix,jx,kx,6,6)) !Allocate the rate array
        call get_col_ion_coeff(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec)
!        call get_col_ion_coeff_aprox(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec)
        Gm_rec_ref=Gm_rec(1,1,1)  !Get the normalisation recombination rates
        !allocate(Nexcite(ix,jx,kx,6)) !Allocate the fractional array (allocated in IC)

        !get the arbitraty heating
        allocate(arb_heat(ix,jx,kx))
        call IRgetionpot(U_m(:,:,:,1),arb_heat) 
    endif
!print*,Gm_rec_ref
!stop
    !Calculate the coefficients
!print*,Te_p(1,1,1),T0,tfac,Te_p(1,1,1)*T0/tfac,f_p_ini,f_p_p_ini
!print*,U_m(1,1,1,1)*n0/n0fac,U_m(1,1,1,1),n0,n0fac
!stop

!Get the dimensional ionisation and recombination rates
    call get_col_ion_coeff(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec) 
!Use the appoximate rates instead
!    call get_col_ion_coeff_aprox(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion,Gm_rec) 

!Normalise the rates based on intial recombination rate
    Gm_ion=Gm_ion/Gm_rec_ref*t_ir
    Gm_rec=Gm_rec/Gm_rec_ref*t_ir

!Get the radiative rates
    if (flag_rad .eq. 1) then
        call get_radrat_fixed(rad_temp,Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion_rad,Gm_rec_rad)

    !Normalise the rates based on intial collisional recombination rate
        Gm_ion_rad=Gm_ion_rad/Gm_rec_ref*t_ir
        Gm_rec_rad=Gm_rec_rad/Gm_rec_ref*t_ir
    endif

    !print*,Gm_rec(1,1,1),Gm_ion(1,1,1)
    end select
  end subroutine set_IR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine expintread
    !Read the exponential integral data calculated using IDL code
    double precision::ytemp,i0temp,i1temp,i2temp
    integer::i

!    open (unit=17, file='exptab2_i0.txt', status='old', action='read')
!    open (unit=18, file='exptab2_i1.txt', status='old', action='read')
!    open (unit=19, file='exptab2_i2.txt', status='old', action='read')
    open (unit=17, file='exptab2_i0_e5.txt', status='old', action='read')
    open (unit=18, file='exptab2_i1_e5.txt', status='old', action='read')
    open (unit=19, file='exptab2_i2_e5.txt', status='old', action='read')

    do i=1,10000
        read(17,*) ytemp, i0temp
        read(18,*) ytemp, i1temp
        read(19,*) ytemp, i2temp
        expinttab(1,i)=ytemp
        expinttab(2,i)=i0temp
        expinttab(3,i)=i1temp
        expinttab(4,i)=i2temp
        !print*,expinttab(:,i)
    enddo

  end subroutine expintread
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_col_ion_coeff_aprox(Telec,Nelec,Gm_ion,Gm_rec)
    !Calculate the approximate excitation and ionisation coefficients from Van Regemorter 1962
    !Assume a 6 level hydrogen atom (1=ground, 2=1st excitation, ...., 6=ionised) 
    double precision,intent(in)::Telec(ix,jx,kx),Nelec(ix,jx,kx)
    double precision,intent(out)::Gm_ion(ix,jx,kx),Gm_rec(ix,jx,kx)
    double precision::Colex(ix,jx,kx,6,6),fosc(6,6) !Array of excitation coefficients
    double precision::nntot(ix,jx,kx)
    double precision::gweight(6) !Statistical weighting of excitation states
    double precision::Eion(6),gfac(6) !Energy to ionise
    double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
    integer::i,j

    Eion=[13.6d0,3.4d0,1.51d0,0.85d0,0.54d0,0.0d0] !in eV
    Eion=Eion/13.6d0*2.18e-18 !Convert to joules (to be dimensionally correct)

    !Gaunt factors
    gfac(1)=2.d0 !ground state
    gfac(2)=8.d0 !1st level excitation
    gfac(3)=18.d0 !2nd level excitation
    gfac(4)=32.d0 !3rd level excitation
    gfac(5)=50.d0 !4th level excitation
    gfac(6)=1.d0 !Ionised hydrogen

    !Oscilator strengths
    !From table 1 in Goldwire 1968
    !http://articles.adsabs.harvard.edu//full/1968ApJS...17..445G/0000450.000.html
    fosc(:,:)=0.0
    fosc(1,2)=4.1619672e-1
    fosc(1,3)=7.9101563e-2
    fosc(1,4)=2.8991029e-2
    fosc(1,5)=1.3938344e-2
    fosc(2,3)=6.4074704e-1
    fosc(2,4)=1.1932114e-1
    fosc(2,5)=4.4670295e-2
    fosc(3,4)=8.4209639e-1
    fosc(3,5)=1.5058408e-1
    fosc(4,5)=1.0377363e0

    !Erase old rates
    Colrat=0.d0

    !Excitation coefficient
    do i=1,5
        do j=i+1,5
            Colrat(:,:,:,i,j)=2.16d0*((Eion(i)-Eion(j))/kboltz/Telec(:,:,:))**(-1.68)*&
                    Telec(:,:,:)**(-3.0/2.0)*dexp(-(Eion(i)-Eion(j))/kboltz/Telec(:,:,:))*fosc(i,j)*Nelec*1.0e-6
            Colrat(:,:,:,j,i)=2.16d0*((Eion(i)-Eion(j))/kboltz/Telec(:,:,:))**(-1.68)*&
                    Telec(:,:,:)**(-3.0/2.0)*fosc(i,j)*gfac(i)/gfac(j)*Nelec*1.0e-6
        enddo
    enddo

    do i=1,5
        !Ionisation coeficients
        colrat(:,:,:,i,6)=2.7d0*(Eion(i)/kboltz/Telec(:,:,:))**(-1)*Telec(:,:,:)**(-3.0/2.0)*&
                dexp(-Eion(i)/kboltz/Telec(:,:,:))*Nelec*1.0e-6
        !Recombination coeficient
        colrat(:,:,:,6,i)=5.6e-16*(Eion(i)/kboltz/Telec(:,:,:))**(-2)*Telec(:,:,:)**(-3)*&
            gfac(i)/gfac(6)*Nelec**2*1.0e-6*1.0e-6
    enddo


    Gm_ion(:,:,:)=0.d0
    Gm_rec(:,:,:)=0.d0

    nntot=0.d0

    do i=1,5
        !Recombination rate
        Gm_rec(:,:,:)=Gm_rec(:,:,:)+Nexcite(:,:,:,6)*colrat(:,:,:,6,i)
        !Ionisation rate
        Gm_ion(:,:,:)=Gm_ion(:,:,:)+Nexcite(:,:,:,i)*colrat(:,:,:,i,6)
        nntot=nntot+Nexcite(:,:,:,i)
    enddo

    !Need to divide by the species density to make consistent with the rest of the code
    Gm_ion=Gm_ion/nntot
    Gm_rec=Gm_rec/Nexcite(:,:,:,6)

!print*,'NOT CODED UP YET'
!stop

  end subroutine get_col_ion_coeff_aprox
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
    double precision::rn(6),bn(6),xrat(5,5),Enn(5,5),yhat,rnn(5,5),zhat,gauntfac(5,5),fnn(5,5)
    double precision::Ann(5,5),Bnn(5,5),E0y,E1y,E2y,E0z,E1z,E2z
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

    rn(1)=0.45d0 !Equation 31
    rn(2)=1.94d0*2.d0**(-1.57d0) !Equation 32
    rn(3)=1.94d0*3.d0**(-1.57d0) !Equation 32
    rn(4)=1.94d0*4.d0**(-1.57d0) !Equation 32
    rn(5)=1.94d0*5.d0**(-1.57d0) !Equation 32
    rn(6)=1.94d0*6.d0**(-1.57d0) !Equation 32

    bn(1)=-0.603d0 !Equation 25
    bn(2)=(1.0d0/2.d0)*(4.0d0-18.63d0/2.d0+36.24d0/(2.d0**2)-&
                    28.09d0/(2.d0**3)) !Equation 26
    bn(3)=(1.0d0/3.d0)*(4.0d0-18.63d0/3.d0+36.24d0/(3.d0**2)-&
                    28.09d0/(3.d0**3)) !Equation 26
    bn(4)=(1.0d0/4.d0)*(4.0d0-18.63d0/4.d0+36.24d0/(4.d0**2)-&
                    28.09d0/(4.d0**3)) !Equation 26
    bn(5)=(1.0d0/5.d0)*(4.0d0-18.63d0/5.d0+36.24d0/(5.d0**2)-&
                    28.09d0/(5.d0**3)) !Equation 26
    bn(6)=(1.0d0/6.d0)*(4.0d0-18.63d0/6.d0+36.24d0/(6.d0**2)-&
                    28.09d0/(6.d0**3)) !Equation 26

    do ii=1,5; do jj=ii+1,5
        xrat(ii,jj)=1.0d0-(dble(ii)/dble(jj))**2 !ratio of transition energy to ionisation energy of lower level
        Enn(ii,jj)=Eion(ii)-Eion(jj) !Difference in ionisation energies
        rnn(ii,jj)=rn(ii)*xrat(ii,jj) !Equation 30

        !Gaunt factors from table 1 and equation 4
        !Constants so can be moved for speed
        if (ii .eq. 1) then 
            gauntfac(ii,jj)=1.1330d0-0.4059d0/xrat(ii,jj)+0.07014d0/(xrat(ii,jj)**2)
        else if (ii .eq. 2) then 
            gauntfac(ii,jj)=1.0785d0-0.2319d0/xrat(ii,jj)+0.02947d0/(xrat(ii,jj)**2)
        else
            gauntfac(ii,jj)=0.9935d0+0.2328d0/dble(ii)-0.1296d0/(dble(ii)**2)&
                -(1.0d0/xrat(ii,jj))*(1.0d0/dble(ii))*(0.6282d0-0.5598d0/dble(ii)+0.5299d0/(dble(ii)**2))&
                +(1.0d0/xrat(ii,jj))**2*(1.0d0/dble(ii)**2)*(0.3887d0-1.181d0/dble(ii)+1.470d0/(dble(ii)**2))
        endif

        fnn(ii,jj)=32.0d0/3.0d0/dsqrt(3.d0)/pi*dble(ii)/dble(jj)**3/(xrat(ii,jj)**3)*gauntfac(ii,jj) !Equation 3

        Ann(ii,jj)=2.0d0*dble(ii)**2/xrat(ii,jj)*fnn(ii,jj) !Equation 11
        Bnn(ii,jj)=4.0d0*dble(ii)**4/(dble(jj)**3)/(xrat(ii,jj)**2)*(1.0d0+4.0d0/3.0d0/xrat(ii,jj)+bn(ii)/(xrat(ii,jj)**2)) !Equation 23
    enddo;enddo

    !set colex to zero initially
    colex(:,:)=0.d0
    !Loop over the grid
    do k=1,kx;do j=1,jx; do i=1,ix
        !loop over the Excitation states
            do ii=1,5

                do jj=ii+1,5

                    yhat=Enn(ii,jj)/kboltz/Telec(i,j,k) !Equation 37

                    zhat=rnn(ii,jj)+Enn(ii,jj)/kboltz/Telec(i,j,k)   !Equation 38

                    !Exponential fits from a table
                    !Linear interpolation
!                    call ionexpfitinterp1(yhat,1,E1y)
!                    call ionexpfitinterp1(zhat,1,E1z)
!                    call ionexpfitinterp1(yhat,2,E2y)
!                    call ionexpfitinterp1(zhat,2,E2z)
                    !Quadratic interpolation
                    call ionexpfitinterp2(yhat,1,E1y)
                    call ionexpfitinterp2(zhat,1,E1z)
                    call ionexpfitinterp2(yhat,2,E2y)
                    call ionexpfitinterp2(zhat,2,E2z)

                    !Exponential fits (THESE NEED WORK)
                    !call ionexpfittest(yhat,1.d0,0.001d0,0.0001d0,E1y)   !Equation 8
                    !call ionexpfittest(zhat,1.d0,0.001d0,0.0001d0,E1z)
                    !call ionexpfittest(yhat,2.d0,0.001d0,0.0001d0,E2y)
                    !call ionexpfittest(zhat,2.d0,0.001d0,0.0001d0,E2z)

!                    print*,E1y,E1z,E2y,E2z

                    !Rate coefficient for excitation Equation 36 using electron mass
                    Colex(ii,jj)=dsqrt(8.d0*kboltz*Telec(i,j,k)/pi/melec)*2.0d0*&
                        dble(ii)**2/xrat(ii,jj)*pi*a0bohr**2*yhat**2*&
                              (Ann(ii,jj)*((1.0d0/yhat+0.5d0)*E1y-(1.0d0/zhat+0.5d0)*E1z)+&
                              (Bnn(ii,jj)-Ann(ii,jj)*dlog(2.d0*dble(ii)**2/xrat(ii,jj)))*(E2y/yhat-E2z/zhat))

                enddo

                !now do the ionisation states

                !Ionisation coefficient
                yn=Eion(ii)/kboltz/Telec(i,j,k) 
                zn=rn(ii)+Eion(ii)/kboltz/Telec(i,j,k)


                !Exponential fits from a table
                !Linear interpolation
!                call ionexpfitinterp1(yn,1,E1y)
!                call ionexpfitinterp1(zn,1,E1z)
!                call ionexpfitinterp1(yn,2,E2y)
!                call ionexpfitinterp1(zn,2,E2z)
!                call ionexpfitinterp1(yn,0,E0y)
!                call ionexpfitinterp1(zn,0,E0z)

                !Quadratic interpolation
                call ionexpfitinterp2(yn,1,E1y)
                call ionexpfitinterp2(zn,1,E1z)
                call ionexpfitinterp2(yn,2,E2y)
                call ionexpfitinterp2(zn,2,E2z)
                call ionexpfitinterp2(yn,0,E0y)
                call ionexpfitinterp2(zn,0,E0z)

                !Call the hokey exponential integral
                !call ionexpfittest(yn,1.d0,0.001d0,0.0001d0,E1y)   !Equation 8
                !call ionexpfittest(zn,1.d0,0.001d0,0.0001d0,E1z)
                !call ionexpfittest(yn,2.d0,0.001d0,0.0001d0,E2y)
                !call ionexpfittest(zn,2.d0,0.001d0,0.0001d0,E2z)
                !call ionexpfittest(yn,0.d0,0.001d0,0.0001d0,E0y)
                !call ionexpfittest(zn,0.d0,0.001d0,0.0001d0,E0z)

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
                Bn0=2.0d0/3.0d0*dble(ii)**2*(5.0d0+bn(ii)) !Equation 24

                !Ionisation coefficients Equation 35 using electron mass
                Colex(ii,6)=dsqrt(8.d0*kboltz*Telec(i,j,k)/pi/melec)&
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
                        colrat(i,j,k,ii,jj)=nelec(i,j,k)*1.0e-6*gweight(ii)/gweight(jj)*Colex(ii,jj)*dsqrt(Telec(i,j,k)) 
                    else
                        colrat(i,j,k,ii,jj)=nelec(i,j,k)*1.0e-6*Colex(jj,ii)*dsqrt(Telec(i,j,k))*&
                            exp(-(Eion(jj)-Eion(ii))/kboltz/Telec(i,j,k))
                    endif
                enddo
                !Ionisation rates
                colrat(i,j,k,ii,6)=nelec(i,j,k)*1.0e-6*Colex(ii,6)*dsqrt(Telec(i,j,k))*exp(-(Eion(6)-Eion(ii))/kboltz/Telec(i,j,k))
            enddo

            do ii=1,5
                colrat(i,j,k,6,ii)=nelec(i,j,k)*1.0e-6*nelec(i,j,k)*1.0e-6*colrat(i,j,k,ii,6)*gweight(ii)/gweight(6)&
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
                Gm_ion(i,j,k)=Gm_ion(i,j,k)+max(Nexcite(i,j,k,ii)*colrat(i,j,k,ii,6),0.d0)
                Gm_rec(i,j,k)=Gm_rec(i,j,k)+max(Nexcite(i,j,k,6)*colrat(i,j,k,6,ii),0.d0)
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
  subroutine get_radrat_fixed(Trad,Telec,nelec,Gm_ion_rad,Gm_rec_rad)
!use the fixed ratiative rates from Sollum 2003 thesis
    double precision,intent(in)::Trad,Telec(ix,jx,kx),nelec(ix,jx,kx)
    double precision,intent(out)::Gm_ion_rad(ix,jx,kx),Gm_rec_rad(ix,jx,kx)
    double precision::nuarr(6,6),fosc(6,6),gfac(6),E(6),dneut(6)
    double precision::nu0,sol,oldsol,s1c,alp0,diff,exf,sahasol
    double precision::dntot
    integer::i,j,k,ii,jj,iii
    double precision,parameter::h=6.62607004e-34 !Planck's constant in m2 kg s^-1
    double precision,parameter::cli=299792458.d0 !Speed of light in m/s
    double precision,parameter::ech=-1.6e-19 !Charge of electron in coulomb.
    double precision,parameter::melec=9.10938356e-31 !Electron mass [kg]
    double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
    integer::nmaxloc
    !Oscilator strengths
    !From table 1 in Goldwire 1968
    !http://articles.adsabs.harvard.edu//full/1968ApJS...17..445G/0000450.000.html
    fosc(:,:)=0.0
    fosc(1,2)=4.1619672e-1
    fosc(1,3)=7.9101563e-2
    fosc(1,4)=2.8991029e-2
    fosc(1,5)=1.3938344e-2
    fosc(2,3)=6.4074704e-1
    fosc(2,4)=1.1932114e-1
    fosc(2,5)=4.4670295e-2
    fosc(3,4)=8.4209639e-1
    fosc(3,5)=1.5058408e-1
    fosc(4,5)=1.0377363e0

    !Transition frequenices [Hz]
    !Lyman series
    nuarr(1,2)=cli/121.57e-9
    nuarr(1,3)=cli/102.57e-9
    nuarr(1,4)=cli/97.254e-9
    nuarr(1,5)=cli/94.973e-9
    nuarr(1,6)=cli/91.175e-9 !Continuum
    !Balmer series
    nuarr(2,3)=cli/656.3e-9
    nuarr(2,4)=cli/486.1e-9
    nuarr(2,5)=cli/434.0e-9
    nuarr(2,6)=cli/364.6e-9 !Continuum
    !Pachen series
    nuarr(3,4)=cli/1875.0e-9
    nuarr(3,5)=cli/1282.0e-9
    nuarr(3,6)=cli/820.4e-9 !Continuum
    !Brackett series
    nuarr(4,5)=cli/4051.0e-9
    nuarr(4,6)=cli/1458.0e-9 !Continuum
    !Pfund series
    nuarr(5,6)=cli/2279.0e-9 !Continuum
    nu0=nuarr(1,6)

    !Gaunt factors
    gfac(1)=2.d0*1.d0**2 !ground state
    gfac(2)=2.d0*2.d0**2 !1st level excitation
    gfac(3)=2.d0*3.d0**2 !2nd level excitation
    gfac(4)=2.d0*4.d0**2 !3rd level excitation
    gfac(5)=2.d0*5.d0**2 !4th level excitation
    gfac(6)=1.d0 !Ionised hydrogen

    !Reverse engineered cross section at level 1
    s1c=2.815e29*1.0**4/1.0**5/nuarr(1,6)**3*gfac(1) !Equation 3.4 in Sollum
    alp0=s1c*((nu0/nuarr(1,6))**(-3)) !in cm^2
    alp0=alp0*1.0e-4 !in m^2

    E=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
    E=E*1.6022e-19 !Convert to joules (to be dimensionally correct)

    radrat(:,:,:,:,:)=0.d0
    !Excitation states
        do ii=1,5 
            do jj=ii+1,5
                radrat(:,:,:,ii,jj)=(4.d0*pi/h/nu0)*(pi*ech**2/melec/cli)*fosc(ii,jj)&
                    *(2.d0*h*nu0**3/cli/cli)/(dexp(h*nu0/kboltz/Trad)-1.d0)
            enddo
        enddo
    !DeExcitation rates (eq 3.3 in Sollum)
        do jj=2,5
            do ii=jj-1,1,-1
                radrat(:,:,:,jj,ii)=radrat(i,j,k,ii,jj)*gfac(ii)/gfac(jj)*exp(h*nu0/kboltz/Trad)
            enddo
        enddo

    !Ionisation rates (Sollum eq 3.5)
    do ii=1,5
        !Evaluate the integral
        sol=0.d0
        oldsol=sol
        diff=1.d0
        iii=1
        do while (diff .GT. 0.0001)
            oldsol=sol
            call expintruttonn1(iii*h*nuarr(ii,6)/kboltz/Trad,exf)
            sol=sol+exf
            diff=abs((sol-oldsol)/sol)
            iii=iii+1
        enddo

        radrat(:,:,:,ii,6)=8.d0*pi*alp0*(1.d0*nu0)**3/cli/cli*sol
    enddo

!print*,'ONLY DONE UP TO HERE!'
    !Recombination rates
    do ii=1,5
        !Evaluate the integral
        sol=0.d0
        oldsol=sol
        diff=1.0d0
        iii=1

        do k=1,kx;do j=1,jx; do i=1,ix
        do while (diff .GT. 0.0001)
            oldsol=sol
            call expintruttonn1(1.d0*iii*h*nu0/kboltz/Trad+1.d0*h*nu0/kboltz/Telec(i,j,k),exf)
            sol=sol+exf
            diff=abs((sol-oldsol)/sol)
            iii=iii+1
        enddo

        sahasol=(2.d0/nelec(i,j,k)*gfac(6)/gfac(ii)*(2.d0*pi*melec*kboltz*Telec(i,j,k)/h/h)**(3.0/2.0)*&
                dexp(-E(ii)/kboltz/Telec(i,j,k)))

        radrat(i,j,k,6,ii)=8.d0*pi*alp0*(1.d0*nu0)**3/cli/cli*sol*sahasol
!print*,radrat(1,1,1,6,*)
        enddo;enddo;enddo
    enddo

    Gm_ion_rad(:,:,:)=0.d0
    Gm_rec_rad(:,:,:)=0.d0
    do ii=1,5
        Gm_ion_rad(:,:,:)=Gm_ion_rad(:,:,:)+max(Nexcite(:,:,:,ii)*radrat(:,:,:,ii,6),0.d0)
        Gm_rec_rad(:,:,:)=Gm_rec_rad(:,:,:)+max(Nexcite(:,:,:,6)*radrat(:,:,:,6,ii),0.d0)
    enddo
!print*,gm_ion(1,1,1),gm_rec(1,1,1)
!stop
    !divide rates by total neutral/plasma density for consistency
    Gm_ion_rad(:,:,:)=Gm_ion_rad(:,:,:)/(Nexcite(:,:,:,1)+Nexcite(:,:,:,2)+Nexcite(:,:,:,3)&
                  +Nexcite(:,:,:,4)+Nexcite(:,:,:,5))        
    Gm_rec_rad(:,:,:)=Gm_rec_rad(:,:,:)/Nexcite(:,:,:,6)

  end subroutine get_radrat_fixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine get_radrat_fixed_slow(Trad,Telec,nelec,Gm_ion_rad,Gm_rec_rad)
!use the fixed ratiative rates from Sollum 2003 thesis
    double precision,intent(in)::Trad,Telec(ix,jx,kx),nelec(ix,jx,kx)
    double precision,intent(out)::Gm_ion_rad(ix,jx,kx),Gm_rec_rad(ix,jx,kx)
    double precision::nuarr(6,6),fosc(6,6),gfac(6),E(6),dneut(6)
    double precision::nu0,sol,oldsol,s1c,alp0,diff,exf,sahasol
    double precision::dntot
    integer::i,j,k,ii,jj,iii
    double precision,parameter::h=6.62607004e-34 !Planck's constant in m2 kg s^-1
    double precision,parameter::cli=299792458.d0 !Speed of light in m/s
    double precision,parameter::ech=-1.6e-19 !Charge of electron in coulomb.
    double precision,parameter::melec=9.10938356e-31 !Electron mass [kg]
    double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
    integer::nmaxloc
    !Oscilator strengths
    !From table 1 in Goldwire 1968
    !http://articles.adsabs.harvard.edu//full/1968ApJS...17..445G/0000450.000.html
    fosc(:,:)=0.0
    fosc(1,2)=4.1619672e-1
    fosc(1,3)=7.9101563e-2
    fosc(1,4)=2.8991029e-2
    fosc(1,5)=1.3938344e-2
    fosc(2,3)=6.4074704e-1
    fosc(2,4)=1.1932114e-1
    fosc(2,5)=4.4670295e-2
    fosc(3,4)=8.4209639e-1
    fosc(3,5)=1.5058408e-1
    fosc(4,5)=1.0377363e0

    !Transition frequenices [Hz]
    !Lyman series
    nuarr(1,2)=cli/121.57e-9
    nuarr(1,3)=cli/102.57e-9
    nuarr(1,4)=cli/97.254e-9
    nuarr(1,5)=cli/94.973e-9
    nuarr(1,6)=cli/91.175e-9 !Continuum
    !Balmer series
    nuarr(2,3)=cli/656.3e-9
    nuarr(2,4)=cli/486.1e-9
    nuarr(2,5)=cli/434.0e-9
    nuarr(2,6)=cli/364.6e-9 !Continuum
    !Pachen series
    nuarr(3,4)=cli/1875.0e-9
    nuarr(3,5)=cli/1282.0e-9
    nuarr(3,6)=cli/820.4e-9 !Continuum
    !Brackett series
    nuarr(4,5)=cli/4051.0e-9
    nuarr(4,6)=cli/1458.0e-9 !Continuum
    !Pfund series
    nuarr(5,6)=cli/2279.0e-9 !Continuum
    nu0=nuarr(1,6)

    !Gaunt factors
    gfac(1)=2.d0*1.d0**2 !ground state
    gfac(2)=2.d0*2.d0**2 !1st level excitation
    gfac(3)=2.d0*3.d0**2 !2nd level excitation
    gfac(4)=2.d0*4.d0**2 !3rd level excitation
    gfac(5)=2.d0*5.d0**2 !4th level excitation
    gfac(6)=1.d0 !Ionised hydrogen

    !Reverse engineered cross section at level 1
    s1c=2.815e29*1.0**4/1.0**5/nuarr(1,6)**3*gfac(1) !Equation 3.4 in Sollum
    alp0=s1c*((nu0/nuarr(1,6))**(-3)) !in cm^2
    alp0=alp0*1.0e-4 !in m^2

    E=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
    E=E*1.6022e-19 !Convert to joules (to be dimensionally correct)

    !Loop over grid
    do k=1,kx;do j=1,jx; do i=1,ix
    radrat(i,j,k,:,:)=0.d0
    !Excitation states
        do ii=1,5 
            do jj=ii+1,5
                radrat(i,j,k,ii,jj)=(4.d0*pi/h/nu0)*(pi*ech**2/melec/cli)*fosc(ii,jj)&
                    *(2.d0*h*nu0**3/cli/cli)/(dexp(h*nu0/kboltz/Trad)-1.d0)
            enddo
        enddo
    !DeExcitation rates (eq 3.3 in Sollum)
        do jj=2,5
            do ii=jj-1,1,-1
                radrat(i,j,k,jj,ii)=radrat(i,j,k,ii,jj)*gfac(ii)/gfac(jj)*exp(h*nu0/kboltz/Trad)
            enddo
        enddo

    !Ionisation rates (Sollum eq 3.5)
    do ii=1,5
        !Evaluate the integral
        sol=0.d0
        oldsol=sol
        diff=1.d0
        iii=1
        do while (diff .GT. 0.0001)
            oldsol=sol
            call expintruttonn1(iii*h*nuarr(ii,6)/kboltz/Trad,exf)
            sol=sol+exf
            diff=abs((sol-oldsol)/sol)
            iii=iii+1
        enddo

        radrat(i,j,k,ii,6)=8.d0*pi*alp0*(1.d0*nu0)**3/cli/cli*sol
    enddo

    !Recombination rates
    do ii=1,5
        !Evaluate the integral
        sol=0.d0
        oldsol=sol
        diff=1.0d0
        iii=1
        do while (diff .GT. 0.0001)
            oldsol=sol
            call expintruttonn1(1.d0*iii*h*nu0/kboltz/Trad+1.d0*h*nu0/kboltz/T0,exf)
            sol=sol+exf
            diff=abs((sol-oldsol)/sol)
            iii=iii+1
        enddo

        sahasol=(2.d0/nelec(i,j,k)*gfac(6)/gfac(j)*(2.d0*pi*melec*kboltz*Telec(i,j,k)/h/h)**(3.0/2.0)*&
                dexp(-E(jj)/kboltz/Telec(i,j,k)))

        radrat(i,j,k,6,jj)=8.d0*pi*alp0*(1.d0*nu0)**3/cli/cli*sol*sahasol
    enddo

    ! Calculate the change in each species
    dntot=0.d0
    do ii=1,6
        do jj=1,6 
            dneut(ii)=Nexcite(i,j,k,jj)*radrat(i,j,k,jj,ii) - Nexcite(i,j,k,ii)*radrat(i,j,k,ii,jj)
        enddo

            if (Nexcite(i,j,k,ii) .ne. maxval(Nexcite(i,j,k,:))) then 
                dntot=dntot+dneut(ii)
            endif
    enddo
    nmaxloc=maxloc(Nexcite(i,j,k,:),DIM=1)
    dneut(nmaxloc)=-dntot

    Gm_ion_rad(i,j,k)=0.d0
    Gm_rec_rad(i,j,k)=0.d0
    do ii=1,5
        Gm_ion_rad(i,j,k)=Gm_ion_rad(i,j,k)+max(Nexcite(i,j,k,ii)*radrat(i,j,k,ii,6),0.d0)
        Gm_rec_rad(i,j,k)=Gm_rec_rad(i,j,k)+max(Nexcite(i,j,k,6)*radrat(i,j,k,6,ii),0.d0)
    enddo
!print*,gm_ion(1,1,1),gm_rec(1,1,1)
!stop
    !divide rates by total neutral/plasma density for consistency
    Gm_ion_rad(i,j,k)=Gm_ion_rad(i,j,k)/(Nexcite(i,j,k,1)+Nexcite(i,j,k,2)+Nexcite(i,j,k,3)&
                  +Nexcite(i,j,k,4)+Nexcite(i,j,k,5))        
    Gm_rec_rad(i,j,k)=Gm_rec_rad(i,j,k)/Nexcite(i,j,k,6)

    enddo;enddo;enddo
  end subroutine get_radrat_fixed_slow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ionexpfitinterp1(zhat,i,sol)
    !Interpolate the exponential integtal using the table of exponential fits 
    !Linear interpolation
    double precision,intent(in)::zhat
    integer,intent(in)::i
    double precision,intent(out)::sol
    double precision::xmin,xmax,xrange,weight
    integer::xindex

    xmin=expinttab(1,1)
    xmax=expinttab(1,size(expinttab(1,:)))
    xrange=xmax-xmin

    xindex=floor((zhat-xmin)/xrange*size(expinttab(1,:)))

    if (xindex .le. 0) xindex=1
    if (xindex .ge. size(expinttab(1,:))) xindex=size(expinttab(1,:))-1

    weight=(zhat - expinttab(1,xindex))/(expinttab(1,xindex+1)-expinttab(1,xindex))
    sol = (1.0-weight)*expinttab(i+2,xindex) + weight*expinttab(i+2,xindex+1);
!print*,xindex,xmin,zhat,xrange,xmax
!print*,expinttab(1,xindex),zhat,expinttab(1,xindex+1)
!print*,expinttab(i+2,xindex),sol,expinttab(i+2,xindex+1)
!stop
  end subroutine ionexpfitinterp1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ionexpfitinterp2(zhat,i,sol)
    !Interpolate the exponential integtal using the table of exponential fits 
    !Quadratic interpolation
    !WORK IN PROGRESS
    double precision,intent(in)::zhat
    integer,intent(in)::i
    double precision,intent(out)::sol
    double precision::xmin,xmax,xrange,weight0,weight1,weight2
    integer::xindex,xindex0,xindex1,xindex2

    xmin=expinttab(1,1)
    xmax=expinttab(1,size(expinttab(1,:)))
    xrange=xmax-xmin

    xindex=floor((zhat-xmin)/xrange*size(expinttab(1,:)))

    if (xindex .le. 0) xindex=1
    if (xindex .ge. size(expinttab(1,:))) xindex=size(expinttab(1,:))-1

    xindex0=xindex
    xindex1=xindex+1
    xindex2=xindex+2

    weight0=((zhat                -expinttab(1,xindex1))*(zhat                - expinttab(1,xindex2)))/&
            ((expinttab(1,xindex0)-expinttab(1,xindex1))*(expinttab(1,xindex0)- expinttab(1,xindex2)))
    weight1=((zhat                -expinttab(1,xindex0))*(zhat                - expinttab(1,xindex2)))/&
            ((expinttab(1,xindex1)-expinttab(1,xindex0))*(expinttab(1,xindex1)- expinttab(1,xindex2)))
    weight2=((zhat                -expinttab(1,xindex0))*(zhat                - expinttab(1,xindex1)))/&
            ((expinttab(1,xindex2)-expinttab(1,xindex0))*(expinttab(1,xindex2)- expinttab(1,xindex1)))

    
    sol = weight0*expinttab(i+2,xindex0)+weight1*expinttab(i+2,xindex1)+weight2*expinttab(i+2,xindex2)
!print*,xindex,xmin,zhat,xrange,xmax
!print*,expinttab(1,xindex),zhat,expinttab(1,xindex+1)
!print*,expinttab(i+2,xindex),sol,expinttab(i+2,xindex+1)
!stop
  end subroutine ionexpfitinterp2

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
subroutine expintruttonn1(x,sol)
!Exponential integral as deffined by Rutton 2003, page 78
!Abramowitz and Stegun(1964), solution 5.1.11 (results pair with Rutton 1999 table 4.1)
    double precision,intent(in)::x
    double precision,intent(out)::sol
    double precision::oldsol,diff
    integer::i
    if (x .LT. 9.5) then 
        sol=-0.5772155-dlog(x)  !constant is Eulars constant
        oldsol=sol
        diff=1.0
        i=1
        do while (diff .GT. 0.0001)
            oldsol=sol
            sol=sol-(-1.d0)**i*x**i/(dble(i)*dgamma(dble(i+1))) !gamma(i+1)=i! I think
            diff=abs((sol-oldsol)/sol)
            i=i+1
        enddo
    else
        sol=dexp(-x)/x !Large value approximation (Rutton 2003, just below eqn 4.13)
    endif

END subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine hydrogen_excitation_update(dt,rom,roh)
! update the hydrogen excitation states
  double precision,intent(in)::dt,rom(ix,jx,kx),roh(ix,jx,kx)
  double precision::dneut(6)
  double precision::dntot
  double precision::dneutv(ix,jx,kx,6)
  double precision::dntotv(ix,jx,kx)
  integer::i,j,k,ii,jj,nmaxloc
!Vector form
!The new number of electrons/protons is:
    Nexcite(:,:,:,6)=rom
!    do k=1,kx;do j=1,jx; do i=1,ix
            ! Calculate the change in each neutral species
            dntotv(:,:,:)=0.d0
            do ii=1,5
                do jj=1,6 
                    dneutv(:,:,:,ii)=Nexcite(:,:,:,jj)*colrat(:,:,:,jj,ii) - Nexcite(:,:,:,ii)*colrat(:,:,:,ii,jj)
                    if (flag_rad .eq. 1) then
                        dneutv(:,:,:,ii)=Nexcite(:,:,:,jj)*radrat(:,:,:,jj,ii) - Nexcite(:,:,:,ii)*radrat(:,:,:,ii,jj)
                    endif
                enddo
            enddo
!print*,dneut
!print*,colrat(i,j,k,:,:)
!print*,Nexcite(i,j,k,:)
            Nexcite(:,:,:,1:5)=Nexcite(:,:,:,1:5)+dt*dneutv(:,:,:,1:5)/n0
            where (Nexcite>=0.d0)
                Nexcite = Nexcite
            elsewhere
                Nexcite = 0.0d0!1.e-16
            end where
!print*,Nexcite(i,j,k,:)
            dntotv(:,:,:)=sum(Nexcite(:,:,:,1:5),DIM=4)
            Nexcite(:,:,:,1)=roh(:,:,:)*Nexcite(:,:,:,1)/dntotv(:,:,:)
            Nexcite(:,:,:,2)=roh(:,:,:)*Nexcite(:,:,:,2)/dntotv(:,:,:)
            Nexcite(:,:,:,3)=roh(:,:,:)*Nexcite(:,:,:,3)/dntotv(:,:,:)
            Nexcite(:,:,:,4)=roh(:,:,:)*Nexcite(:,:,:,4)/dntotv(:,:,:)
            Nexcite(:,:,:,5)=roh(:,:,:)*Nexcite(:,:,:,5)/dntotv(:,:,:)

!    enddo;enddo;enddo
!print*,Nexcite(1,1,1,6),rom(1,1,1),sum(Nexcite(1,1,1,1:5)),roh(1,1,1)
  end subroutine hydrogen_excitation_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IRgetionpot(nde,enloss)
! calculate the ionisation potential loss
  double precision,intent(in)::nde(ix,jx,kx)
  double precision,intent(out)::enloss(ix,jx,kx)
  double precision::ieloss(ix,jx,kx),Eev(6)

Eev=[13.6,3.4,1.51,0.85,0.54,0.0]

!Ionisation potential energy
ieloss(:,:,:)=colrat(:,:,:,1,6)*13.6*nexcite(:,:,:,1)+&
        colrat(:,:,:,2,6)*3.4*nexcite(:,:,:,2)+&
        colrat(:,:,:,3,6)*1.51*nexcite(:,:,:,3)+&
        colrat(:,:,:,4,6)*0.85*nexcite(:,:,:,4)+&
        colrat(:,:,:,5,6)*0.54*nexcite(:,:,:,5)

!Excitation energy loss
ieloss=ieloss+nexcite(:,:,:,1)*(colrat(:,:,:,1,2)*(Eev(1)-Eev(2))+colrat(:,:,:,1,3)*(Eev(1)-Eev(3))+&
                    colrat(:,:,:,1,4)*(Eev(1)-Eev(4))+colrat(:,:,:,1,5)*(Eev(1)-Eev(5)))+&
              nexcite(:,:,:,2)*(colrat(:,:,:,2,3)*(Eev(2)-Eev(3))+colrat(:,:,:,2,4)*(Eev(2)-Eev(4))+&
                    colrat(:,:,:,2,5)*(Eev(2)-Eev(5)))+&
              nexcite(:,:,:,3)*(colrat(:,:,:,3,4)*(Eev(3)-Eev(4))+colrat(:,:,:,3,5)*(Eev(3)-Eev(5)))+&
              nexcite(:,:,:,4)*(colrat(:,:,:,4,5)*(Eev(4)-Eev(5)))

    if(mod(flag_col,2) .eq. 1) then
!	    enloss=nde/gm/T0/8.6173e-5*(&
!                (13.6d0*max(Nexcite(:,:,:,1)*colrat(:,:,:,1,6),0.d0))+& !ground state
!                ( 3.4d0*max(Nexcite(:,:,:,2)*colrat(:,:,:,2,6),0.d0))+& !1st excited state
!                (1.51d0*max(Nexcite(:,:,:,3)*colrat(:,:,:,3,6),0.d0))+& !2nd excited  state
!                (0.85d0*max(Nexcite(:,:,:,4)*colrat(:,:,:,4,6),0.d0))+& !3rd excited  state
!                (0.54d0*max(Nexcite(:,:,:,5)*colrat(:,:,:,5,6),0.d0)))/& !4th excited  state
!                (Nexcite(:,:,:,1)+Nexcite(:,:,:,2)+Nexcite(:,:,:,3)+Nexcite(:,:,:,4)+Nexcite(:,:,:,5)) 
	    enloss=nde/gm/T0/8.6173e-5*ieloss/&
                (Nexcite(:,:,:,1)+Nexcite(:,:,:,2)+Nexcite(:,:,:,3)+Nexcite(:,:,:,4)+Nexcite(:,:,:,5)) 
    elseif(mod(flag_col,2) .eq. 0) then
!	    enloss=nde*(beta/T0/2.d0/8.6173e-5)*(&
!                (13.6d0*max(Nexcite(:,:,:,1)*colrat(:,:,:,1,6),0.d0))+& !ground state
!                ( 3.4d0*max(Nexcite(:,:,:,2)*colrat(:,:,:,2,6),0.d0))+& !1st excited state
!                (1.51d0*max(Nexcite(:,:,:,3)*colrat(:,:,:,3,6),0.d0))+& !2nd excited  state
!                (0.85d0*max(Nexcite(:,:,:,4)*colrat(:,:,:,4,6),0.d0))+& !3rd excited  state
!                (0.54d0*max(Nexcite(:,:,:,5)*colrat(:,:,:,5,6),0.d0)))/& !4th excited  state
!                (Nexcite(:,:,:,1)+Nexcite(:,:,:,2)+Nexcite(:,:,:,3)+Nexcite(:,:,:,4)+Nexcite(:,:,:,5))
	    enloss=nde*(beta/T0/2.d0/8.6173e-5)*ieloss/& 
                (Nexcite(:,:,:,1)+Nexcite(:,:,:,2)+Nexcite(:,:,:,3)+Nexcite(:,:,:,4)+Nexcite(:,:,:,5))  
    else
	    print*,'option not included!'
	    stop
    endif
!print*,'EXCITATION NEEDS TO BE ADDED'
    !Normalise based on parameters    
    enloss=enloss/Gm_rec_ref*t_ir

  end subroutine IRgetionpot
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
        if (IR_type .eq. 4) then
            call irgetionpot(nde,ion_pot)
        else
		    if(mod(flag_col,2) .eq. 1) then
			    ion_pot=Gm_ion*nde*(13.6d0/gm/T0/8.6173e-5)
		    elseif(mod(flag_col,2) .eq. 0) then
			    ion_pot=Gm_ion*nde*(13.6d0*beta/T0/2.d0/8.6173e-5)
		    else
			    print*,'option not included!'
			    stop
		    endif
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
	endif
    !Radiative ionisation
    if (flag_rad .eq. 1) then
       ds(:,:,:,1)=Gm_rec_rad*de-Gm_ion_rad*nde
       ds(:,:,:,2)=Gm_rec_rad*de*vx-Gm_ion_rad*nde*nvx
       ds(:,:,:,3)=Gm_rec_rad*de*vy-Gm_ion_rad*nde*nvy
       ds(:,:,:,4)=Gm_rec_rad*de*vz-Gm_ion_rad*nde*nvz
	   ds(:,:,:,5)=0.5d0*(Gm_rec_rad*de*(vx*vx+vy*vy+vz*vz)- &
	        Gm_ion_rad*nde*(nvx*nvx+nvy*nvy+nvz*nvz)) -&
	        Gm_ion_rad*npr/(gm-1.d0) + 0.5d0*Gm_rec_rad*pr/(gm-1.d0)
	   S_h(:,:,:,1:5)=S_h(:,:,:,1:5)+ds(:,:,:,1:5)
	   S_m(:,:,:,1:5)=S_m(:,:,:,1:5)-ds(:,:,:,1:5)
    endif

    endif    

    return
  end subroutine source_PIP

end module PIP_rot
