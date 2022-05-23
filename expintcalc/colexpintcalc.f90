PROGRAM radexpintcalc
	IMPLICIT NONE
	double precision::sol,oldsol,diff
	double precision::nuarr,exf,dy,yn,ymax,ymin,E0y,E1y,E2y
	double precision,parameter::cli=299792458.d0 !Speed of light in m/s
	double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
	double precision,parameter::h=6.62607004e-34 !Planck's constant in m2 kg s^-1
	integer::iii,ti,nsamps
	
	! output data into a file 
    open(1, file = 'colexp.dat', status = 'replace')  
   
    nsamps=100000
    
    ymin=0.54d0/13.6d0*2.18e-18/kboltz/100000.d0 !Assuming a maximum electron temperature of 100000 K, and 0.54 eV ionisation energy
    !ymax=0.45d0+2.18e-18/kboltz/100.d0 !Minimum simulation temperature 100 K,
	ymax=108.0 !Values are practially zero after this (10^-50).

	write(1,*) nsamps+1, ymin,ymax

	dy=(ymax-ymin)/(nsamps)
    
	do ti=0,nsamps

		yn=ymin+ti*dy

        call ionexpfittest(yn,0.d0,0.001d0,0.0001d0,E0y)   !Equation 8
        call ionexpfittest(yn,1.d0,0.001d0,0.0001d0,E1y)
        call ionexpfittest(yn,2.d0,0.001d0,0.0001d0,E2y)
				
				
	   write(1,*) yn, E0y,E1y,E2y
	 enddo
	 
	 close(1) 

END PROGRAM radexpintcalc
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



