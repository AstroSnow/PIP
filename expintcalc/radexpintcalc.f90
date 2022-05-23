PROGRAM radexpintcalc
	IMPLICIT NONE
	double precision::sol,oldsol,diff,Trad,sol2,Telec
	double precision::nuarr,exf,dion,ionval
	double precision,parameter::cli=299792458.d0 !Speed of light in m/s
	double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
	double precision,parameter::h=6.62607004e-34 !Planck's constant in m2 kg s^-1
	integer::iii,ti,nsamps
	
	! output data into a file 
    open(1, file = 'radexp.dat', status = 'replace')  
   
    Trad=5000.0
    
    ionmin=cli/121.57e-9*h/kboltz/100000.d0
    ionmax=cli/2279.0e-9*h/kboltz/100.d0
    nsamps=1000

	dion=(ionmax-ionmin)/(nsamps)
    
	do ti=0,nsamps
! 		Trad=ti*100.d0
! 		Telec=Trad*10.0
	 	!nuarr=cli/121.57e-9
!	 	nuarr=cli/2279.0e-9

		ionval=ionmin+ti*dion
			 
	 
	 	!IONISATION VALUES
		sol=0.d0
		oldsol=sol
		diff=1.d0
		iii=1
		do while (diff .GT. 0.0000001)
			oldsol=sol
!			call expintrutton1(dble(iii)*h*nuarr/kboltz/Trad,exf)
			call expintrutton1(dble(iii)*ionval,exf)
			sol=sol+exf
			diff=abs((sol-oldsol)/sol)
			iii=iii+1
		enddo


!WHAT ARE THE MAX AN MIN VALUES?
!USE EQUATION 3.12 in Sollum to simplify!	 
		!Recombination rates    
		sol2=0.d0
		call expintrutton1(h*nuarr/kboltz/Trad,sol2)
		oldsol=sol2
		diff=1.0d0
		iii=1
		do while (diff .GT. 0.0001)
		    oldsol=sol2
		    call expintrutton1(dble(iii)*h*nuarr/kboltz/Trad+h*nuarr/kboltz/Telec,exf)
		    sol2=sol2+exf
		    diff=abs((sol2-oldsol)/sol2)
		    iii=iii+1
		enddo
				
				
	   write(1,*) h*nuarr/kboltz/Trad, sol,h*nuarr/kboltz/Trad+h*nuarr/kboltz/Telec,sol2
	 enddo
	 
	 close(1) 

END PROGRAM radexpintcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine expintrutton1(x,sol)
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
        do while (diff .GT. 0.0000001)
            oldsol=sol
            sol=sol-(-1.d0)**i*x**i/(dble(i)*dgamma(dble(i+1))) !gamma(i+1)=i! I think
            diff=abs((sol-oldsol)/sol)
            i=i+1
        enddo
    else
        sol=dexp(-x)/x !Large value approximation (Rutton 2003, just below eqn 4.13)
    endif
ENDsubroutine expintrutton1
