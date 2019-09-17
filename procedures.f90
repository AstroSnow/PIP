module procedures
implicit none
contains


!=======================================================================================!
	function interpol(x,y,xint) result(yint)
!=======================================================================================!
		!Interpolation routine
		!1D
		implicit none
		double precision, intent(in) :: xint !Desired interpolation point
		double precision, dimension(:), intent(in) :: x,y !Initial arrays
		double precision :: yint !Interpolated value
		integer :: lx,ly, lxi, xex
		lx=size(x)
		ly=size(y)
		lxi=1
		xex=0
		do while ((x(lxi) .LE. xint) .AND. (lxi .LE. lx-1) .AND. (xex .EQ. 0))
		  	if (x(lxi) .EQ. xint) then !.OR. (abs(x(lxi)-xint) .LE. 1.0e-6)) then
				yint=y(lxi)
				xex=1
			else
				lxi=lxi+1
				yint=y(lxi-1)+(xint-x(lxi-1))*(y(lxi)-y(lxi-1))/(x(lxi)-x(lxi-1))
			endif
		enddo
	end function interpol
!=======================================================================================!


!=======================================================================================!
	function integrate(x,y,xint,h,order,locint) result(yint)
!=======================================================================================!
		!integrate y(x) at xint. Result is stored as yint.
		!11 IS BROKEN. DO NOT USE!!!!
		implicit none
		integer, intent(in) :: order, locint !linear=1
		double precision, dimension(:), intent(in) :: x,y !Initial arrays
		double precision, intent(in) :: xint,h !Height to integrate at
		double precision :: yint !Integrated value
		double precision :: yip1,yim1,yi1,yip2,yim2
		double precision :: yip3,yim3,yip4,yim4,yip5,yim5
		integer :: lx,lxi
		if (locint .EQ. 1) then
		!use edge-points
			yi1=interpol(x,y,xint-0.5d0*h)
			yip1=interpol(x,y,xint)
			yim1=interpol(x,y,xint-h)
			yip2=interpol(x,y,xint-0.25d0*h)
			yim2=interpol(x,y,xint-0.75d0*h)
		elseif (locint .EQ. 2) then
		!use mid-points
			yi1=interpol(x,y,xint)
			yip1=interpol(x,y,xint+0.5d0*h)
			yim1=interpol(x,y,xint-0.5d0*h)
			yip2=interpol(x,y,xint+0.25d0*h)
			yim2=interpol(x,y,xint-0.25d0*h)
		endif
		!Linear interpolation
		if (order .EQ. 1) then
			yint=0.5d0*h*(yip1+yim1)
		!Simpson's Rule
		elseif (order .EQ. 2) then
			yint=h/6.0d0*(yim1+yip1+4.d0*yi1)
		!Boole's Rule
		elseif (order .EQ. 3) then
			yint=2.d0/45.d0*h/4.d0*(7.d0*yim1+32.d0*yim2+12.d0*yi1+32.d0*yip2+7.d0*yip1)
		!11th order rule (BROKEN yim3 is cursed)
		elseif (order .EQ. 11) then
			yip5=interpol(x,y,xint+0.1d0*h)
			yip4=interpol(x,y,xint+0.2d0*h)
			yip3=interpol(x,y,xint+0.3d0*h)
			yip2=interpol(x,y,xint+0.4d0*h)
			yim5=interpol(x,y,xint-0.1d0*h)
			yim4=interpol(x,y,xint-0.2d0*h)
			yim3=interpol(x,y,xint-0.3d0*h) !This breaks it somehow
			yim2=interpol(x,y,xint-0.4d0*h)
			yint=16067.d0*(yip1+yim1)+106300.d0*(yip2+yim2)-48525.d0*(yip3+yim3)
			yint=yint+272400.d0*(yip4+yim4)-260550.d0*(yip5+yim5)+427368.d0*yi1
			yint=5.d0/299376.d0*h/10.d0*yint
		endif
	end function integrate

end module

