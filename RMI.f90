subroutine RMI
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin,debug_direction
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
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
  double precision Atwood,mach,rcom,rpres,v_c(8),v_l(8),v_r(8)
  double precision b0,wtr
  integer i,j,k


!RMI initial condition from https://link.springer.com/article/10.1007/s00193-007-0104-z

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-10.0d0 ;end(1)=4.0d0
  start(2)= 0.0d0 ;end(2)=1.0d0
  start(3)=-0.001d0 ;end(3)=0.001d0
  call set_coordinate(start,end)
  !---------------------------------------

  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !atwood number
  atwood=3.0
  !atwood=1.0

  !Define shock states
	mach=2.d0
    B0=0.0d0
!B0=1.0d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)
    v_l=(/rcom,rpres/gm,-mach/rcom,0.0d0,0.0d0,B0*dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) 
    v_c=(/1.0d0,1.0d0/gm,-mach,0.0d0,0.0d0,B0*dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) 
  !right state
    v_r=(/atwood,1.0d0/gm,-mach,0.0d0,0.0d0,B0*dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) 
  !!!========================================================
  do k=1,kx
     do j=1,jx
        do i=1,ix
        !Define density
        if (x(i) .lt. 0.2d0) then
	    wtr=5.d0*dx(i)
            ro_m(i,j,k)=f_p*(v_l(1)+(v_c(1)-v_l(1))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            p_m(i,j,k)=f_p_p*(v_l(2)+(v_c(2)-v_l(2))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vx_m(i,j,k)=(v_l(3)+(v_c(3)-v_l(3))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vy_m(i,j,k)=(v_l(4)+(v_c(4)-v_l(4))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vz_m(i,j,k)=(v_l(5)+(v_c(5)-v_l(5))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            b_x(i,j,k)=(v_l(6)+(v_c(6)-v_l(6))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            b_y(i,j,k)=(v_l(7)+(v_c(7)-v_l(7))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            b_z(i,j,k)=(v_l(8)+(v_c(8)-v_l(8))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)

            ro_h(i,j,k)=f_n*(v_l(1)+(v_c(1)-v_l(1))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            p_h(i,j,k)=f_p_n*(v_l(2)+(v_c(2)-v_l(2))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vx_h(i,j,k)=(v_l(3)+(v_c(3)-v_l(3))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vy_h(i,j,k)=(v_l(4)+(v_c(4)-v_l(4))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
            vz_h(i,j,k)=(v_l(5)+(v_c(5)-v_l(5))*(1.0d0+dtanh(pi*(x(i))/wtr))*0.5d0)
        else
        if(debug_direction.eq.1) then
        if ((k.eq.1).and.(j.eq.1).and.(i.eq.1)) print*,'linear contact discontinuity'
	wtr=10.0
        ro_m(i,j,k)=f_p*(v_c(1)+(v_r(1)-v_c(1))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        p_m(i,j,k)=f_p_p*(v_c(2)+(v_r(2)-v_c(2))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vx_m(i,j,k)=(v_c(3)+(v_r(3)-v_c(3))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vy_m(i,j,k)=(v_c(4)+(v_r(4)-v_c(4))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vz_m(i,j,k)=(v_c(5)+(v_r(5)-v_c(5))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        b_x(i,j,k)=(v_c(6)+(v_r(6)-v_c(6))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        b_y(i,j,k)=(v_c(7)+(v_r(7)-v_c(7))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        b_z(i,j,k)=(v_c(8)+(v_r(8)-v_c(8))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)

        ro_h(i,j,k)=f_n*(v_c(1)+(v_r(1)-v_c(1))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        p_h(i,j,k)=f_p_n*(v_c(2)+(v_r(2)-v_c(2))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vx_h(i,j,k)=(v_c(3)+(v_r(3)-v_c(3))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vy_h(i,j,k)=(v_c(4)+(v_r(4)-v_c(4))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        vz_h(i,j,k)=(v_c(5)+(v_r(5)-v_c(5))*(1.0d0+dtanh(pi*(x(i)-2.d0-y(j))/10.0/dx(1)))*0.5d0)
        else if(debug_direction.eq.2) then
        if ((k.eq.1).and.(j.eq.1).and.(i.eq.1)) print*,'tanh contact discontinuity'
        ro_m(i,j,k)=f_p*(v_c(1)+(v_r(1)-v_c(1))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                            dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        p_m(i,j,k)=f_p_p*v_c(2)!f_p_p*(v_c(2)+(v_r(2)-v_c(2))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                          !  dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vx_m(i,j,k)=v_c(3)!(v_c(3)+(v_r(3)-v_c(3))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                           ! dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vy_m(i,j,k)=v_c(4)!(v_c(4)+(v_r(4)-v_c(4))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                           ! dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vz_m(i,j,k)=v_c(5)!(v_c(5)+(v_r(5)-v_c(5))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                           ! dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        b_x(i,j,k)=v_c(7)!(v_c(6)+(v_r(6)-v_c(6))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                          !  dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        b_y(i,j,k)=v_c(7)!(v_c(7)+(v_r(7)-v_c(7))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                          !  dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        b_z(i,j,k)=v_c(8)!(v_c(8)+(v_r(8)-v_c(8))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                          !  dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))

        ro_h(i,j,k)=f_n*(v_c(1)+(v_r(1)-v_c(1))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                            dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        p_h(i,j,k)=f_p_n*v_c(2)!f_p_n*(v_c(2)+(v_r(2)-v_c(2))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                            !dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vx_h(i,j,k)=v_c(3)!(v_c(3)+(v_r(3)-v_c(3))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                           ! dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vy_h(i,j,k)=v_c(4)!(v_c(4)+(v_r(4)-v_c(4))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                           ! dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        vz_h(i,j,k)=v_c(5)!(v_c(5)+(v_r(5)-v_c(5))*(1.0d0+dtanh(pi*(x(i)-2.d0-& 
                          !  dtanh(1.0d0*pi*(y(j)-0.5)**2/0.3d0)*0.5d0)/10.0/dx(1))))
        endif
!        elseif (x(i) .lt. 2) then
!            ro_m(i,j,k)=v_c(1)
!            p_m(i,j,k)=v_c(2)
!            vx_m(i,j,k)=v_c(3)
!            vy_m(i,j,k)=v_c(4)
!            vz_m(i,j,k)=v_c(5)
!        elseif (x(i) .ge. 2) then
!            ro_m(i,j,k)=v_r(1)
!            p_m(i,j,k)=v_r(2)
!            vx_m(i,j,k)=v_r(3)
!            vy_m(i,j,k)=v_r(4)
!            vz_m(i,j,k)=v_r(5)
        endif
        enddo
     enddo
  enddo
!print*,flag_mpi,my_rank
!stop
  !!!========================================================
  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=4.0d0
     dtout=tend/4.0/5.d0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine RMI
