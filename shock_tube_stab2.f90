subroutine shock_tube_stab2
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option,nout,time
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use IO_rot,only:restart,output,mk_config,save_coordinates 
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision ::mask(ix,jx,kx)
  double precision ::v_para(ix,jx,kx),v_perp(ix,jx,kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),B0
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr,wtrf
  double precision mach,rcom,rpres,alf,ang, byrat,vyrat,vu,bxu,byu,rou,pru,vxu
  double precision rorand1,rorand2
  integer i,j,k,dpl
integer::read_size=3004,tmpr
double precision :: ro_read(3004),mx_read(3004),my_read(3004),mz_read(3004)
double precision :: en_read(3004),bx_read(3004),by_read(3004),bz_read(3004)
double precision :: ron_read(3004),mnx_read(3004),mny_read(3004),mnz_read(3004)
double precision :: enn_read(3004)
Character(len=40) :: fname,ftime
integer::rds,rde
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
  start(1)=0.0d0 ;end(1)=200.d0!10000.0d0          !400.0d0/f_p
  start(2)=0.0d0 ;end(2)=50.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

!Read in old snapshots
!call restart  
!!	nout=0    
!if(flag_mpi.eq.0 .or. my_rank.eq.0) then
!	call mk_config
!endif 
!time=0.d0
!call save_coordinates 
!call output(0)
!print*,ix

!fname='Data_MHD_1cpu'; ftime='0100'; rds=3549;rde=3650;read_size=3004 !switch-off shock
!fname='Data_MHD_slow'; ftime='0030';rds=4539;rde=4640;read_size=3004 !slow-shock
!fname='Data_MHD_slow2'; ftime='0100';rds=4399;rde=4500;read_size=3004 !slow-shock2
!fname='Data_PIP_slow2'; ftime='0050';rds=1980;rde=2080;read_size=3002 !slow-shock2 PIP

!fname='Data_MHD_par_beta_0.1_m_2'; ftime='0100';rds=1949;rde=2050;read_size=3002 !parallel MHD shock
!fname='Data_PIP_par_beta_0.1_m_2_xn_0.1'; ftime='0090';rds=2639;rde=2740;read_size=3002 !parallel PIP shock
!fname='Data_PIP_par_beta_0.1_m_2_xn_0.9'; ftime='0090';rds=2639;rde=2740;read_size=3002 !parallel PIP shock

fname='MHD_1D_xr_600_xg_3000'; ftime='0030';rds=2580;rde=2679;read_size=3000 !parallel MHD shock

print*,fname

!Gather the data to root
!open the density
    open(101,file=trim(fname)//'/'//trim(ftime)//'ro_p.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)ro_read(1:read_size)
    close(101)
!open the mx
    open(101,file=trim(fname)//'/'//trim(ftime)//'mx_p.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)mx_read(1:read_size)
    close(101)
!open the my
    open(101,file=trim(fname)//'/'//trim(ftime)//'my_p.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)my_read(1:read_size)
    close(101)
!open the mz
    open(101,file=trim(fname)//'/'//trim(ftime)//'mz_p.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)mz_read(1:read_size)
    close(101)
!open the en
    open(101,file=trim(fname)//'/'//trim(ftime)//'en_p.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)en_read(1:read_size)
    close(101)
!open the bx
    open(101,file=trim(fname)//'/'//trim(ftime)//'bx.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)bx_read(1:read_size)
    close(101)
!open the by
    open(101,file=trim(fname)//'/'//trim(ftime)//'by.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)by_read(1:read_size)
    close(101)
!open the bz
    open(101,file=trim(fname)//'/'//trim(ftime)//'bz.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)bz_read(1:read_size)
    close(101)
!print*,ro_read(3550:3650)

if(flag_pip.eq.1) then
    open(101,file=trim(fname)//'/'//trim(ftime)//'ro_n.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)ron_read(1:read_size)
    close(101)
!open the mx
    open(101,file=trim(fname)//'/'//trim(ftime)//'mx_n.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)mnx_read(1:read_size)
    close(101)
!open the my
    open(101,file=trim(fname)//'/'//trim(ftime)//'my_n.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)mny_read(1:read_size)
    close(101)
!open the mz
    open(101,file=trim(fname)//'/'//trim(ftime)//'mz_n.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)mnz_read(1:read_size)
    close(101)
!open the en
    open(101,file=trim(fname)//'/'//trim(ftime)//'en_n.dac.0000',form="unformatted",status="old")
    do i=1,5
       read(101)tmpr
    enddo
    i=tmpr
    read(101)enn_read(1:read_size)
    close(101)
endif

do k=1,kx;do j=1,jx;do i=1,ix
	if (x(1) .le. 0) then
		if(i .lt. 100) then
			U_m(i,j,k,1)=ro_read(rds+i)
			U_m(i,j,k,2)=mx_read(rds+i)
			U_m(i,j,k,3)=my_read(rds+i)
			U_m(i,j,k,4)=mz_read(rds+i)
			U_m(i,j,k,5)=en_read(rds+i)
			U_m(i,j,k,6)=bx_read(rds+i)
			U_m(i,j,k,7)=by_read(rds+i)
			U_m(i,j,k,8)=bz_read(rds+i)
			if(flag_pip.eq.1) then
			U_h(i,j,k,1)=ron_read(rds+i)
			U_h(i,j,k,2)=mnx_read(rds+i)
			U_h(i,j,k,3)=mny_read(rds+i)
			U_h(i,j,k,4)=mnz_read(rds+i)
			U_h(i,j,k,5)=enn_read(rds+i)
			endif
		else
			U_m(i,j,k,1)=ro_read(rde)
			U_m(i,j,k,2)=mx_read(rde)
			U_m(i,j,k,3)=my_read(rde)
			U_m(i,j,k,4)=mz_read(rde)
			U_m(i,j,k,5)=en_read(rde)
			U_m(i,j,k,6)=bx_read(rde)
			U_m(i,j,k,7)=by_read(rde)
			U_m(i,j,k,8)=bz_read(rde)
			if(flag_pip.eq.1) then
			U_h(i,j,k,1)=ron_read(rde)
			U_h(i,j,k,2)=mnx_read(rde)
			U_h(i,j,k,3)=mny_read(rde)
			U_h(i,j,k,4)=mnz_read(rde)
			U_h(i,j,k,5)=enn_read(rde)
			endif
		endif
	else
		U_m(i,j,k,1)=ro_read(rde)
		U_m(i,j,k,2)=mx_read(rde)
		U_m(i,j,k,3)=my_read(rde)
		U_m(i,j,k,4)=mz_read(rde)
		U_m(i,j,k,5)=en_read(rde)
		U_m(i,j,k,6)=bx_read(rde)
		U_m(i,j,k,7)=by_read(rde)
		U_m(i,j,k,8)=bz_read(rde)
		if(flag_pip.eq.1) then
		U_h(i,j,k,1)=ron_read(rde)
		U_h(i,j,k,2)=mnx_read(rde)
		U_h(i,j,k,3)=mny_read(rde)
		U_h(i,j,k,4)=mnz_read(rde)
		U_h(i,j,k,5)=enn_read(rde)
		endif
	endif

call srand(81728)
wtr=1.d0

	if ((x(i) .GE. 50.0d0) .AND. (x(i) .LE. 100.0d0)) then	
!		ro_m(i,j,k)=ro_m(i,j,k)+ro_m(i,j,k)*0.1d0*dsin((x(i)-200.0d0)*3.14d0/100.0d0)* dcos(y(j)*3.14d0/50.0d0*2.d0)
!		ro_m(i,j,k)=ro_m(i,j,k)+0.2d0*dsin((x(i)-4000.d0)/100.0/3.14)*dcos((y(j))/20.d0/3.14d0)
		do dpl=1,10
		rorand1=rand()
		rorand2=rand()
		U_m(i,j,k,1)=U_m(i,j,k,1)+f_p*0.1d0*rorand1*&
 dsin((x(i)-50.0d0)*3.14d0/50.0d0)* dcos((y(j)-rorand2*50.d0)*3.14d0/50.0d0*2.d0*wtr)
		if(flag_pip.eq.1) then
			U_h(i,j,k,1)=U_h(i,j,k,1)+f_n*0.1d0*rorand1*&
 	dsin((x(i)-50.0d0)*3.14d0/50.0d0)* dcos((y(j)-rorand2*50.d0)*3.14d0/50.0d0*2.d0*wtr)
		endif
		wtr=wtr+1.d0
		enddo

	endif
enddo;enddo;enddo

!print*,U_m(3,1,1,:)


!stop
!  do i=1,ix
!	if (i .gt. 7000) then
!		U_m(i,:,:,:)=U_m(7000,:,:,:)
!	else
!		U_m(i,:,:,:)=U_m(i+3500,:,:,:)
!	endif
!  enddo  
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=0.4d0
     dtout=tend/10.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine shock_tube_stab2
