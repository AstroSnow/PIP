module model_rot
!====================================================================
! This is the model module for the PIP code.
!Module for initial settings of the simulation (distributions of physical parameters are set in the initial module)
! Author A.Hillier
! First version - 2013/04/23
! Second version - 2013/04/24 NN
! Third version - 2013/04/26 NN
! Fourth version - 2013/04/29 AH - 
! Fifth version - 2013/06/4 NN
!====================================================================
  use globalvar,only:ndim,ix,jx,kx,flag_mhd,flag_pip,ig,margin,&
       x,y,z,dx,dy,dz,dxc,dyc,dzc,flag_resi,flag_ir,flag_grav,&
       col,db_clean,debug_direction,debug_option,debug_parameter,&
       flag_sch,flag_artvis,flag_amb,eta_0,dtout,beta,my_rank,flag_mpi,&
       flag_divb,flag_heat,flag_ini,flag_mpi_split,flag_pip_imp,&
       U_h,U_m,s_order,flag_col,nvar_h,nvar_m,safety,t_ir,sor_omega,&
       pip_imp_factor,t_order,tend,theta,mpi_siz,outdir,kapper,hc_split,indir,&
       n_fraction,hc_sch,hc_max,gm,flag_restart,flag_debug,mpi_siz,mpi_pos,&
       flag_bnd,dsc,output_type,flag_hll,flag_time,hc_integ,hc_type,ro_lim,pr_lim, &
       flag_hc_test,safety_cnd,nsub_max,b_cr,flag_ps,flag_cyl, &
       vd_cri,j_cri,flag_damp,damp_time,flag_rad, T0, n0, L0,flag_IR_type, &
       flag_visc, nu_0, esav, emsavtime, &
	ac_sav, xi_sav, ion_sav, rec_sav, col_sav, gr_sav, vs_sav, heat_sav, et_sav, ps_sav,&
	rad_ts,radrhoref
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use HC_rot,only:initialize_HC
  use Res_rot,only:initialize_resistivity
  use Gra_rot,only:initialize_gravity
  use visc_rot,only:initialize_visc
  use PIP_rot,only:initialize_collisional,initialize_IR,initialize_xin
  use Util_rot,only:get_word
  use Rad_cooling,only:initialize_radloss
  implicit none
contains 
!------------------------------------------------------------------------



subroutine get_parameters
! this reads the parameters from params.txt
  character*15 key
  character*8 flag_eqs
  character*200  tmp
  integer ind_e,ix0,jx0,kx0,i_tmp,r_count
  double precision gm_d,gm_n
  open(11,file='setting.txt',form='formatted',status='old')
  r_count=0
  do 
     read(11,"(A)",end=999)tmp
     call get_word(tmp,key,ind_e)         
     if(key.eq.'flag_eqs') then
        flag_eqs=tmp(1:ind_e-1)        
     else if(key.eq.'outdir') then
        outdir=tmp(1:ind_e-1)
     else if(key.eq.'indir') then
        indir=tmp(1:ind_e-1)        
     else if(key.eq.'ndim') then
        read(tmp(1:ind_e-1),*)ndim
     else if(key.eq.'ix') then
        read(tmp(1:ind_e-1),*)ix0
     else if(key.eq.'jx') then
        read(tmp(1:ind_e-1),*)jx0
     else if(key.eq.'kx') then
        read(tmp(1:ind_e-1),*)kx0
     else if(key.eq.'gm_d') then
        read(tmp(1:ind_e-1),*)gm_d
     else if(key.eq.'gm_n') then
        read(tmp(1:ind_e-1),*)gm_n
     else if(key.eq.'safety') then
        read(tmp(1:ind_e-1),*)safety
     else if(key.eq.'safety_cnd') then
        read(tmp(1:ind_e-1),*)safety_cnd
     else if(key.eq.'nsub_max') then
        read(tmp(1:ind_e-1),*)nsub_max
     else if(key.eq.'flag_ini') then
        flag_ini=tmp(1:ind_e-1)
!        read(tmp(1:ind_e-1),*)flag_ini
     else if(key.eq.'flag_sch') then
        read(tmp(1:ind_e-1),*)flag_sch
     else if(key.eq.'flag_cyl') then
        read(tmp(1:ind_e-1),*)flag_cyl
     else if(key.eq.'s_order') then
        read(tmp(1:ind_e-1),*)s_order
     else if(key.eq.'t_order') then
        read(tmp(1:ind_e-1),*)t_order
     else if(key.eq.'flag_grav') then
        read(tmp(1:ind_e-1),*)flag_grav
     else if(key.eq.'flag_visc') then
        read(tmp(1:ind_e-1),*)flag_visc
     else if(key.eq.'nu_0') then
        read(tmp(1:ind_e-1),*)nu_0
     else if(key.eq.'flag_bnd_x_l') then
        read(tmp(1:ind_e-1),*)flag_bnd(1)
     else if(key.eq.'flag_bnd_x_r') then
        read(tmp(1:ind_e-1),*)flag_bnd(2)
     else if(key.eq.'flag_bnd_y_l') then
        read(tmp(1:ind_e-1),*)flag_bnd(3)
     else if(key.eq.'flag_bnd_y_r') then
        read(tmp(1:ind_e-1),*)flag_bnd(4)
     else if(key.eq.'flag_bnd_z_l') then
        read(tmp(1:ind_e-1),*)flag_bnd(5)
     else if(key.eq.'flag_bnd_z_r') then
        read(tmp(1:ind_e-1),*)flag_bnd(6)
     else if(key.eq.'flag_restart') then
        read(tmp(1:ind_e-1),*)flag_restart
     else if(key.eq.'flag_artvis') then
        read(tmp(1:ind_e-1),*)flag_artvis
     else if(key.eq.'dtout')then
        read(tmp(1:ind_e-1),*)dtout
     else if(key.eq.'tend') then
        read(tmp(1:ind_e-1),*)tend        
     else if(key.eq.'output_type') then
        read(tmp(1:ind_e-1),*)output_type       
     else if(key.eq.'beta') then
        read(tmp(1:ind_e-1),*)beta
!	if(beta .LT. 0.0) then
!		read(*,*) beta
!	endif
     else if(key.eq.'flag_resi') then
        read(tmp(1:ind_e-1),*)flag_resi
     else if(key.eq.'eta_0') then
        read(tmp(1:ind_e-1),*)eta_0
     else if(key.eq.'vd_cri') then
        read(tmp(1:ind_e-1),*)vd_cri
     else if(key.eq.'j_cri') then
        read(tmp(1:ind_e-1),*)j_cri
     else if(key.eq.'flag_amb') then
        read(tmp(1:ind_e-1),*)flag_amb
     else if(key.eq.'flag_divb') then
        read(tmp(1:ind_e-1),*)flag_divb
     else if(key.eq.'n_fraction') then
        read(tmp(1:ind_e-1),*)n_fraction
!	if(n_fraction .LT. 0.0) then
!		read(*,*) n_fraction
!	endif
     else if(key.eq.'col') then
        read(tmp(1:ind_e-1),*)col
     else if(key.eq.'flag_pip_imp') then
        read(tmp(1:ind_e-1),*)flag_pip_imp
     else if(key.eq.'flag_col') then
        read(tmp(1:ind_e-1),*)flag_col
     else if(key.eq.'flag_IR') then
        read(tmp(1:ind_e-1),*)flag_IR
     else if(key.eq.'flag_IR_type') then
        read(tmp(1:ind_e-1),*)flag_IR_type
     else if(key.eq.'t_IR') then
        read(tmp(1:ind_e-1),*)t_ir
     else if(key.eq.'pip_imp_factor') then
        read(tmp(1:ind_e-1),*)pip_imp_factor
     else if(key.eq.'db_clean') then
        read(tmp(1:ind_e-1),*)db_clean
     else if(key.eq.'b_cr') then
        read(tmp(1:ind_e-1),*)b_cr
     else if(key.eq.'theta') then
        read(tmp(1:ind_e-1),*)theta
     else if(key.eq.'flag_heat') then
        read(tmp(1:ind_e-1),*)flag_heat
     else if(key.eq.'hc_split') then
        read(tmp(1:ind_e-1),*)hc_split
     else if(key.eq.'hc_sch') then
        read(tmp(1:ind_e-1),*)hc_sch
     else if(key.eq.'kapper') then
        read(tmp(1:ind_e-1),*)kapper
     else if(key.eq.'sor_omega') then
        read(tmp(1:ind_e-1),*)sor_omega
     else if(key.eq.'hc_max') then
        read(tmp(1:ind_e-1),*)hc_max
     else if(key.eq.'flag_mpi') then
        read(tmp(1:ind_e-1),*)flag_mpi
     else if(key.eq.'flag_mpi_split') then
        read(tmp(1:ind_e-1),*)flag_mpi_split
     else if(key.eq.'mpi_x') then
        read(tmp(1:ind_e-1),*)mpi_siz(1)
     else if(key.eq.'mpi_y') then
        read(tmp(1:ind_e-1),*)mpi_siz(2)
     else if(key.eq.'mpi_z') then
        read(tmp(1:ind_e-1),*)mpi_siz(3)
     else if(key.eq.'flag_debug') then
        read(tmp(1:ind_e-1),*)flag_debug
     else if(key.eq.'debug_option') then
        read(tmp(1:ind_e-1),*)debug_option
     else if(key.eq.'debug_parameter') then
        read(tmp(1:ind_e-1),*)debug_parameter
     else if(key.eq.'debug_direction') then
        read(tmp(1:ind_e-1),*)debug_direction
     else if(key.eq.'flag_hll') then
        read(tmp(1:ind_e-1),*)flag_hll
     else if(key.eq.'flag_time') then
        read(tmp(1:ind_e-1),*)flag_time
     else if(key.eq.'hc_integ') then
        read(tmp(1:ind_e-1),*)hc_integ
     else if(key.eq.'hc_type') then
        read(tmp(1:ind_e-1),*)hc_type
     else if(key.eq.'ro_lim') then
        read(tmp(1:ind_e-1),*)ro_lim
     else if(key.eq.'pr_lim') then
        read(tmp(1:ind_e-1),*)pr_lim        
     else if(key.eq.'flag_ps') then
        read(tmp(1:ind_e-1),*)flag_ps
    else if(key.eq.'flag_damp') then
        read(tmp(1:ind_e-1),*)flag_damp
     else if(key.eq.'damp_time') then
        read(tmp(1:ind_e-1),*)damp_time
     else if(key.eq.'flag_rad') then
        read(tmp(1:ind_e-1),*)flag_rad
     else if(key.eq.'rad_ts') then
        read(tmp(1:ind_e-1),*)rad_ts
     else if(key.eq.'radrhoref') then
        read(tmp(1:ind_e-1),*)radrhoref
     else if(key.eq.'T_norm') then
        read(tmp(1:ind_e-1),*)T0
     else if(key.eq.'n_norm') then
        read(tmp(1:ind_e-1),*)n0
     else if(key.eq.'L_norm') then
        read(tmp(1:ind_e-1),*)L0
     else if(key.eq.'esav') then
        read(tmp(1:ind_e-1),*)esav
     else if(key.eq.'emsavtime') then
        read(tmp(1:ind_e-1),*)emsavtime
     else if(key.eq.'ac_sav') then
        read(tmp(1:ind_e-1),*)ac_sav
     else if(key.eq.'xi_sav') then
        read(tmp(1:ind_e-1),*)xi_sav
     else if(key.eq.'ion_sav') then
        read(tmp(1:ind_e-1),*)ion_sav
     else if(key.eq.'rec_sav') then
        read(tmp(1:ind_e-1),*)rec_sav
     else if(key.eq.'col_sav') then
        read(tmp(1:ind_e-1),*)col_sav
     else if(key.eq.'gr_sav') then
        read(tmp(1:ind_e-1),*)gr_sav
     else if(key.eq.'vs_sav') then
        read(tmp(1:ind_e-1),*)vs_sav
     else if(key.eq.'heat_sav') then
        read(tmp(1:ind_e-1),*)heat_sav
     else if(key.eq.'et_sav') then
        read(tmp(1:ind_e-1),*)et_sav
     else if(key.eq.'ps_sav') then
        read(tmp(1:ind_e-1),*)ps_sav
     endif    
     !Make config is delegated to mod.IO_rot sub.mk_config
!     if(flag_mpi.eq.0 .or. my_rank.eq.0) then
!        if(r_count.eq.0) then
!           open(99,file=trim(outdir) // "/config.txt",status="replace",form="formatted")
!        endif
!        if(ind_e.gt.1) write(99,"(A)")key//":"//tmp(1:ind_e-1)
!     endif 
!     r_count=r_count+1
   enddo
 999 continue
!   if(flag_mpi.eq.0 .or. my_rank.eq.0)write(99,"(A)")"ENDSETTING"
   close(11)
  !Set Basic equation
   select case(trim(adjustl(flag_eqs)))
   case('HD')
      flag_mhd=0
      flag_pip=0
      nvar_h=5
      nvar_m=0
      flag_col=0
      flag_IR=0
      flag_amb=0
      flag_resi=0
   case('MHD')
      flag_mhd=1
      flag_pip=0
      nvar_h=0
      nvar_m=8
!      if (flag_amb.eq.0) flag_col=0 !! by Tak
      flag_IR=0
   case('PIP')
      flag_mhd=1
      flag_pip=1
      nvar_h=5
      nvar_m=8
      flag_amb=0
   end select

   !Set Scheme depend parameter for each scheme
   select case(flag_sch)
   case(0) ! HLL
!      t_order=s_order
      if(flag_time.eq.0) t_order=s_order ! by Tak
      flag_artvis=0 ! by Tak
      margin(1)=s_order+min(flag_resi,1)+min(flag_amb,1)+flag_artvis*2
      if(flag_heat.ge.1) then
         margin(1)=margin(1)+1
      endif
   case(1) ! SLW
!      margin(1)=4+mod(max(flag_resi,flag_amb),2)*2
      s_order=4
!      t_order=s_order
      if(flag_time.eq.0) t_order=s_order ! by Tak
      if(flag_artvis.eq.2) then
         margin(1)=margin(1)+2         
         i_tmp=mod(flag_resi,10)
         flag_resi=flag_resi-i_tmp+2
      endif
      if(flag_artvis.eq.-1)flag_artvis=1
      margin(1)=s_order+min(flag_resi,1)+min(flag_amb,1)+flag_artvis*2
      if(flag_divb.gt.2) flag_divb=1	!ORIGINALLY 1, SHOULD BE 2 for iterative
   case(2)
      margin(1)=4+flag_amb*2
      if(flag_divb.gt.2) flag_divb=1	!ORIGINALLY 1, SHOULD BE 2 for iterative
      t_order=4
   end select
   !modification for divb cleaning
   if(ndim.eq.1) flag_divb=0
   select case(flag_divb)
   case(1)
      nvar_m=nvar_m+1
   case(2)
      nvar_m=nvar_m+1
   end select

   !set dimension
   select case(ndim)
   case(1)
      margin(2:3)=0
      jx0=1
      kx0=1
   case(2)
      margin(2)=margin(1)
      margin(3)=0
      kx0=1
   case(3)
      margin(2:3)=margin(1)
   end select
   ix=ix0+2*margin(1)
   jx=jx0+2*margin(2)
   kx=kx0+2*margin(3)
   gm=gm_d/dble(gm_n)


   !! To avoid errors...
   if(flag_heat.eq.1 .and. flag_eqs.eq.'PIP') then
      print*,'heat conduction for PIP is still under construction...'
      print*,'STOP'
      stop
   endif
   if(flag_heat.eq.1 .and. hc_integ.eq.1 .and. hc_type.eq.0) then
      print*,'constant heat conductivity for implicit method has not been implemented yet'
      print*,'STOP'
      stop
   endif
   ! if(flag_heat.eq.1 .and. hc_integ.eq.1 .and. flag_mpi.eq.1) then
   !    print*,'implicit method has not been parallerized yet.'
   !    print*,'STOP'
   !    stop
   ! endif
   if(flag_heat.eq.1 .and. hc_integ.eq.0) then
      print*,'Explicit integration method for heat conduction has not been implemented yet'
      print*,'STOP'
      stop
   endif
   if(flag_cyl.eq.1 .and. ndim.ne.2) then
      print*,'Cylindrical coord is available only when ndim=2'
      print*,'STOP'
   endif

   
   !! Default setting is flag_hc_test = 0
   flag_hc_test = 0
   
!   ! copy the setting file (setting.txt) to the input and output directories
   if(flag_mpi.eq.0 .or. my_rank.eq.0) then
!      call system('cp setting.txt '//trim(adjustl(indir)))
      call system('cp setting.txt '//trim(adjustl(outdir)))
   endif
   
 end subroutine get_parameters


 !Allocate the array sizes
 subroutine allocate_vars
   ig(1)=ix ; ig(2)=jx ; ig(3)=kx   
   if(flag_mpi.eq.0 .or.my_rank.eq.0) print *,"total grid :",ix,jx,kx
   allocate(x(ix),dx(ix),y(jx),dy(jx),z(kx),dz(kx),dsc(max(ix,jx,kx),3),&
        dxc(ix),dyc(jx),dzc(kx))
   if(ndim.le.2) dz(1)=1.0d10
   if(ndim.le.1) dy(1)=1.0d10
   
   allocate(u_h(ix,jx,kx,nvar_h),u_m(ix,jx,kx,nvar_m))   
   call initialize_xin(flag_amb,flag_col)
   call initialize_collisional(flag_col)
   call initialize_IR(flag_IR)
   call initialize_resistivity(flag_resi)
   call initialize_HC(flag_heat)   
   call initialize_gravity(flag_grav)   
   call initialize_visc(flag_visc)   
   call initialize_radloss(flag_rad)      
 end subroutine allocate_vars
 
 subroutine set_coordinate(start,end)
   double precision,intent(inout)::start(3),end(3)
   integer i,j,k,n,zeros(3)
   double precision  ds,length
   double precision dr(3)

   !set region for MPI domain decomposition
   if(ndim.le.2) then
      start(3)=0.0d0
      end(3)=0.0d0
      if(ndim.le.1) then
         start(2)=0.0d0
         end(2)=0.0d0
      endif
   endif

   if(flag_mpi.eq.1) then 
      do n=1,3
         length=end(n)-start(n)
         ds=1.0d0/mpi_siz(n)
         start(n)=start(n)+ds*mpi_pos(n)*length
         end(n)=start(n)+ds*length
      enddo
   endif

   !!set grid width   
   dr(1)=(end(1)-start(1))/(ix-2*margin(1)) 
   dr(2)=(end(2)-start(2))/(jx-2*margin(2))
   dr(3)=(end(3)-start(3))/(kx-2*margin(3))
   dx(:)=dr(1);dy(:)=dr(2);dz(:)=dr(3)
   if (ndim.le.2) then
      dz(1)=100
      if(ndim.le.1) then
         dy(1)=100
      endif
   endif
   dxc=dx
   dyc=dy
   dzc=dz   
   !!set origin
   zeros(:)=margin(:)+1
   x(zeros(1))=start(1)+0.5d0*dr(1)
   y(zeros(2))=start(2)+0.5d0*dr(2)
   z(zeros(3))=start(3)+0.5d0*dr(3)

   do i=zeros(1)+1,ix
      x(i)=x(i-1)+dx(i-1)
   enddo
   do i=zeros(1)-1,1,-1
      x(i)=x(i+1)-dx(i)
   enddo   
   do j=zeros(2)+1,jx
      y(j)=y(j-1)+dy(j-1)
   enddo
   do j=zeros(2)-1,1,-1
      y(j)=y(j+1)-dy(j)
   enddo
   
   do k=zeros(3)+1,kx
      z(k)=z(k-1)+dz(k-1)
   enddo
   do k=zeros(3)-1,1,-1
      z(k)=z(k+1)-dz(k)
   enddo
 end subroutine set_coordinate

 subroutine set_dsc
   integer i,j,k
   if (s_order.eq.4) then
      do i=margin(1)+1,ix-margin(1)
         dxc(i)=0.25d0*(dx(i+1)+dx(i)+dx(i-1)+dx(i-2))
      enddo
      do i=1,margin(1)
         dxc(i)=dxc(2*margin(1)+1-i)
         dxc(ix+1-i)=dxc(ix-2*margin(1)+i)
      enddo
      if(ndim.ge.2) then
         do j=margin(2)+1,jx-margin(2)
            dyc(j)=0.25d0*(dy(j+1)+dy(j)+dy(j-1)+dy(j-2))
         enddo
         do j=1,margin(2)
            dyc(j)=dyc(2*margin(2)+1-j)
            dyc(jx+1-j)=dxc(jx-2*margin(2)+j)
         enddo         
      else
         dyc=dy(1)
         dzc=dz(1)
      endif 
   else
   endif
 end subroutine set_dsc

 subroutine setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
          ro_h,vx_h,vy_h,vz_h,p_h)
   double precision,intent(inout) :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
   double precision,intent(inout) :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
   double precision,intent(inout) :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
   double precision,intent(inout) :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
   double precision,intent(inout) :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
   double precision ,intent(inout):: B_x (1:ix,1:jx,1:kx)
   double precision,intent(inout) :: B_y (1:ix,1:jx,1:kx)
   double precision,intent(inout) :: B_z (1:ix,1:jx,1:kx)
   
   if(flag_mhd.eq.1) then      
      call pv2cq_mhd(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z,U_m)
      if(flag_pip.eq.1) then
         call pv2cq_hd(ro_h,vx_h,vy_h,vz_h,p_h,U_h)
      endif
   else
      call pv2cq_hd(ro_h,vx_h,vy_h,vz_h,p_h,U_h)
   endif
 end subroutine setcq
end module model_rot
