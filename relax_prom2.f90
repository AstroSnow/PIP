subroutine relax_prom2
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,flag_divb &
       ,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,mpi_pos,mpi_siz,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
  use boundary_rot,only:PIPbnd    ! here all the boundary routines are called
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
  double precision pr_val
  real*4 :: lenx,lenz,znull,zaxis
  integer :: num,nx,nz,i,j,k,opt,div_x,div_y,start_x(2),start_y(2),o_height(1)
  character*40 :: filename
  double precision, allocatable:: b_mag(:,:)
  real*4 :: bx(1000,1000),by(1000,1000),bz(1000,1000)
  double precision, allocatable:: p_tot(:),ro_tot(:),temp(:),y_1(:)
  double precision t0,t1,trw,trp,b_norm
  integer, allocatable:: new(:), old(:)
  integer size, seed(2), gseed(2), hiseed(2), zseed(2)
  real harvest(ix*jx*kx)
    data seed /123456789, 987654321/
    data hiseed /-1, -1/
    data zseed /0, 0/
   call random_seed(SIZE=size)

   ALLOCATE (new(size))
   ALLOCATE (old(size))
   CALL RANDOM_SEED(GET=old(1:size))
   new = old*(my_rank+1)
   CALL RANDOM_SEED(PUT=new(1:size))
   call random_number(HARVEST=harvest)


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
  ! load the file
      num=1  ! 1 is the x-point, 3 is the no-x-point
      if(num.eq.0) stop 'Quit'

     nx=500
     nz=500

      write (filename,fmt='("bx_shear",i1,"_lowres.dat")') num
      if (my_rank.eq.1) print *,'Reading data from '//filename
      open(unit=1,file=filename,form='unformatted',status='old')

      read (1) ((bx(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.3) print *, 'loaded b_x'

     close(1)
      if (my_rank.eq.2) print *,'OK'

      write (filename,fmt='("by_shear",i1,"_lowres.dat")') num
      if (my_rank.eq.1) print *,'Reading data from '//filename
      open(unit=1,file=filename,form='unformatted',status='old')

      read (1) ((by(i,j),i=1,nx),j=1,nz)
      if (my_rank.eq.3) print *, 'loaded b_y'

     close(1)

      write (filename,fmt='("bz_shear",i1,"_lowres.dat")') num
      if (my_rank.eq.1) print *,'Reading data from '//filename
      open(unit=1,file=filename,form='unformatted',status='old')

      read (1) ((bz(i,j),i=1,nx),j=1,nz)
!      if (my_rank.eq.3) print *, 'loaded bz',(bz(251,i),i=1,nz)

     close(1)
      if (my_rank.eq.2) print *, 'OK'
     ! NEED TO SET UP B FOR DIFFERENT CORES!!

    div_x=ix-2*margin(1)
    div_y=jx-2*margin(2)
    if (my_rank.eq.3) print *, div_x,div_y, my_rank


    start_x(1)=1+mpi_pos(1)*div_x
    start_x(2)=(mpi_pos(1)+1)*div_x+2*margin(1)

    start_y(1)=1+mpi_pos(2)*div_y
    start_y(2)=(mpi_pos(2)+1)*div_y+2*margin(2)

!   if (my_rank.eq.0) print *, ix,jx,kx
!  if (my_rank.eq.1) 
!   print *, 'start',start_x,start_y,my_rank
    ! deal with Margin!!

   !-----------------------------------------------------
  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-2.5d0 ;end(1)=2.5d0
  start(2)=0.0d0 ;end(2)=5.0d0
  start(3)=0.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=2
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=2
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=4
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=2
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !density of lower fluid is unity
  nz=500
  nx=500

  flag_grav=1
  gra(:,:,:,:)=0.0d0
  scl_height=1.d0
  gra(:,:,:,2)=-1.0d0/(gm*scl_height)

  if (mpi_pos(2).eq.0) then
  gra(:,1:margin(2),:,2)=1.0d0/(gm*scl_height)
  endif
  if (mpi_pos(2).eq.(mpi_siz(2)-1)) then
  gra(:,jx-margin(2)+1:jx,:,2)=1.0d0/(gm*scl_height)
  endif

  allocate(p_tot(nz+2*margin(2)),ro_tot(nz+2*margin(2)) &
     ,temp(nz),y_1(nz))


!      write (filename,fmt='("temp",i1,"_lowres.dat")')
!      if (my_rank.eq.1) print *,'Reading data from '//filename
!      open(unit=1,file=filename,form='unformatted',status='old')

!      read (1) (temp(i),i=1,nz)

!      close(1)

  y_1(1)=dy(1)/2.d0
  do i=2,nz
  y_1(i)=y_1(i-1)+dy(1)
  enddo
!  if (my_rank.eq.2) print *, 'dz=',dy(1)
!  if (my_rank.eq.1) print *, 'y =',y_1
  temp(:)=1.d0/25.d0+0.5d0*(1.d0-1.d0/25.d0)*(tanh((y_1-1.d0/2.5d0)/0.07d0) +1.d0)
!  if (my_rank.eq.0)  print *, 'temp=',temp
  pr_val=1.d0/gm

  do  i=margin(2)+1,nz+margin(2)
    if (i .eq. margin(2)+1) then
      p_tot(i)=pr_val
    else
      p_tot(i)=p_tot(i-1) - dy(1) * p_tot(i-1)/temp(i-1-margin(2))
    endif  
  enddo
    p_tot(1:margin(2))=p_tot(margin(2)*2:margin(2)+1:-1)
    p_tot((nz+margin(2)+1):(nz+2*margin(2)))=p_tot(nz+margin(2):nz+1:-1)
!  if (my_rank.eq.0)  print *, 'p=',p_tot
!    p_tot=p_tot/p_tot(150)/gm
!  if (my_rank.eq.0)  print *, 'p=',p_tot 
    ro_tot(margin(2)+1:nz+margin(2))=p_tot(margin(2)+1:nz+margin(2)) &
        *gm/temp(:)
    ro_tot(1:margin(2))=ro_tot(margin(2)*2:margin(2)+1:-1)
    ro_tot((nz+margin(2)+1):(nz+2*margin(2)))=ro_tot(nz+margin(2):nz+1:-1)


!---------------------



!    if (my_rank.eq.0) print *, ro_tot
!    print *, 2*margin(2)+start_y   
!    if (my_rank.eq.0) print *, p_tot
  ! MAKW MPI READY
!    P_m= pr_val*f_p_p 
    P_m=spread(spread(p_tot(start_y(1):start_y(2)) &
      *f_p_p,1,ix),3,kx)
!    ro_m= 1.d0*f_p 
     ro_m= spread(spread(ro_tot(start_y(1):start_y(2)) & 
      *f_p ,1,ix),3,kx)         
!    P_h= pr_val*f_p_n 
     P_h= spread(spread(p_tot(start_y(1):start_y(2)) & 
      *f_p_n ,1,ix),3,kx)
!    ro_h= 1.d0*f_n 
     ro_h= spread(spread(ro_tot(start_y(1):start_y(2)) &
      *f_n ,1,ix),3,kx)


   o_height(1)=maxloc(abs(bz(250,60:300)),1)+59
!   print *, 'o_height',o_height, bz(250,o_height(1)),my_rank
   b_norm=sqrt(p_tot(o_height(1)-margin(2))*2.d0/beta)/abs(bz(250,o_height(1)))

!sqrt(p_tot(164)/2.d0/beta)/bz(499,164)
!   if (my_rank.eq.0) print *, b_norm,bz(499,328),maxloc(abs(bz(499,:)))

!  STOP
  vx_h=0.0d0;vy_h=0.0d0;vz_h=0.0d0
  vx_m=0.0d0;vy_m=0.0d0;vz_m=0.0d0
!  B_x=0.0d0;B_y=0.0d0;B_z=0.0d0

   allocate(b_mag(nx+2*margin(1),nz+2*margin(2)))
   b_mag(margin(1)+1:nx+margin(1),margin(2)+1:nz+margin(2))= &
    bx

   b_mag(1:margin(1),:)=-b_mag(2*margin(1):margin(1)+1:-1,:)
   b_mag(nx+margin(1)+1:nx+2*margin(1),:)= &
     -b_mag(nx+1:nx+margin(1):-1,:)

   b_mag(:,1:margin(2))=-b_mag(:,2*margin(2):margin(2)+1:-1)
   b_mag(:,nz+margin(2)+1:nz+2*margin(2))= &
     b_mag(:,nz+margin(2):nz+1:-1)

    B_x=DBLE(spread(b_mag(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))*b_norm
    if (my_rank.eq.3) print *, 'done bx', my_rank

   b_mag(margin(1)+1:nx+margin(1),margin(2)+1:nz+margin(2))= &
    by

    b_mag(1:margin(1),:)=b_mag(2*margin(1):margin(1)+1:-1,:)
   b_mag(nx+margin(1)+1:nx+2*margin(1),:)= &
     b_mag(nx+1:nx+margin(1):-1,:)

   b_mag(:,1:margin(2))=b_mag(:,2*margin(2):margin(2)+1:-1)
   b_mag(:,nz+margin(2)+1:nz+2*margin(2))= &
     -b_mag(:,nz+margin(2):nz+1:-1)


    B_y=DBLE(spread(b_mag(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))*b_norm
!    print *, b_norm,my_rank
    if (my_rank.eq.2) print *, 'done by',my_rank

   b_mag(margin(1)+1:nx+margin(1),margin(2)+1:nz+margin(2))= &
    bz

   b_mag(1:margin(1),:)=-b_mag(2*margin(1):margin(1)+1:-1,:)
   b_mag(nx+margin(1)+1:nx+2*margin(1),:)= &
     -b_mag(nx+1:nx+margin(1):-1,:)

   b_mag(:,1:margin(2))=-b_mag(:,2*margin(2):margin(2)+1:-1)
   b_mag(:,nz+margin(2)+1:nz+2*margin(2))= &
     b_mag(:,nz+margin(2):nz+1:-1)


    B_z=DBLE(spread(b_mag(start_x(1):start_x(2),start_y(1) &
      :start_y(2)),3,kx))*b_norm
    if (my_rank.eq.1) print *, 'done bz',my_rank

!    print *, B_z(200,1:12,1)
!    STOP

!   B_x=0.0d0;B_y=0.0d0;B_z=0.0d0  

  !!!=================
  ! Add mass
  
  ro_m=ro_m+15.d0*f_p*spread(0.5d0*(1.d0-dtanh((abs(spread(x,2,jx))-0.1d0 )/0.1d0)) &
     *0.5d0*(1.d0-dtanh((abs(spread(y-1.1d0,1,ix))-0.3d0 )/0.05d0)),3,kx)* &
!     spread(spread(1.d0+cos(z*2.d0*pi/8.d0),1,ix),2,jx)* &
     ro_tot(o_height(1)-margin(2))
  ro_h=ro_h+15.d0*f_n*spread(0.5d0*(1.d0-dtanh((abs(spread(x,2,jx))-0.1d0 )/0.1d0)) &
     *0.5d0*(1.d0-dtanh((abs(spread(y-1.1d0,1,ix))-0.3d0 )/0.05d0)),3,kx)* &
!     spread(spread(1.d0+cos(z*2.d0*pi/8.d0),1,ix),2,jx)* &
     ro_tot(o_height(1)-margin(2))
!   do k=1,kx
!     do j=1,jx
!        do i=1,ix
!        if (abs(x(i)).lt.0.15d0 .and. abs(y(j)-1.1d0).lt.0.4d0) then
!        ro_m(i,j,k)=ro_m(i,j,k)*(1.d0+0.2d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0))
!        ro_h(i,j,k)=ro_h(i,j,k)*(1.d0+0.2d0*(harvest((k-1)*jx*ix+(j-1)*ix+i)-0.5d0))
!        endif
!        enddo
!     enddo
!  enddo

  deallocate(ro_tot,temp,y_1)

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=40.0d0
     dtout=tend/40.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine relax_prom2
