module HC_rot
  use globalvar,only:ix,jx,kx,kapper,dx,dy,dz,dxc,dyc,dzc,gm,margin,ndim,&
       flag_mhd,dt,dt_cnd,nt,hc_split,hc_sch,hc_max,x,y,z,sor_omega,total_iter,&
       neighbor,mpi_siz,total_prc,my_rank,x,y,z,tiny,hc_type,flag_mpi,flag_bnd,&
       safety_cnd,nsub_max,time,pr_lim,ro_lim,flag_cyl
  use scheme_rot,only:cq2pv_hd,pv2cq_hd,cq2pv_mhd,pv2cq_mhd
  use MPI_rot,only:mpi_double_interface
  use boundary_rot,only:PIPbnd,bnd_energy
  implicit none
  include "mpif.h"
  integer,save::n_hc
  integer,save::is,ie,js,je,ks,ke
!  integer,parameter::hc_split=1,hc_sch=2,hc_max=500
  integer,parameter::hc_unit=82
  integer,allocatable,save::di(:),dj(:),dk(:)
  integer ierr,status(6)
  double precision,parameter::error=1.0d-5,epsi=1.0d-5
contains
  subroutine initialize_HC(flag_hc)
    integer,intent(inout)::flag_hc
    if(flag_hc.eq.0) return

    total_iter=0
    is=margin(1)+1;ie=ix-margin(1)
    js=margin(2)+1;je=jx-margin(2)        
    ks=margin(3)+1;ke=kx-margin(3)        
    select case(ndim)
    case(1)
       n_hc=3
       allocate(di(n_hc),dj(n_hc),dk(n_hc))
       di=(/0,-1,1/)
       dj=(/0,0,0/)
       dk=(/0,0,0/)
    case(2)
       select case(hc_split)
       case(1)
          n_hc=5
          allocate(di(n_hc),dj(n_hc),dk(n_hc))
          di=(/0,-1,1,0,0/)
          dj=(/0,0,0,-1,1/)
          dk=(/0,0,0,0,0/)
       case(2)
          n_hc=9
          allocate(di(n_hc),dj(n_hc),dk(n_hc))          
          di=(/0,-1, 1, 0,0,-1, 1,-1,1/)
          dj=(/0, 0, 0,-1,1,-1,-1, 1,1/)
          dk=(/0, 0, 0, 0,0, 0, 0, 0,0/)
       end select
    case(3)
       select case(hc_split)
       case(1)
          n_hc=7
          allocate(di(n_hc),dj(n_hc),dk(n_hc))
          di=(/0,-1,1,0,0,0,0/)
          dj=(/0,0,0,-1,1,0,0/)
          dk=(/0,0,0,0,0,-1,1/)
       case(2)
          n_hc=19
          allocate(di(n_hc),dj(n_hc),dk(n_hc))           
          di=(/0,-1,1,0,0,0,0,-1, 1,-1,1,-1, 1,-1,1, 0, 0, 0,0/)
          dj=(/0,0,0,-1,1,0,0,-1,-1, 1,1, 0, 0, 0,0,-1, 1,-1,1/)
          dk=(/0,0,0,0,0,-1,1, 0, 0, 0,0,-1,-1, 1,1,-1,-1, 1,1/)
       end select
    end select
  end subroutine initialize_HC

  subroutine HC(U,nvar,mhd)
    integer,intent(in)::nvar,mhd
    double precision,intent(inout)::U(ix,jx,kx,nvar)
    double precision ro(ix,jx,kx),pr(ix,jx,kx),te(ix,jx,kx)
    double precision vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    double precision bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision matrix(ix,jx,kx,n_hc),source(ix,jx,kx)
    integer count


    if(mhd.eq.1) then
       call cq2pv_mhd(ro,vx,vy,vz,pr,bx,by,bz,U)
    else
       call cq2pv_hd(ro,vx,vy,vz,pr,U)
    endif

    te=gm*(pr/ro)    

    !set matrix and source term  
    
    if(total_prc.ge.2) then       
       call mpiBND(ro,2)
       call mpiBND(te,2)
       call mpiBND(bx,2)
       call mpiBND(by,2)
       call mpiBND(bz,2)
    endif

    if(flag_mhd.eq.1.and.ndim.ge.2) then
       select case(hc_split)
       case(1)
          if(ndim.eq.1) then
!             call set_matrix_source_MHD2D5(matrix,source,ro,te,bx,by,bz)
             call set_matrix_source_1D(matrix,source,ro,te)
          else if(ndim.eq.2) then
             call set_matrix_source_MHD2D5(matrix,source,ro,te,bx,by,bz)
          else 
!             call set_matrix_source_MHD3D7(matrix,source,ro,te,bx,by,bz)
          endif
       case(2)
          if(ndim.eq.2) then             
             call set_matrix_source_MHD2D9(matrix,source,ro,te,bx,by,bz)
!             call set_matrix_source_MHD2D9_limiter(matrix,source,ro,te,bx,by,bz)
          else 
!             call set_matrix_source_MHD3D19(matrix,source,ro,te,bx,by,bz)
          endif
       end select
    else
       if(ndim.eq.1) then
          call set_matrix_source_1D(matrix,source,ro,te)       
       endif
    endif
!    call test_matrix(matrix,source,te) 
    !set boundary condition of matrix
    call set_matrix_boundary(matrix) 
 

 
    !iteration schemes    
    select case(hc_sch)
    case(1)
       call SOR(matrix,te,source,count)
    case(2)
       if(total_prc.le.1) then
          call BiCGstab(matrix,te,source,count)
       else
          call BiCGstabMPI(matrix,te,source,count)
       endif
    end select    
    total_iter=total_iter+count
    !back to conservative
    pr=(te*ro)/gm
    call pv2cq_mhd(ro,vx,vy,vz,pr,bx,by,bz,U)        
  end subroutine HC

  subroutine SOR(matrix,te,source,count)
    integer,intent(inout)::count
    double precision,intent(in)::matrix(ix,jx,kx,n_hc)
    double precision,intent(inout)::te(ix,jx,kx),source(ix,jx,kx)
    double precision residual(ix,jx,kx)
    integer i,j,k,n
    double precision total_resi,total_source,tmp
    residual(is:ie,js:je,ks:ke)=source(is:ie,js:je,ks:ke) 
    do n=1,n_hc
       residual(is:ie,js:je,ks:ke)=residual(is:ie,js:je,ks:ke) &
            -matrix(is:ie,js:je,ks:ke,n)* &
            te(is+di(n):ie+di(n),js+dj(n):je+dj(n),ks+dk(n):ke:dk(n))
    enddo
    total_resi=sum(residual(is:ie,js:je,ks:ke)**2)
    total_source=sum(source(is:ie,js:je,ks:ke)**2)
    count=0    
    select case(hc_split)
    case(1)
       do while(sqrt(total_resi/total_source).gt.error.and.count.le.hc_max)
          do k=ks,ke; do j=js,je;do i=is,ie
             residual(i,j,k)=source(i,j,k)-(matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k))   
             te(i,j,k)=te(i,j,k)+sor_omega*residual(i,j,k)/matrix(i,j,k,1)  
          enddo;enddo;enddo          
          total_resi=sum(residual(is:ie,js:je,ks:ke)**2)
          count=count+1
       enddo
    case(2)
       do while(sqrt(total_resi/total_source).gt.error.and.count.le.hc_max)
          do k=ks,ke; do j=js,je;do i=is,ie
             residual(i,j,k)=source(i,j,k)-(matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k) &
                  +matrix(i,j,k,6)*te(i-1,j-1,k)+matrix(i,j,k,7)*te(i+1,j-1,k) &
                  +matrix(i,j,k,8)*te(i-1,j+1,k)+matrix(i,j,k,9)*te(i+1,j+1,k))
             te(i,j,k)=te(i,j,k)+sor_omega*residual(i,j,k)/matrix(i,j,k,1)   
          enddo;enddo;enddo
          total_resi=sum(residual(is:ie,js:je,ks:ke)**2)
          count=count+1
       enddo       
    end select
  end subroutine SOR


  subroutine BiCGstab(matrix,te,source,count)
    integer,intent(inout)::count
    double precision,intent(in)::matrix(ix,jx,kx,n_hc)
    double precision,intent(inout)::te(ix,jx,kx),source(ix,jx,kx)
    double precision residual(ix,jx,kx)
    integer i,j,k,n
    double precision total_resi,total_source,tmp
    double precision alpha,alphau,alphad
    double precision beta,betau,betad
    double precision omega,omegau,omegad
    double precision pk(ix,jx,kx),r0(ix,jx,kx),Apk(ix,jx,kx),tkp(ix,jx,kx)
    double precision atk
    !    residual(is:ie,js:je,ks:ke)=source(is:ie,js:je,ks:ke)
    residual=source
    do n=1,n_hc
       residual(is:ie,js:je,ks:ke)=residual(is:ie,js:je,ks:ke) &
            -matrix(is:ie,js:je,ks:ke,n)* &
            te(is+di(n):ie+di(n),js+dj(n):je+dj(n),ks+dk(n):ke:dk(n))
    enddo

    total_resi=sum(residual(is:ie,js:je,ks:ke)**2)
    total_source=sum(source(is:ie,js:je,ks:ke)**2)
    count=0    
    r0=residual
    pk=r0
    apk=0
    
!    print *,"TOTAL_RESI",total_resi,total_source
!    stop
    select case(hc_split)
    case(1)
       do while(sqrt(total_resi/total_source).gt.error.and.count.le.hc_max)
          alphau=0.0d0;alphad=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             apk(i,j,k)=matrix(i,j,k,1)*pk(i,j,k) &
                  +matrix(i,j,k,2)*pk(i-1,j,k)+matrix(i,j,k,3)*pk(i+1,j,k)  &
                  +matrix(i,j,k,4)*pk(i,j-1,k)+matrix(i,j,k,5)*pk(i,j+1,k)
             alphau=alphau+(residual(i,j,k)*r0(i,j,k))
             alphad=alphad+(r0(i,j,k)*apk(i,j,k))
          enddo;enddo;enddo
          alpha=alphau/alphad
          tkp=residual-alpha*apk
          omegau=0.0d0;omegad=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             atk=matrix(i,j,k,1)*tkp(i,j,k) &
                  +matrix(i,j,k,2)*tkp(i-1,j,k)+matrix(i,j,k,3)*tkp(i+1,j,k)  &
                  +matrix(i,j,k,4)*tkp(i,j-1,k)+matrix(i,j,k,5)*tkp(i,j+1,k)
             omegau=omegau+tkp(i,j,k)*atk
             omegad=omegad+atk*atk
          enddo;enddo;enddo
          omega=omegau/omegad
          te=te+alpha*pk+omega*tkp
          betau=0.0d0;betad=0.0d0;total_resi=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             residual(i,j,k)=source(i,j,k)-( &
                  matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k))  
             betau=betau+r0(i,j,k)*residual(i,j,k)
             total_resi=total_resi+residual(i,j,k)*residual(i,j,k)
          enddo;enddo;enddo
          beta=betau/alphau*alpha/omega
          pk=residual+beta*(pk-omega*apk)    
          count=count+1
       enddo
    case(2)
       do while(sqrt(total_resi/total_source).gt.error.and.count.le.hc_max)
          alphau=0.0d0;alphad=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             apk(i,j,k)=matrix(i,j,k,1)*pk(i,j,k) &
                  +matrix(i,j,k,2)*pk(i-1,j,k)+matrix(i,j,k,3)*pk(i+1,j,k)  &
                  +matrix(i,j,k,4)*pk(i,j-1,k)+matrix(i,j,k,5)*pk(i,j+1,k) &
                  +matrix(i,j,k,6)*pk(i-1,j-1,k)+matrix(i,j,k,7)*pk(i+1,j-1,k) &
                  +matrix(i,j,k,8)*pk(i-1,j+1,k)+matrix(i,j,k,9)*pk(i+1,j+1,k)
             alphau=alphau+(residual(i,j,k)*r0(i,j,k))
             alphad=alphad+(r0(i,j,k)*apk(i,j,k))
          enddo;
       enddo;enddo
          alpha=alphau/alphad
!          print *,count,total_resi,total_source,alphau,alphad,alpha   
          tkp=residual-alpha*apk
          omegau=0.0d0;omegad=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             atk=matrix(i,j,k,1)*tkp(i,j,k) &
                  +matrix(i,j,k,2)*tkp(i-1,j,k)+matrix(i,j,k,3)*tkp(i+1,j,k)  &
                  +matrix(i,j,k,4)*tkp(i,j-1,k)+matrix(i,j,k,5)*tkp(i,j+1,k) &
                  +matrix(i,j,k,6)*tkp(i-1,j-1,k)+matrix(i,j,k,7)*tkp(i+1,j-1,k)&
                  +matrix(i,j,k,8)*tkp(i-1,j+1,k)+matrix(i,j,k,9)*tkp(i+1,j+1,k)
             omegau=omegau+tkp(i,j,k)*atk
             omegad=omegad+atk*atk
          enddo;enddo;enddo
          omega=omegau/omegad
          te=te+alpha*pk+omega*tkp
          betau=0.0d0;betad=0.0d0;total_resi=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             residual(i,j,k)=source(i,j,k)-( &
                  matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k)  &
                  +matrix(i,j,k,6)*te(i-1,j-1,k)+matrix(i,j,k,7)*te(i+1,j-1,k) &
                  +matrix(i,j,k,8)*te(i-1,j+1,k)+matrix(i,j,k,9)*te(i+1,j+1,k))
             betau=betau+r0(i,j,k)*residual(i,j,k)
             total_resi=total_resi+residual(i,j,k)*residual(i,j,k)
          enddo;enddo;enddo
          beta=betau/alphau*alpha/omega
          pk=residual+beta*(pk-omega*apk)          
          count=count+1
!          print *,count,total_resi,alpha,omega,beta,betau
       enddo
!       stop
    end select

  end subroutine BiCGstab



  subroutine BiCGstabMPI(matrix,te,source,count)
    integer,intent(inout)::count
    double precision,intent(in)::matrix(ix,jx,kx,n_hc)
    double precision,intent(inout)::te(ix,jx,kx),source(ix,jx,kx)
    double precision residual(ix,jx,kx)
    integer i,j,k,n
    double precision total_resi,total_source,tmp
    double precision alpha,alphas(2),g_alphas(2)
    double precision beta,betas(2),g_betas(2)
    double precision omega,omegas(2),g_omegas(2)
    double precision pk(ix,jx,kx),r0(ix,jx,kx),Apk(ix,jx,kx),tkp(ix,jx,kx)
    double precision atk,g_total_resi,g_total_source
    integer,parameter::g_margin=1
    !    residual(is:ie,js:je,ks:ke)=source(is:ie,js:je,ks:ke)
    residual=source

    do n=1,n_hc
       residual(is:ie,js:je,ks:ke)=residual(is:ie,js:je,ks:ke) &
            -matrix(is:ie,js:je,ks:ke,n)* &
            te(is+di(n):ie+di(n),js+dj(n):je+dj(n),ks+dk(n):ke:dk(n))       
!       print *,n,matrix(is:ie,js:js,ks:ke,n)
    enddo


    
    total_resi=sum(residual(is:ie,js:je,ks:ke)**2)
    total_source=sum(source(is:ie,js:je,ks:ke)**2)

    if(total_prc.ge.2) then
       call mpi_allreduce(total_resi,g_total_resi,1,mpi_double_precision, &
            MPI_SUM,mpi_comm_world,ierr)
       call mpi_allreduce(total_source,g_total_source,1,mpi_double_precision, &
            MPI_SUM,mpi_comm_world,ierr)
       call mpiBND(residual,g_margin)
    endif
    count=0    
    r0=residual
    pk=r0
    apk=0
    select case(hc_split)
    case(1)
       do while(sqrt(g_total_resi/g_total_source).gt.error.and.count.le.hc_max)
          alphas(:)=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             apk(i,j,k)=matrix(i,j,k,1)*pk(i,j,k) &
                  +matrix(i,j,k,2)*pk(i-1,j,k)+matrix(i,j,k,3)*pk(i+1,j,k)  &
                  +matrix(i,j,k,4)*pk(i,j-1,k)+matrix(i,j,k,5)*pk(i,j+1,k)
             alphas(1)=alphas(1)+(residual(i,j,k)*r0(i,j,k))
             alphas(2)=alphas(2)+(r0(i,j,k)*apk(i,j,k))
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpi_allreduce(alphas,g_alphas,2,mpi_double_precision,&
                  mpi_sum,mpi_comm_world,ierr)
             call mpiBND(apk,g_margin)
          endif
          alpha=g_alphas(2)/g_alphas(1)
          tkp=residual-alpha*apk
          omegas(:)=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             atk=matrix(i,j,k,1)*tkp(i,j,k) &
                  +matrix(i,j,k,2)*tkp(i-1,j,k)+matrix(i,j,k,3)*tkp(i+1,j,k)  &
                  +matrix(i,j,k,4)*tkp(i,j-1,k)+matrix(i,j,k,5)*tkp(i,j+1,k)
             omegas(1)=omegas(1)+tkp(i,j,k)*atk
             omegas(2)=omegas(2)+atk*atk
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpi_allreduce(omegas,g_omegas,2,mpi_double_precision,&
                  mpi_sum,mpi_comm_world,ierr)
!             call mpiBND(atk,2)
          endif          
          omega=omegas(2)/omegas(1)
          te=te+alpha*pk+omega*tkp
          if(total_prc.ge.2) then
             call mpiBND(te,g_margin)
          endif          
          betas(:)=0.0d0;total_resi=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             residual(i,j,k)=source(i,j,k)-( &
                  matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k))  
             betas(1)=betas(1)+r0(i,j,k)*residual(i,j,k)
             total_resi=total_resi+residual(i,j,k)*residual(i,j,k)
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpiBND(residual,g_margin)
             call mpi_allreduce(betas,g_betas,2,mpi_double_precision,mpi_sum &
     ,mpi_comm_world,ierr)
             call mpi_allreduce(total_resi,g_total_resi,1, &
                  mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
          endif          
          beta=g_betas(1)/g_alphas(1)*alpha/omega
          pk=residual+beta*(pk-omega*apk)    
          count=count+1

       enddo
       
    case(2)
       do while(sqrt(g_total_resi/g_total_source).gt.error.and.count.le.hc_max)
          alphas(:)=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             apk(i,j,k)=matrix(i,j,k,1)*pk(i,j,k) &
                  +matrix(i,j,k,2)*pk(i-1,j,k)+matrix(i,j,k,3)*pk(i+1,j,k)  &
                  +matrix(i,j,k,4)*pk(i,j-1,k)+matrix(i,j,k,5)*pk(i,j+1,k) &
                  +matrix(i,j,k,6)*pk(i-1,j-1,k)+matrix(i,j,k,7)*pk(i+1,j-1,k) &
                  +matrix(i,j,k,8)*pk(i-1,j+1,k)+matrix(i,j,k,9)*pk(i+1,j+1,k)
             alphas(1)=alphas(1)+(residual(i,j,k)*r0(i,j,k))
             alphas(2)=alphas(2)+(r0(i,j,k)*apk(i,j,k))
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpi_allreduce(alphas,g_alphas,2,mpi_double_precision,&
                  mpi_sum,mpi_comm_world,ierr)
             call mpiBND(apk,2)
          endif
          alpha=g_alphas(1)/g_alphas(2)
          
          tkp=residual-alpha*apk
          omegas(:)=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             atk=matrix(i,j,k,1)*tkp(i,j,k) &
                  +matrix(i,j,k,2)*tkp(i-1,j,k)+matrix(i,j,k,3)*tkp(i+1,j,k)  &
                  +matrix(i,j,k,4)*tkp(i,j-1,k)+matrix(i,j,k,5)*tkp(i,j+1,k) &
                  +matrix(i,j,k,6)*tkp(i-1,j-1,k)+matrix(i,j,k,7)*tkp(i+1,j-1,k)&
                  +matrix(i,j,k,8)*tkp(i-1,j+1,k)+matrix(i,j,k,9)*tkp(i+1,j+1,k)
             omegas(1)=omegas(1)+tkp(i,j,k)*atk
             omegas(2)=omegas(2)+atk*atk
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpi_allreduce(omegas,g_omegas,2,mpi_double_precision,&
                  mpi_sum,mpi_comm_world,ierr)
          endif
          omega=omegas(1)/omegas(2)
          te=te+alpha*pk+omega*tkp
          if(total_prc.ge.2) then
             call mpiBND(te,2)
          endif
          betas(:)=0.0d0;total_resi=0.0d0
          do k=ks,ke; do j=js,je; do i=is,ie
             residual(i,j,k)=source(i,j,k)-( &
                  matrix(i,j,k,1)*te(i,j,k) &
                  +matrix(i,j,k,2)*te(i-1,j,k)+matrix(i,j,k,3)*te(i+1,j,k)  &
                  +matrix(i,j,k,4)*te(i,j-1,k)+matrix(i,j,k,5)*te(i,j+1,k)  &
                  +matrix(i,j,k,6)*te(i-1,j-1,k)+matrix(i,j,k,7)*te(i+1,j-1,k) &
                  +matrix(i,j,k,8)*te(i-1,j+1,k)+matrix(i,j,k,9)*te(i+1,j+1,k))
             
             betas(1)=betas(1)+r0(i,j,k)*residual(i,j,k)
             total_resi=total_resi+residual(i,j,k)*residual(i,j,k)
          enddo;enddo;enddo
          if(total_prc.ge.2) then
             call mpiBND(residual,2)
             call mpi_allreduce(betas,g_betas,2,mpi_double_precision,mpi_sum &
                  ,mpi_comm_world,ierr)
             call mpi_allreduce(total_resi,g_total_resi,1, &
                  mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

          endif
          beta=g_betas(1)/g_alphas(1)*alpha/omega
          pk=residual+beta*(pk-omega*apk)    
          count=count+1
!          if(my_rank.eq.0)print*,count,g_total_resi,alpha,omega,beta,g_betas(1)
       enddo
    end select
!    stop
  end subroutine BiCGstabMPI
  
  subroutine set_matrix_boundary(matrix)
    double precision,intent(inout):: matrix(ix,jx,kx,n_hc)
    select case(ndim)
    case(1)
       
    case(2)
       select case(hc_split)
       case(1)
          if(neighbor(1).eq.-1) then
             matrix(is,:,:,1)=matrix(is,:,:,2)+matrix(is,:,:,1)
             matrix(is,:,:,2)=0.0d0       
             matrix(ie,:,:,1)=matrix(ie,:,:,3)+matrix(ie,:,:,1)
             matrix(ie,:,:,3)=0.0d0
          endif
          if(neighbor(2).eq.-1) then
             matrix(:,js,:,1)=matrix(:,js,:,4)+matrix(:,js,:,1)
             matrix(:,js,:,4)=0.0d0       
             matrix(:,je,:,1)=matrix(:,je,:,5)+matrix(:,je,:,1)
             matrix(:,je,:,5)=0.0d0              
          endif
       case(2)
          !x-minus
          if(neighbor(1).eq.-1) then
             matrix(is,:,:,1)=matrix(is,:,:,2)+matrix(is,:,:,1)
             matrix(is,:,:,2)=0.0d0       
             matrix(is,:,:,4)=matrix(is,:,:,4)+matrix(is,:,:,6)
             matrix(is,:,:,6)=0.0d0       
             matrix(is,:,:,5)=matrix(is,:,:,5)+matrix(is,:,:,8)
             matrix(is,:,:,8)=0.0d0       
          endif
          !x-plus
          if(neighbor(2).eq.-1) then
             matrix(ie,:,:,1)=matrix(ie,:,:,3)+matrix(ie,:,:,1)
             matrix(ie,:,:,3)=0.0d0              
             matrix(ie,:,:,4)=matrix(ie,:,:,4)+matrix(ie,:,:,7)
             matrix(ie,:,:,7)=0.0d0              
             matrix(ie,:,:,5)=matrix(ie,:,:,5)+matrix(ie,:,:,9)
             matrix(ie,:,:,9)=0.0d0              
          endif
          !y-minus
          if(neighbor(3).eq.-1) then
             matrix(:,js,:,1)=matrix(:,js,:,4)+matrix(:,js,:,1)
             matrix(:,js,:,4)=0.0d0       
             matrix(:,js,:,2)=matrix(:,js,:,6)+matrix(:,js,:,2)
             matrix(:,js,:,6)=0.0d0       
             matrix(:,js,:,3)=matrix(:,js,:,7)+matrix(:,js,:,3)
             matrix(:,js,:,7)=0.0d0       
          endif
          if(neighbor(4).eq.-1) then
             matrix(:,je,:,1)=matrix(:,je,:,5)+matrix(:,je,:,1)
             matrix(:,je,:,5)=0.0d0              
             matrix(:,je,:,2)=matrix(:,je,:,2)+matrix(:,je,:,8)
             matrix(:,je,:,8)=0.0d0              
             matrix(:,je,:,3)=matrix(:,je,:,3)+matrix(:,je,:,9)
             matrix(:,je,:,9)=0.0d0              
          endif
       end select
    case(3)
    end select
    

  end subroutine set_matrix_boundary

  subroutine test_matrix(matrix,source,te)
    double precision,intent(inout):: matrix(ix,jx,n_hc),source(ix,jx),te(ix,jx)
    double precision eta
    character*2 c_sch
    character*4 c_ome
    integer i,j
    write(c_sch,"(i2.2)")hc_sch
    write(c_ome,"(f4.2)")sor_omega    
    open(hc_unit,file=trim("tmp_hc_"//c_sch//"_"//c_ome//".txt"), &
         status="replace",form="formatted")
    
    eta=-kapper
!    eta=1.5d0
    do j=js,je  
       do i=is,ie
          matrix(i,j,2)=eta*2
          matrix(i,j,3)=eta*2
          matrix(i,j,4)=eta
          matrix(i,j,5)=eta
          matrix(i,j,1)=1.0d0-matrix(i,j,2)-matrix(i,j,3)-matrix(i,j,4)-matrix(i,j,5)     
          source(i,j)=1.0d0+10.0d0*exp(-((x(i)**2)/0.1d0))
          te(i,j)=1.0d0
       enddo
    enddo
  end subroutine test_matrix
  subroutine set_matrix_source_1D(matrix,source,ro,te)
    double precision,intent(inout):: matrix(ix,n_hc),source(ix)
    double precision,intent(in)::ro(ix),te(ix)
    double precision eta
    integer i
    source=te
    do i=is,ie
       eta=kapper*dt*(gm-1.0d0)/ro(i)/dxc(i)
       matrix(i,2)=-eta*(0.5d0*(te(i-1)+te(i)))**2.5/dx(i-1)
       matrix(i,3)=-eta*(0.5d0*(te(i+1)+te(i)))**2.5/dx(i)
       matrix(i,1)=1.0d0-matrix(i,2)-matrix(i,3)
    enddo    
  end subroutine set_matrix_source_1D

  subroutine set_matrix_source_2D(matrix,source,ro,te)
    double precision,intent(inout):: matrix(ix,jx,n_hc),source(ix,jx)
    double precision,intent(in)::ro(ix,jx),te(ix,jx)
    double precision eta
    integer i,j
    source=te
    do j=js,je;do i=is,ie
!       eta=dt/(dx(i)*dx(i)*kapper
       eta=kapper*dt*(gm-1.0d0)/ro(i,j)*gm
       matrix(i,j,2)=-eta*(0.5d0*(te(i-1,j)+te(i,j)))**2.5/(dx(i-1)*dxc(i))
       matrix(i,j,3)=-eta*(0.5d0*(te(i+1,j)+te(i,j)))**2.5/(dx(i)*dxc(i))
       matrix(i,j,4)=-eta*(0.5d0*(te(i,j-1)+te(i,j)))**2.5/(dy(j-1)*dyc(j))
       matrix(i,j,5)=-eta*(0.5d0*(te(i,j+1)+te(i,j)))**2.5/(dy(j)*dyc(j))
       matrix(i,j,1)=1.0d0-matrix(i,j,2)-matrix(i,j,3)&
            -matrix(i,j,4)-matrix(i,j,5)
    enddo;enddo
  end subroutine set_matrix_source_2D
  
  subroutine set_matrix_source_MHD2D5(matrix,source,ro,te,bx,by,bz)
    double precision,intent(inout):: matrix(ix,jx,n_hc),source(ix,jx)
    double precision,intent(in)::ro(ix,jx),te(ix,jx)
    double precision,intent(in)::bx(ix,jx),by(ix,jx),bz(ix,jx)
    double precision bypx,bymx,bxpy,bxmy,bxp,byp,bxm,bym,b2n,b2e,b2s,b2w,eta
    integer i,j
    !set matrix and source term-----------------------------------------
    eta=gm*(gm-1.0d0)*dt*kapper
    do j=js,je ; do i=is,ie
!       bypx=(by(i,j-1)+by(i,j)+by(i+1,j)+by(i+1,j-1))/4.0d0
!       bymx=(by(i-1,j-1)+by(i-1,j)+by(i,j)+by(i,j-1))/4.0d0
!       bxpy=(bx(i-1,j)+bx(i,j)+bx(i,j+1)+bx(i-1,j+1))/4.0d0
!       bxmy=(bx(i-1,j-1)+bx(i,j-1)+bx(i,j)+bx(i-1,j))/4.0d0
!       bxp=bx(i,j) ; bxm=bx(i-1,j)
!       byp=by(i,j) ; bym=by(i,j-1)
       bxp=0.5d0*(bx(i,j)+bx(i+1,j))
       bxm=0.5d0*(bx(i-1,j)+bx(i,j))
       bxpy=0.5d0*(bx(i,j+1)+bx(i,j))
       bxmy=0.5d0*(bx(i,j-1)+bx(i,j))

       bypx=0.5d0*(by(i,j)+by(i+1,j))
       bymx=0.5d0*(by(i,j)+by(i-1,j))
       byp=0.5d0*(by(i,j)+by(i,j+1))
       bym=0.5d0*(by(i,j-1)+by(i,j))
       


       b2n=byp**2+bxpy**2+((bz(i,j)+bz(i,j+1))/2.0d0)**2+epsi
       b2e=bxp**2+bypx**2+((bz(i,j)+bz(i+1,j))/2.0d0)**2+epsi
       b2s=bym**2+bxmy**2+((bz(i,j)+bz(i,j-1))/2.0d0)**2+epsi
       b2w=bxm**2+bymx**2+((bz(i,j)+bz(i-1,j))/2.0d0)**2+epsi       
       source(i,j)=te(i,j)+eta/ro(i,j)/dxc(i)/dyc(j)*(  &
            ((te(i+1,j)+te(i,j))/2.0)**2.5*bxp*bypx/b2e* &
            ((te(i,j+1)+te(i+1,j+1)-te(i,j-1)-te(i+1,j-1))/4.0d0)  - &
            (((te(i,j)+te(i-1,j))/2.0)**2.5*bxm*bymx/b2w*  &
            ((te(i-1,j+1)+te(i,j+1)-te(i-1,j-1)-te(i,j-1))/4.0d0))  + &
            ((te(i,j+1)+te(i,j))/2.0)**2.5*byp*bxpy/b2n*  &
            ((te(i+1,j)+te(i+1,j+1)-te(i-1,j)-te(i-1,j+1))/4.0d0)  - &
            (((te(i,j)+te(i,j-1))/2.0)**2.5*bym*bxmy/b2s*  &
            ((te(i+1,j-1)+te(i+1,j)-te(i-1,j-1)-te(i-1,j))/4.0d0)))
       !set Matrix
       matrix(i,j,2)=-eta/(dx(i-1)*dxc(i))/ro(i,j) &
            *((te(i-1,j)+te(i,j))/2)**2.5d0*bxm**2/b2w
       matrix(i,j,3)=-eta/(dx(i)*dxc(i))/ro(i,j)  &
            *((te(i+1,j)+te(i,j))/2)**2.5d0*bxp**2/b2e
       matrix(i,j,4)=-eta/(dy(j-1)*dyc(j))/ro(i,j) &
            *((te(i,j)+te(i,j-1))/2)**2.5d0*bym**2/b2s
       matrix(i,j,5)=-eta*dt/(dy(j)*dyc(j))/ro(i,j) &
            *((te(i,j)+te(i,j+1))/2)**2.5d0*byp**2/b2n
       matrix(i,j,1)=1.0-(matrix(i,j,2)+matrix(i,j,3)+matrix(i,j,4)+matrix(i,j,5))
    enddo;enddo
    
    !--------------------------------------------------------------------
  end subroutine set_matrix_source_MHD2D5

  subroutine set_matrix_source_MHD2D9(matrix,source,ro,te,bx,by,bz)
    double precision,intent(inout):: matrix(ix,jx,kx,n_hc),source(ix,jx,kx)
    double precision,intent(in)::ro(ix,jx,kx),te(ix,jx,kx)
    double precision,intent(in)::bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision ::bxl,bxr,bxd,bxu,byl,byr,byd,byu,tnorm,nyl,nyr,nyu,nyd
    double precision ::bxxl,bxxr,byyd,byyu,bxyl,bxyr,byxd,byxu,nxl,nxr,nxu,nxd
    double precision ::bzl,bzr,bzd,bzu,res0,eta
    double precision ::te_xm,te_xp,te_ym,te_yp
    double precision ::te_xm25,te_xp25,te_ym25,te_yp25
    integer i,j,k
    !set matrix and source term-----------------------------------------
    source=te
    
    do k=ks,ke;do j=js,je ;do i=is,ie          

       eta=kapper*dt*(gm-1.0d0)/ro(i,j,k)*gm
!       bxl=bx(i-1,j,k) ; bxr=bx(i,j,k)
!       bxd=(bx(i-1,j-1,k)+bx(i,j-1,k)+bx(i-1,j,k)+bx(i,j,k))/4.0d0
!       bxu=(bx(i-1,j,k)+bx(i,j,k)+bx(i-1,j+1,k)+bx(i,j+1,k))/4.0d0

!       byl=(by(i-1,j-1,k)+by(i,j-1,k)+by(i-1,j,k)+by(i,j,k))/4.0d0
!       byr=(by(i,j-1,k)+by(i+1,j-1,k)+by(i,j,k)+by(i+1,j,k))/4.0d0
!       byd=by(i,j-1,k) ; byu=by(i,j,k)
       
       bxl=(bx(i-1,j,k)+bx(i,j,k))/2.0d0
       bxr=(bx(i+1,j,k)+bx(i,j,k))/2.0d0
       bxd=(bx(i,j-1,k)+bx(i,j,k))/2.0d0
       bxu=(bx(i,j+1,k)+bx(i,j,k))/2.0d0

       byl=(by(i-1,j,k)+by(i,j,k))/2.0d0
       byr=(by(i+1,j,k)+by(i,j,k))/2.0d0
       byd=(by(i,j-1,k)+by(i,j,k))/2.0d0
       byu=(by(i,j+1,k)+by(i,j,k))/2.0d0

       bzl=(bz(i-1,j,k)+bz(i,j,k))/2.0d0
       bzr=(bz(i+1,j,k)+bz(i,j,k))/2.0d0
       bzd=(bz(i,j-1,k)+bz(i,j,k))/2.0d0
       bzu=(bz(i,j+1,k)+bz(i,j,k))/2.0d0
       
       nxl=bxl/(max(sqrt(bxl**2+byl**2+bzl**2),epsi))
       nxr=bxr/(max(sqrt(bxr**2+byr**2+bzr**2),epsi))
       nxd=bxd/(max(sqrt(bxd**2+byd**2+bzd**2),epsi))
       nxu=bxu/(max(sqrt(bxu**2+byu**2+bzu**2),epsi))

       nyl=byl/(max(sqrt(bxl**2+byl**2+bzl**2),epsi))
       nyr=byr/(max(sqrt(bxr**2+byr**2+bzr**2),epsi))
       nyd=byd/(max(sqrt(bxd**2+byd**2+bzd**2),epsi))
       nyu=byu/(max(sqrt(bxu**2+byu**2+bzu**2),epsi))   

       te_xm=0.5d0*(te(i,j,k)+te(i-1,j,k))
       te_xp=0.5d0*(te(i,j,k)+te(i+1,j,k))
       te_ym=0.5d0*(te(i,j,k)+te(i,j-1,k))
       te_yp=0.5d0*(te(i,j,k)+te(i,j+1,k))
       
       te_xm25=te_xm*te_xm*sqrt(te_xm)
       te_xp25=te_xp*te_xp*sqrt(te_xp)
       te_ym25=te_ym*te_ym*sqrt(te_ym)
       te_yp25=te_yp*te_yp*sqrt(te_yp)
    
       bxxl=-eta*te_xm25*nxl*nxl/dxc(i)/dx(i-1)
       bxxr=-eta*te_xp25*nxr*nxr/dxc(i)/dx(i)
       byyd=-eta*te_ym25*nyd*nyd/dyc(j)/dy(j-1)
       byyu=-eta*te_yp25*nyu*nyu/dyc(j)/dy(j)

       bxyl=-eta/4.0d0/dxc(i)/dyc(j)*nxl*nyl*te_xm25
       bxyr=-eta/4.0d0/dxc(i)/dyc(j)*nxr*nyr*te_xp25
       byxd=-eta/4.0d0/dxc(i)/dyc(j)*nyd*nxd*te_ym25
       byxu=-eta/4.0d0/dxc(i)/dyc(j)*nyu*nxu*te_yp25

!       bxxl=-eta*(0.5d0*(te(i,j,k)+te(i-1,j,k))* &
!            (0.5d0*(te(i,j,k)+te(i-1,j,k)))* &
!            sqrt(0.5d0*(te(i,j,k)+te(i-1,j,k)))) &
!            *nxl*nxl/dxc(i)/dx(i-1)
!       bxxr=-eta*(0.5d0*(te(i,j,k)+te(i+1,j,k)))**2.5d0*nxr**2/dxc(i)/dx(i)
!       byyd=-eta*(0.5d0*(te(i,j,k)+te(i,j-1,k)))**2.5d0*nyd**2/dyc(j)/dy(j-1)
!       byyu=-eta*(0.5d0*(te(i,j,k)+te(i,j+1,k)))**2.5d0*nyu**2/dyc(j)/dy(j)
       
!       bxyl=-eta/4.0d0/dxc(i)/dyc(j)*nxl*nyl*(0.5d0*(te(i,j,k)+te(i-1,j,k)))**2.5d0
!       bxyr=-eta/4.0d0/dxc(i)/dyc(j)*nxr*nyr*(0.5d0*(te(i,j,k)+te(i+1,j,k)))**2.5d0
!       byxd=-eta/4.0d0/dxc(i)/dyc(j)*nyd*nxd*(0.5d0*(te(i,j,k)+te(i,j-1,k)))**2.5d0
!       byxu=-eta/4.0d0/dxc(i)/dyc(j)*nyu*nxu*(0.5d0*(te(i,j,k)+te(i,j+1,k)))**2.5d0

       matrix(i,j,k,2)=bxxl-(byxu-byxd)
       matrix(i,j,k,3)=bxxr+(byxu-byxd)
       matrix(i,j,k,4)=byyd-(bxyr-bxyl)
       matrix(i,j,k,5)=byyu+(bxyr-bxyl)

       matrix(i,j,k,6)=bxyl+byxd
       matrix(i,j,k,7)=-(bxyr+byxd)
       matrix(i,j,k,8)=-(bxyl+byxu)
       matrix(i,j,k,9)=bxyr+byxu
       matrix(i,j,k,1)=1.0-(matrix(i,j,k,2)+matrix(i,j,k,3)+matrix(i,j,k,4)+matrix(i,j,k,5))

    enddo;enddo;enddo

    !--------------------------------------------------------------------
  end subroutine set_matrix_source_MHD2D9


  subroutine set_matrix_source_MHD2D9_limiter(matrix,source,ro,te,bx,by,bz)
    double precision,intent(inout):: matrix(ix,jx,kx,n_hc),source(ix,jx,kx)
    double precision,intent(in)::ro(ix,jx,kx),te(ix,jx,kx)
    double precision,intent(in)::bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision ::bxl,bxr,bxd,bxu,byl,byr,byd,byu,tnorm,nyl,nyr,nyu,nyd
    double precision ::bxxl,bxxr,byyd,byyu,bxyl,bxyr,byxd,byxu,nxl,nxr,nxu,nxd
    double precision ::bzl,bzr,bzd,bzu,res0,eta
    double precision ::te_xm,te_xp,te_ym,te_yp
    double precision ::te_xm25,te_xp25,te_ym25,te_yp25
    double precision ::lim_l,lim_r,lim_d,lim_u
    double precision ::F_l,F_r,F_d,F_u
    double precision ::F_sat_l,F_sat_r,F_sat_d,F_sat_u
    double precision ::dTdx_l,dTdx_r,dTdx_d,dTdx_u
    double precision ::dTdy_l,dTdy_r,dTdy_d,dTdy_u
    integer i,j,k
    !set matrix and source term-----------------------------------------
    source=te
    
    do k=ks,ke;do j=js,je ;do i=is,ie          

       eta=kapper*dt*(gm-1.0d0)/ro(i,j,k)*gm
       
       bxl=(bx(i-1,j,k)+bx(i,j,k))/2.0d0
       bxr=(bx(i+1,j,k)+bx(i,j,k))/2.0d0
       bxd=(bx(i,j-1,k)+bx(i,j,k))/2.0d0
       bxu=(bx(i,j+1,k)+bx(i,j,k))/2.0d0

       byl=(by(i-1,j,k)+by(i,j,k))/2.0d0
       byr=(by(i+1,j,k)+by(i,j,k))/2.0d0
       byd=(by(i,j-1,k)+by(i,j,k))/2.0d0
       byu=(by(i,j+1,k)+by(i,j,k))/2.0d0

       bzl=(bz(i-1,j,k)+bz(i,j,k))/2.0d0
       bzr=(bz(i+1,j,k)+bz(i,j,k))/2.0d0
       bzd=(bz(i,j-1,k)+bz(i,j,k))/2.0d0
       bzu=(bz(i,j+1,k)+bz(i,j,k))/2.0d0
       
       nxl=bxl/(max(sqrt(bxl**2+byl**2+bzl**2),epsi))
       nxr=bxr/(max(sqrt(bxr**2+byr**2+bzr**2),epsi))
       nxd=bxd/(max(sqrt(bxd**2+byd**2+bzd**2),epsi))
       nxu=bxu/(max(sqrt(bxu**2+byu**2+bzu**2),epsi))

       nyl=byl/(max(sqrt(bxl**2+byl**2+bzl**2),epsi))
       nyr=byr/(max(sqrt(bxr**2+byr**2+bzr**2),epsi))
       nyd=byd/(max(sqrt(bxd**2+byd**2+bzd**2),epsi))
       nyu=byu/(max(sqrt(bxu**2+byu**2+bzu**2),epsi))   

       te_xm=0.5d0*(te(i,j,k)+te(i-1,j,k))
       te_xp=0.5d0*(te(i,j,k)+te(i+1,j,k))
       te_ym=0.5d0*(te(i,j,k)+te(i,j-1,k))
       te_yp=0.5d0*(te(i,j,k)+te(i,j+1,k))
       

!    kc_sat = 5.d0*ro*te*sqrt(te) !! for saturation flux (defined at i,j,k)
       F_sat_l=2.5d0*(ro(i-1,j,k)+ro(i,j,k))*te_xm*sqrt(te_xm)
       F_sat_r=2.5d0*(ro(i+1,j,k)+ro(i,j,k))*te_xp*sqrt(te_xp)
       F_sat_d=2.5d0*(ro(i,j-1,k)+ro(i,j,k))*te_ym*sqrt(te_ym)
       F_sat_u=2.5d0*(ro(i,j+1,k)+ro(i,j,k))*te_yp*sqrt(te_yp)       
       
       te_xm25=te_xm*te_xm*sqrt(te_xm)
       te_xp25=te_xp*te_xp*sqrt(te_xp)
       te_ym25=te_ym*te_ym*sqrt(te_ym)
       te_yp25=te_yp*te_yp*sqrt(te_yp)

       dTdx_l=(te(i,j,k)-te(i-1,j,k))/dx(i-1)
       dTdx_r=(te(i+1,j,k)-te(i,j,k))/dx(i)
       dTdx_d=(te(i+1,j-1,k)+te(i+1,j,k)-te(i-1,j-1,k)-te(i-1,j,k))/(4.0d0*dxc(i))
       dTdx_u=(te(i+1,j+1,k)+te(i+1,j,k)-te(i-1,j+1,k)-te(i-1,j,k))/(4.0d0*dxc(i))

       dTdy_l=(te(i-1,j+1,k)+te(i,j+1,k)-te(i-1,j-1,k)-te(i,j-1,k))/(4.0d0*dyc(j))
       dTdy_r=(te(i+1,j+1,k)+te(i,j+1,k)-te(i+1,j-1,k)-te(i,j-1,k))/(4.0d0*dyc(j))
       dTdy_d=(te(i,j,k)-te(i,j-1,k))/dy(j-1)
       dTdy_u=(te(i,j+1,k)-te(i,j,k))/dy(j)
    
       F_l=te_xm25*kapper*sqrt(dTdx_l**2+dTdy_l**2)
       F_r=te_xp25*kapper*sqrt(dTdx_r**2+dTdy_r**2)
       F_d=te_ym25*kapper*sqrt(dTdx_d**2+dTdy_d**2)
       F_u=te_yp25*kapper*sqrt(dTdx_u**2+dTdy_u**2)
       
       lim_l=max(0.0d0,min(1.0d0,1.0d0-0.5d0*F_l/F_sat_l))
       lim_r=max(0.0d0,min(1.0d0,1.0d0-0.5d0*F_r/F_sat_r))
       lim_d=max(0.0d0,min(1.0d0,1.0d0-0.5d0*F_d/F_sat_d))
       lim_u=max(0.0d0,min(1.0d0,1.0d0-0.5d0*F_u/F_sat_u))

       

       F_sat_l=F_sat_l*sign(1.0d0,-nxl*dTdx_l-nyl*dTdy_l)
       F_sat_r=F_sat_r*sign(1.0d0,-nxr*dTdx_r-nyr*dTdy_r)
       F_sat_d=F_sat_d*sign(1.0d0,-nxd*dTdx_d-nyd*dTdy_d)
       F_sat_u=F_sat_u*sign(1.0d0,-nxu*dTdx_u-nyu*dTdy_u)


       source(i,j,k)=source(i,j,k)&
            -dt*(((1.0d0-lim_r)*F_sat_r*nxr-(1.0d0-lim_l)*F_sat_l*nxl)/dxc(i) &
            +((1.0d0-lim_u)*F_sat_u*nyu-(1.0d0-lim_d)*F_sat_d*nyd)/dyc(j))
       

       bxxl=-eta*te_xm25*nxl*nxl/dxc(i)/dx(i-1)*lim_l
       bxxr=-eta*te_xp25*nxr*nxr/dxc(i)/dx(i)*lim_r
       byyd=-eta*te_ym25*nyd*nyd/dyc(j)/dy(j-1)*lim_d
       byyu=-eta*te_yp25*nyu*nyu/dyc(j)/dy(j)*lim_u

       bxyl=-eta/4.0d0/dxc(i)/dyc(j)*nxl*nyl*te_xm25*lim_l
       bxyr=-eta/4.0d0/dxc(i)/dyc(j)*nxr*nyr*te_xp25*lim_r
       byxd=-eta/4.0d0/dxc(i)/dyc(j)*nyd*nxd*te_ym25*lim_d
       byxu=-eta/4.0d0/dxc(i)/dyc(j)*nyu*nxu*te_yp25*lim_u

       matrix(i,j,k,2)=bxxl-(byxu-byxd)
       matrix(i,j,k,3)=bxxr+(byxu-byxd)
       matrix(i,j,k,4)=byyd-(bxyr-bxyl)
       matrix(i,j,k,5)=byyu+(bxyr-bxyl)

       matrix(i,j,k,6)=bxyl+byxd
       matrix(i,j,k,7)=-(bxyr+byxd)
       matrix(i,j,k,8)=-(bxyl+byxu)
       matrix(i,j,k,9)=bxyr+byxu
       matrix(i,j,k,1)=1.0-(matrix(i,j,k,2)+matrix(i,j,k,3)+matrix(i,j,k,4)+matrix(i,j,k,5))

    enddo;enddo;enddo
  
    !--------------------------------------------------------------------
  end subroutine set_matrix_source_MHD2D9_limiter

  
  subroutine mpiBND(var,n_margin)
    integer,intent(in)::n_margin
    double precision,intent(inout)::var(ix,jx,kx)
    double precision recvx(n_margin,jx,kx),sendx(n_margin,jx,kx)
    double precision recvy(ix,n_margin,kx),sendy(ix,n_margin,kx)
    double precision recvz(ix,jx,n_margin),sendz(ix,jx,n_margin)    
    integer nsend,nrecv,xcount,ycount,zcount
    xcount=jx*kx*n_margin
    ycount=ix*kx*n_margin
    zcount=ix*jx*n_margin

    !x-direction
    if(mpi_siz(1).gt.1) then
       !x-plus
       sendx=var(ix-margin(1)-n_margin+1:ix-margin(1),:,:)
       nsend=neighbor(2)
       nrecv=neighbor(1)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null

       call mpi_SendRecv(sendx,xcount,mpi_double_precision,nsend,0, &
            recvx,xcount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)
       if(neighbor(1).ne.-1) then
          var(margin(1)-n_margin+1:margin(1),:,:)=recvx
       endif
       !x-minus
       sendx=var(margin(1)+1:margin(1)+n_margin,:,:)
       nsend=neighbor(1)
       nrecv=neighbor(2)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendx,xcount,mpi_double_precision,nsend,0, &
            recvx,xcount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)
       if(neighbor(2).ne.-1) then
          var(ix-margin(1)+1:ix-margin(1)+n_margin,:,:)=recvx
       endif
    endif
    !y-direction
    if(mpi_siz(2).gt.1) then
       !=============================================================
       !y-plus
       sendy=var(:,jx-margin(2)-n_margin+1:jx-margin(2),:)      
       nsend=neighbor(4)
       nrecv=neighbor(3)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ycount,mpi_double_precision,nsend,0, &
            recvy,ycount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(3).ne.-1) then          
          var(:,margin(2)-n_margin+1:margin(2),:)=recvy
       endif

       !=============================================================
       !y-minus
       sendy=var(:,margin(2)+1:margin(2)+n_margin,:)
       nsend=neighbor(3)
       nrecv=neighbor(4)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ycount,mpi_double_precision,nsend,0, &
            recvy,ycount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(4).ne.-1) then
          var(:,jx-margin(2)+1:jx-margin(2)+n_margin,:)=recvy
       endif
    endif
    !for z-direction
    if(mpi_siz(3).gt.1) then
       !=============================================================
       !z-plus
       sendz=var(:,:,kx-margin(3)-n_margin+1:kx-margin(3))
       nsend=neighbor(6)
       nrecv=neighbor(5)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,zcount,mpi_double_precision,nsend,0, &
            recvz,zcount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(5).ne.-1) then
          var(:,:,1+margin(3)-n_margin:margin(3))=recvz(:,:,:)
       endif

       !=============================================================
       !right-boundary of z-direction
       sendz=var(:,:,margin(3)+1:margin(3)+n_margin)
       nsend=neighbor(5)
       nrecv=neighbor(6)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,zcount,mpi_double_precision,nsend,0, &
            recvz,zcount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)
       
       if(neighbor(6).ne.-1) then
          var(:,:,kx-margin(3)+1:kx-margin(3)+n_margin)=recvz
       endif
    endif
  end subroutine mpiBND

  !! A second-order accurate in time Super Timestepping method
  !! for (anisotropic) heat conduction (Ref: Meyer+2012, MNRS)
  subroutine HC_STS(U,nvar,dt,mhd,flag_sts)
    double precision,intent(inout) :: U(ix,jx,kx,nvar)
    double precision,intent(inout) :: dt
    integer, intent(in) :: nvar,mhd,flag_sts
    integer :: sub,j
    double precision :: tau
    double precision :: sol0(ix,jx,kx),sol1(ix,jx,kx),sol2(ix,jx,kx) &
         ,sol_new(ix,jx,kx)
    double precision, allocatable :: b_sts(:),mu_sts(:),mut_sts(:),nu_sts(:),gmt_sts(:)
    double precision :: bb(ix,jx,kx),vv(ix,jx,kx)
    double precision :: ro(ix,jx,kx),vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx),pr(ix,jx,kx) &
         ,bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx),te(ix,jx,kx)
    double precision :: jd,subd,cond_func0(ix,jx,kx)
    double precision :: coef_hc(ix,jx,kx,3,3),nbvec(ix,jx,kx,3)
    double precision :: nbx_int(ix,jx,kx,3),nby_int(ix,jx,kx,3),nbz_int(ix,jx,kx,3)
    
    select case(mhd)
    case(0) !! HD (isotropic conduction)
       call cq2pv_hd(ro,vx,vy,vz,pr,U)       
       te=gm*(pr/ro)
       bb = 0.d0
       vv = vx**2 + vy**2 + vz**2
       bx = 0.d0; by = 0.d0; bz = 0.d0
    case(1) !! MHD (anisotropic. Switch to isotropic cond when B is very small)
       call cq2pv_mhd(ro,vx,vy,vz,pr,bx,by,bz,U)
       te=gm*(pr/ro)    
       bb = bx**2 + by**2 + bz**2
       vv = vx**2 + vy**2 + vz**2
    end select

    !! calculate the coefficients that depend on the magnetic field    
    call cal_coef_hc(coef_hc,nbvec,nbx_int,nby_int,nbz_int,bx,by,bz,bb)
    
    !! timestep    
    tau = 0.5d0*dt    
    call cal_dt_cnd(dt_cnd,ro,te)
    sub = int(0.5d0 * (sqrt(9.d0+16.d0*tau/dt_cnd)-1.d0) +1.d0) !! +1 is necessary
    if(mod(sub,2)==0) sub = sub + 1
    !! avoid unphysical solutions resulting from
    !! the large difference in size between dt_cnd and dt_mhd
    if(flag_sts .eq. 0 .and. sub .ge. nsub_max) then
       sub = nsub_max
!       dt  = dble(nsub_max) * dt_cnd !! decrease dt_mhd
       tau = 0.25d0*dble(nsub_max*nsub_max+nsub_max-2) * dt_cnd !! Meyer+2012 eq(18)
       dt = 2.d0 * tau !! recalculate dt
    endif
    
    !! Preparation of the coefficients for STS
    allocate(b_sts(0:sub),mu_sts(0:sub),mut_sts(0:sub),nu_sts(0:sub) &
         ,gmt_sts(0:sub))
    b_sts = 0.d0; mu_sts = 0.d0; mut_sts = 0.d0; nu_sts = 0.d0; gmt_sts = 0.d0

    b_sts(0) = 1.d0/3.d0
    b_sts(1) = 1.d0/3.d0
    subd = dble(sub)
    do j = 2,sub
       jd = dble(j)
       b_sts(j) = (jd**2+jd-2.d0)/(2.d0*jd*(jd+1.d0))
    enddo
    mut_sts(1) = 4.d0 / (3.d0*(subd**2+subd-2.d0))
    do j = 2,sub
       jd = dble(j)
       mu_sts(j)  = (2.d0*jd-1.d0)*b_sts(j)/(jd*b_sts(j-1))
       nu_sts(j)  = -(jd-1.d0)*b_sts(j)/(jd*b_sts(j-2))
       mut_sts(j) = (4.d0*(2.d0*jd-1.d0)*b_sts(j)) &
            / (jd*(subd**2+subd-2.d0)*b_sts(j-1))
       gmt_sts(j) = -(1.d0-b_sts(j-1))*mut_sts(j)
    enddo
    
    !! Start subcycle (Total energy is conserved)
    sol0 = U(:,:,:,5)
    j = 1
    cond_func0 = cond_func(sol0,ro,vv,bb,coef_hc,nbvec,nbx_int,nby_int,nbz_int)
    sol1 = sol0 + mut_sts(j) * tau * cond_func0    
    call bnd_energy(sol1)
    
    
    j = 2
    sol2 = mu_sts(j)*sol1 + nu_sts(j)*sol0 &
         + (1.d0-mu_sts(j)-nu_sts(j))*sol0 &
         + mut_sts(j)*tau*cond_func(sol1,ro,vv,bb,coef_hc,nbvec,nbx_int,nby_int,nbz_int) &
         + gmt_sts(j)*tau*cond_func0
    call bnd_energy(sol2)

    do j=3,sub
       sol_new = mu_sts(j)*sol2 + nu_sts(j)*sol1 &
            + (1.d0-mu_sts(j)-nu_sts(j))*sol0 &
            + mut_sts(j)*tau*cond_func(sol2,ro,vv,bb,coef_hc,nbvec,nbx_int,nby_int,nbz_int) &
            + gmt_sts(j)*tau*cond_func0
       call bnd_energy(sol_new)

       sol1 = sol2
       sol2 = sol_new
    enddo

    U(:,:,:,5) = sol_new
    
    deallocate(b_sts,mu_sts,mut_sts,nu_sts,gmt_sts)

    if(flag_sts.eq.0) total_iter = sub
  end subroutine HC_STS



  !! Subcycling
  subroutine HC_SUBCYCLE(U,nvar,dt,mhd)
    double precision,intent(inout) :: U(ix,jx,kx,nvar)
    double precision,intent(inout) :: dt
    integer, intent(in) :: nvar,mhd
    integer :: sub,j
    double precision :: tau,dt_sub
    double precision :: sol0(ix,jx,kx),sol1(ix,jx,kx),sol2(ix,jx,kx) &
         ,sol_new(ix,jx,kx)
    double precision, allocatable :: b_sts(:)
    double precision :: bb(ix,jx,kx),vv(ix,jx,kx)
    double precision :: ro(ix,jx,kx),vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx),pr(ix,jx,kx) &
         ,bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx),te(ix,jx,kx)
    double precision :: dummy(ix,jx,kx,nvar)    
    double precision :: coef_hc(ix,jx,kx,3,3),nbvec(ix,jx,kx,3)
    double precision :: nbx_int(ix,jx,kx,3),nby_int(ix,jx,kx,3),nbz_int(ix,jx,kx,3)
    
    select case(mhd)
    case(0) !! HD (isotropic conduction)
       call cq2pv_hd(ro,vx,vy,vz,pr,U)       
       te=gm*(pr/ro)
       bb = 0.d0
       vv = vx**2 + vy**2 + vz**2
       bx = 0.d0; by = 0.d0; bz = 0.d0
    case(1) !! MHD (anisotropic. Switch to isotropic cond when B is very small)
       call cq2pv_mhd(ro,vx,vy,vz,pr,bx,by,bz,U)
       te=gm*(pr/ro)    
       bb = bx**2 + by**2 + bz**2
       vv = vx**2 + vy**2 + vz**2
    end select

    !! calculate the coefficients that depend on the magnetic field    
    call cal_coef_hc(coef_hc,nbvec,nbx_int,nby_int,nbz_int,bx,by,bz,bb)

    call cal_dt_cnd(dt_cnd,ro,te)
    sub = ceiling(dt/dt_cnd)
    !! avoid unphysical solutions resulting from
    !! the large difference in size between dt_cnd and dt_mhd
    if(sub .ge. nsub_max) then
       sub = nsub_max
       dt  = dble(nsub_max) * dt_cnd !! decrease dt_mhd
       tau = dt !! recalculate tau
    endif
    dt_sub = dt/dble(sub)

    !! Start subcycle (Total energy is conserved)
    sol_new = U(:,:,:,5)
    do j=1,sub
       sol_new = sol_new &
            + dt_sub * cond_func(sol_new,ro,vv,bb,coef_hc,nbvec,nbx_int,nby_int,nbz_int)
       call bnd_energy(sol_new)
    enddo

    U(:,:,:,5) = sol_new
    
    total_iter = sub       
  end subroutine HC_SUBCYCLE

  

  subroutine cal_dt_cnd(dt_cnd,ro,te)
    double precision, intent(in) :: ro(ix,jx,kx)
    double precision, intent(inout) :: te(ix,jx,kx)
    double precision, intent(out) :: dt_cnd
    double precision :: kc(ix,jx,kx),kappa0(ix,jx,kx)
    double precision :: gmin,dtq,c_p
    integer :: i,j,k
    
    dt_cnd = 999999.d0
    c_p = 1.d0 / (gm-1.d0)
    
    select case(hc_type)
    case(0) ! Constant-type
       kc(:,:,:) = kapper
    case(1) ! Spitzer-type
       call cal_kc_spitzer(kc,te)
    end select


    kappa0 = kc / c_p / max(ro,tiny)

    do k=margin(3)+1,kx-margin(3)
       do j=margin(2)+1,jx-margin(2)
          do i=margin(1)+1,ix-margin(1)
             gmin=min(dx(i),dy(j),dz(k))
             dtq = safety_cnd * 0.5d0 * gmin**2 / kappa0(i,j,k)
             dt_cnd = min(dt_cnd,dtq)
          enddo
       enddo
    enddo

    if(flag_mpi.eq.1) then
       dt_cnd=mpi_double_interface(dt_cnd,1)
    endif
    
  end subroutine cal_dt_cnd


  subroutine cal_kc_spitzer(kc1,te)
    double precision, intent(in) :: te(ix,jx,kx)
    double precision, intent(out) :: kc1(ix,jx,kx)

    !    kc1 = kapper * te**2.5
    kc1 = kapper * te*te*sqrt(te)    
    
  end subroutine cal_kc_spitzer


  !! Calculate -div F_cond [ Delta dot (kappa_para*Delta_para*T) ]
  function cond_func(sol,ro,vv,bb,coef_hc,nbvec,nbx_int,nby_int,nbz_int) result(res)
    double precision, intent(in) :: sol(ix,jx,kx) !! total energy
    double precision, intent(in) :: ro(ix,jx,kx),vv(ix,jx,kx),bb(ix,jx,kx) &
         ,coef_hc(ix,jx,kx,3,3),nbvec(ix,jx,kx,3)
    double precision, intent(in) :: nbx_int(ix,jx,kx,3),nby_int(ix,jx,kx,3),nbz_int(ix,jx,kx,3)
    double precision :: res(ix,jx,kx)    
    double precision :: te(ix,jx,kx),kc(ix,jx,kx)
    double precision :: flux_hc(ix,jx,kx,3)
    integer :: i,j,k

    te = gm*(gm-1.d0)/ro*(sol-0.5d0*ro*vv-0.5d0*bb)
    te = max(te,tiny)    

    select case(hc_type)
    case(0) ! Constant-type
       kc = kapper
    case(1) ! Spitzer-type
       call cal_kc_spitzer(kc,te)
    end select
    
    !! NB:
    !! x-component of flux_hc is defined at i+1/2,j,k
    !! y-component of flux_hc is defined at i,j+1/2,k
    !! z-component of flux_hc is defined at i,j,k+1/2
    call cal_hc_flux(flux_hc,ro,te,kc,coef_hc,nbvec,nbx_int,nby_int,nbz_int)

    do k=1,kx; do j=1,jx; do i=2,ix
       res(i,j,k) = - (flux_hc(i,j,k,1)-flux_hc(i-1,j,k,1))/dx(i)
    enddo;enddo;enddo
    if(flag_cyl.eq.1) then
       do k=1,kx; do j=1,jx; do i=2,ix
          res(i,j,k) = res(i,j,k) - 0.5d0*(flux_hc(i-1,j,k,1)+flux_hc(i,j,k,1))/x(i)
       enddo;enddo;enddo
    endif

    if (ndim .ge. 2) then
       do k=1,kx; do j=2,jx; do i=1,ix
          res(i,j,k) = res(i,j,k) - (flux_hc(i,j,k,2)-flux_hc(i,j-1,k,2))/dy(j)
       enddo;enddo;enddo
       if (ndim .ge. 3) then
          do k=2,kx; do j=1,jx; do i=1,ix
             res(i,j,k) = res(i,j,k) - (flux_hc(i,j,k,3)-flux_hc(i,j,k-1,3))/dz(k)
          enddo;enddo;enddo
       endif
    endif
    
  end function cond_func


  !! Currently second-order accurate in space 
  !! Isotropic conduction is used in the region where B is very small or zero
  !! Limiter proposed by Sharma & Hammett 2007 is adopted (minmod is used)
  !! [[ Set upper limit on the heat conduction flux (Ref: Cowie & McKee 1977) ]]
  subroutine cal_hc_flux(flux_hc_net,ro,te,kc,coef_hc,nbvec,nbx_int,nby_int,nbz_int)
    double precision, parameter :: bb_lim = 1d-14
    double precision, intent(inout) :: flux_hc_net(ix,jx,kx,3)
    double precision, intent(in) :: ro(ix,jx,kx),te(ix,jx,kx),kc(ix,jx,kx) &
         ,coef_hc(ix,jx,kx,3,3),nbvec(ix,jx,kx,3)
    double precision, intent(in) :: nbx_int(ix,jx,kx,3),nby_int(ix,jx,kx,3),nbz_int(ix,jx,kx,3)
    double precision :: flux_hc(ix,jx,kx,3)
    double precision :: dtedx(ix,jx,kx),dtedx_ys(ix,jx,kx),dtedx_zs(ix,jx,kx)
    double precision :: dtedy_xs(ix,jx,kx),dtedy(ix,jx,kx),dtedy_zs(ix,jx,kx)
    double precision :: dtedz_xs(ix,jx,kx),dtedz_ys(ix,jx,kx),dtedz(ix,jx,kx)    
    double precision :: kc_xs(ix,jx,kx),kc_ys(ix,jx,kx),kc_zs(ix,jx,kx)
    double precision :: kc_sat(ix,jx,kx)
    double precision :: kc_sat_xs(ix,jx,kx),kc_sat_ys(ix,jx,kx),kc_sat_zs(ix,jx,kx)
    double precision :: ro_xs(ix,jx,kx),ro_ys(ix,jx,kx),ro_zs(ix,jx,kx)
    double precision :: te_xs(ix,jx,kx),te_ys(ix,jx,kx),te_zs(ix,jx,kx)
    double precision :: flux_sat(ix,jx,kx,3),ratio_flux(ix,jx,kx)
    double precision :: sign_sat(ix,jx,kx),coef_sat(ix,jx,kx),g_lim0(ix,jx,kx)
    double precision :: temp1(ix,jx,kx),temp2(ix,jx,kx)
    integer :: i,j,k,n

!    kc_sat = 5.d0*ro*te*sqrt(te) !! for saturation flux (defined at i,j,k)

    !! ndim=1
    !! defined at i+1/2,j,k
    do k=1,kx; do j=1,jx; do i=1,ix-1
       dtedx(i,j,k) = (te(i+1,j,k)-te(i,j,k))/dx(i)
       kc_xs(i,j,k) = 2.d0*kc(i+1,j,k)*kc(i,j,k) &
            /max((kc(i+1,j,k)+kc(i,j,k)),tiny) !! harmonic mean
       ! ro_xs(i,j,k) = 2.d0*ro(i+1,j,k)*ro(i,j,k) &
       !      /max(ro(i+1,j,k)+ro(i,j,k),tiny) !! harmonic mean
       ! te_xs(i,j,k) = 2.d0*te(i+1,j,k)*te(i,j,k) &
       !      /max(te(i+1,j,k)+te(i,j,k),tiny) !! harmonic mean
       ! kc_sat_xs(i,j,k) = 2.d0*kc_sat(i+1,j,k)*kc_sat(i,j,k) &
       !      /max(kc_sat(i+1,j,k)+kc_sat(i,j,k),tiny) !! harmonic mean
    enddo; enddo; enddo

    flux_hc(:,:,:,1) = - kc_xs * coef_hc(:,:,:,1,1) * dtedx !! defined at i+1/2,j,k

    ! !! Preparation to set upper limit
    ! kc_sat_xs = 5.d0*ro_xs*te_xs*sqrt(te_xs)
    ! flux_sat(:,:,:,1) = -kc_sat_xs * nbvec(:,:,:,1) !! defined at i+1/2,j,k
    
    !! ndim=2 or 3
    if (ndim .ge. 2) then
       !! defined at i,j+1/2,k
       do k=1,kx; do j=1,jx-1; do i=1,ix
          dtedy(i,j,k) = (te(i,j+1,k)-te(i,j,k))/dy(j) 
          kc_ys(i,j,k) = 2.d0*kc(i,j+1,k)*kc(i,j,k) &
               /max((kc(i,j+1,k)+kc(i,j,k)),tiny) !! harmonic mean
          ! ro_ys(i,j,k) = 2.d0*ro(i,j+1,k)*ro(i,j,k) &
          !      /max(ro(i,j+1,k)+ro(i,j,k),tiny) !! harmonic mean
          ! te_ys(i,j,k) = 2.d0*te(i,j+1,k)*te(i,j,k) &
          !      /max(te(i,j+1,k)+te(i,j,k),tiny) !! harmonic mean
          ! ! kc_sat_ys(i,j,k) = 2.d0*kc_sat(i,j+1,k)*kc_sat(i,j,k) &
          ! !      /max(kc_sat(i,j+1,k)+kc_sat(i,j,k),tiny) !! harmonic mean
       enddo;enddo;enddo

       !! defined at i+1/2,j,k
       temp1 = 0.d0 ; temp2 = 0.d0
       temp1(:,2:jx,:)      = minmod_two(dtedy(:,1:jx-1,:),dtedy(:,2:jx,:),ix,jx-1,kx)
       temp2(1:ix-1,2:jx,:) = minmod_two(dtedy(2:ix,1:jx-1,:),dtedy(2:ix,2:jx,:),ix-1,jx-1,kx)
       dtedy_xs = minmod_two(temp1,temp2,ix,jx,kx)       
       ! do k=1,kx; do j=2,jx-1; do i=1,ix-1
       !    dtedy_xs(i,j,k) = 0.25d0*( dtedy(i,j,k)   + dtedy(i+1,j,k) &
       !                              +dtedy(i,j-1,k) + dtedy(i+1,j-1,k) )
       ! enddo;enddo;enddo

       !! defined at i,j+1/2,k
       temp1 = 0.d0 ; temp2 = 0.d0
       temp1(2:ix,:,:)      = minmod_two(dtedx(1:ix-1,:,:),dtedx(2:ix,:,:),ix-1,jx,kx)
       temp2(2:ix,1:jx-1,:) = minmod_two(dtedx(1:ix-1,2:jx,:),dtedx(2:ix,2:jx,:),ix-1,jx-1,kx)
       dtedx_ys = minmod_two(temp1,temp2,ix,jx,kx)
       ! do k=1,kx; do j=1,jx-1; do i=2,ix
       !    dtedx_ys(i,j,k) = 0.25d0*( dtedx(i-1,j+1,k) + dtedx(i,j+1,k) &
       !                              +dtedx(i-1,j,k)   + dtedx(i,j,k))
       ! enddo;enddo;enddo
       
       flux_hc(:,:,:,1) = flux_hc(:,:,:,1)     &
            - kc_xs * coef_hc(:,:,:,2,1) * dtedy_xs  !! defined at i+1/2,j,k
       flux_hc(:,:,:,2) = &
            - kc_ys * ( coef_hc(:,:,:,1,2) * dtedx_ys + coef_hc(:,:,:,2,2) * dtedy )   !! defined at i,j+1/2,k
            
       ! !! Preparation to set upper limit
       ! kc_sat_ys = 5.d0*ro_ys*te_ys*sqrt(te_ys)
       ! flux_sat(:,:,:,2) = -kc_sat_ys * nbvec(:,:,:,2)
       
       if (ndim .ge. 3) then
          !! defined at i,j,k+1/2
          do k=1,kx-1
             dtedz(:,:,k) = (te(:,:,k+1)-te(:,:,k))/dz(k) 
             kc_zs(:,:,k) = 2.d0*kc(:,:,k+1)*kc(:,:,k) &
                  /max((kc(:,:,k+1)+kc(:,:,k)),tiny) !! harmonic mean
             ! ro_zs(:,:,k) = 2.d0*ro(:,:,k+1)*ro(:,:,k) &
             !      /max((ro(:,:,k+1)+ro(:,:,k)),tiny) !! harmonic mean
             ! te_zs(:,:,k) = 2.d0*te(:,:,k+1)*te(:,:,k) &
             !      /max((te(:,:,k+1)+te(:,:,k)),tiny) !! harmonic mean
             ! ! kc_sat_zs(:,:,k) = 2.d0*kc_sat(:,:,k+1)*kc_sat(:,:,k) &
             ! !      /max((kc_sat(:,:,k+1)+kc_sat(:,:,k)),tiny) !! harmonic mean
          enddo

          !! defined at i+1/2,j,k
          temp1 = 0.d0 ; temp2 = 0.d0
          temp1(:,:,2:kx)      = minmod_two(dtedz(:,:,1:kx-1),dtedz(:,:,2:kx),ix,jx,kx-1)
          temp2(1:ix-1,:,2:kx) = minmod_two(dtedz(2:ix,:,1:kx-1),dtedz(2:ix,:,2:kx),ix-1,jx,kx-1)
          dtedz_xs = minmod_two(temp1,temp2,ix,jx,kx)
          ! do k=2,kx; do j=1,jx; do i=1,ix-1
          !    dtedz_xs(i,j,k) = 0.25d0*( dtedz(i+1,j,k-1) + dtedz(i+1,j,k) &
          !                              +dtedz(i,j,k-1)   + dtedz(i,j,k))
          ! enddo;enddo;enddo

          !! defined at i,j+1/2,k
          temp1 = 0.d0 ; temp2 = 0.d0
          temp1(:,:,2:kx)      = minmod_two(dtedz(:,:,1:kx-1),dtedz(:,:,2:kx),ix,jx,kx-1)
          temp2(:,1:jx-1,2:kx) = minmod_two(dtedz(:,2:jx,1:kx-1),dtedz(:,2:jx,2:kx),ix,jx-1,kx-1)
          dtedz_ys = minmod_two(temp1,temp2,ix,jx,kx)
          ! do k=2,kx; do j=1,jx; do i=1,ix-1
          !    dtedz_ys(i,j,k) = 0.25d0*( dtedz(i,j+1,k-1) + dtedz(i,j+1,k) &
          !                              +dtedz(i,j,k-1)   + dtedz(i,j,k))
          ! enddo;enddo;enddo

          !! defined at i,j,k+1/2
          temp1 = 0.d0 ; temp2 = 0.d0
          temp1(2:ix,:,:)      = minmod_two(dtedx(1:ix-1,:,:),dtedx(2:ix,:,:),ix-1,jx,kx)
          temp2(2:ix,:,1:kx-1) = minmod_two(dtedx(1:ix-1,:,2:kx),dtedx(2:ix,:,2:kx),ix-1,jx,kx-1)
          dtedx_zs = minmod_two(temp1,temp2,ix,jx,kx)
          ! do k=1,kx-1; do j=1,jx; do i=2,ix
          !    dtedx_zs(i,j,k) = 0.25d0*( dtedx(i-1,j,k)   + dtedx(i,j,k) &
          !                              +dtedx(i-1,j,k+1) + dtedx(i,j,k+1))
          ! enddo;enddo;enddo

          !! defined at i,j,k+1/2
          temp1 = 0.d0 ; temp2 = 0.d0
          temp1(:,2:jx,:)      = minmod_two(dtedy(:,1:jx-1,:),dtedy(:,2:jx,:),ix,jx-1,kx)
          temp2(:,2:jx,1:kx-1) = minmod_two(dtedy(:,1:jx-1,2:kx),dtedy(:,2:jx,2:kx),ix,jx-1,kx-1)
          dtedy_zs = minmod_two(temp1,temp2,ix,jx,kx)
          ! do k=1,kx-1; do j=2,jx; do i=1,ix
          !    dtedy_zs(i,j,k) = 0.25d0*( dtedy(i,j-1,k)   + dtedy(i,j,k) &
          !                              +dtedy(i,j-1,k+1) + dtedy(i,j,k+1))
          ! enddo;enddo;enddo
          
          
          flux_hc(:,:,:,1) = flux_hc(:,:,:,1) &
               - kc_xs * coef_hc(:,:,:,3,1) * dtedz_xs  !! defined at i+1/2,j,k
          flux_hc(:,:,:,2) = flux_hc(:,:,:,2) &
               - kc_ys * coef_hc(:,:,:,3,2) * dtedz_ys  !! defined at i,j+1/2,k
          flux_hc(:,:,:,3) = &
               - kc_zs * ( coef_hc(:,:,:,1,3) * dtedx_zs &
                         + coef_hc(:,:,:,2,3) * dtedy_zs &
                         + coef_hc(:,:,:,3,3) * dtedz ) !! defined at i,j,k+1/2

          ! !! Preparation to set upper limit
          ! kc_sat_zs = 5.d0*ro_zs*te_zs*sqrt(te_zs)
          ! flux_sat(:,:,:,3) = -kc_sat_zs * nbvec(:,:,:,3)
          
       endif
    endif

    ! !! set upper limit
    ! select case(ndim)
    ! case(1)
    !    sign_sat = sign(1.d0,nbx_int(:,:,:,1)*dtedx)
    !    flux_sat(:,:,:,1) = flux_sat(:,:,:,1) * sign_sat
    ! case(2)
    !    sign_sat = sign(1.d0,nbx_int(:,:,:,1)*dtedx &
    !                        +nby_int(:,:,:,1)*dtedy_xs)
    !    flux_sat(:,:,:,1) = flux_sat(:,:,:,1) * sign_sat

    !    sign_sat = sign(1.d0,nbx_int(:,:,:,2)*dtedx_ys &
    !                        +nby_int(:,:,:,2)*dtedy)
    !    flux_sat(:,:,:,2) = flux_sat(:,:,:,2) * sign_sat
    ! case(3)
    !    sign_sat = sign(1.d0,nbx_int(:,:,:,1)*dtedx &
    !                        +nby_int(:,:,:,1)*dtedy_xs &
    !                        +nbz_int(:,:,:,1)*dtedz_xs)          
    !    flux_sat(:,:,:,1) = flux_sat(:,:,:,1) * sign_sat

    !    sign_sat = sign(1.d0,nbx_int(:,:,:,2)*dtedx_ys &
    !                        +nby_int(:,:,:,2)*dtedy &
    !                        +nbz_int(:,:,:,2)*dtedz_ys)          
    !    flux_sat(:,:,:,2) = flux_sat(:,:,:,2) * sign_sat
       
    !    sign_sat = sign(1.d0,nbx_int(:,:,:,3)*dtedx_zs &
    !                        +nby_int(:,:,:,3)*dtedy_zs &
    !                        +nbz_int(:,:,:,3)*dtedz)          
    !    flux_sat(:,:,:,3) = flux_sat(:,:,:,3) * sign_sat
    ! end select

    ! !! NB:
    ! !! Fclas = | Fclas_vec*b_vec | = kc*(b_vec*nabla(T))
    ! !! Fsat  = | Fsat_vec*b_vec |  = 5*rho*Cs^3
    ! !! => Fclas/Fsat = kc*(b_vec*nabla(T)) / ( 5*rho*Cs^3 )
    ! !! The number of the non-zero terms in nabla(T) depends on ndim.    
    ! select case(ndim)
    ! case(1)
    !    !! x-component (defined at i+1/2,j,k)
    !    ratio_flux = abs(kc_xs * nbx_int(:,:,:,1)*dtedx) / kc_sat_xs
    !    g_lim0 = g_lim(ratio_flux) 
    !    flux_hc_net(:,:,:,1) = g_lim0 * flux_hc(:,:,:,1) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,1)
    ! case(2)
    !    !! x-component (defined at i+1/2,j,k)
    !    ratio_flux = abs(kc_xs * ( nbx_int(:,:,:,1)*dtedx &
    !                             + nby_int(:,:,:,1)*dtedy_xs )) &
    !                             / kc_sat_xs
    !    g_lim0 = g_lim(ratio_flux)
    !    flux_hc_net(:,:,:,1) = g_lim0 * flux_hc(:,:,:,1) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,1)
    !    !! y-component (defined at i,j+1/2,k)
    !    ratio_flux = abs(kc_ys * ( nbx_int(:,:,:,2)*dtedx_ys &
    !                             + nby_int(:,:,:,2)*dtedy )) &
    !                             / kc_sat_ys
    !    g_lim0 = g_lim(ratio_flux)
    !    flux_hc_net(:,:,:,2) = g_lim0 * flux_hc(:,:,:,2) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,2)       
    ! case(3)
    !    !! x-component (defined at i+1/2,j,k)
    !    ratio_flux = abs(kc_xs * ( nbx_int(:,:,:,1)*dtedx &
    !                             + nby_int(:,:,:,1)*dtedy_xs &
    !                             + nbz_int(:,:,:,1)*dtedz_xs )) &            
    !                             / kc_sat_xs
    !    g_lim0 = g_lim(ratio_flux)
    !    flux_hc_net(:,:,:,1) = g_lim0 * flux_hc(:,:,:,1) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,1)
    !    !! y-component (defined at i,j+1/2,k)
    !    ratio_flux = abs(kc_ys * ( nbx_int(:,:,:,2)*dtedx_ys &
    !                             + nby_int(:,:,:,2)*dtedy  &
    !                             + nbz_int(:,:,:,2)*dtedz_ys )) &            
    !                             / kc_sat_ys
    !    g_lim0 = g_lim(ratio_flux)
    !    flux_hc_net(:,:,:,2) = g_lim0 * flux_hc(:,:,:,2) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,2)       
    !    !! z-component (defined at i,j,k+1/2)
    !    ratio_flux = abs(kc_zs * ( nbx_int(:,:,:,3)*dtedx_zs &
    !                             + nby_int(:,:,:,3)*dtedy_zs  &
    !                             + nbz_int(:,:,:,3)*dtedz )) &            
    !                             / kc_sat_zs
    !    g_lim0 = g_lim(ratio_flux)
    !    flux_hc_net(:,:,:,3) = g_lim0 * flux_hc(:,:,:,3) &
    !         + (1.d0-g_lim0) * flux_sat(:,:,:,3)       
    ! end select

    flux_hc_net = flux_hc !! for test
    !    flux_hc_net = flux_sat !! for test    
    
  end subroutine cal_hc_flux
  

  function g_lim(ratio) result(res)
    double precision, intent(in) :: ratio(ix,jx,kx)
    double precision :: res(ix,jx,kx)

    res = max(0.d0,min(1.d0,1.d0-0.5d0*ratio))

  end function g_lim


  !! NB:
  !! coef_hc(:,:,:,1:3,1) : defined at i+1/2,j,k
  !! coef_hc(:,:,:,1:3,2) : defined at i,j+1/2,k
  !! coef_hc(:,:,:,1:3,3) : defined at i,j,k+1/2
  !! nbvec(:,:,:,1) : defined at i+1/2,j,k
  !! nbvec(:,:,:,2) : defined at i,j+1/2,k
  !! nbvec(:,:,:,3) : defined at i,j,k+1/2
  subroutine cal_coef_hc(coef_hc,nbvec,nbx_int,nby_int,nbz_int,bx,by,bz,bb)
    double precision, intent(out) :: coef_hc(ix,jx,kx,3,3),nbvec(ix,jx,kx,3)
    double precision, intent(out) :: nbx_int(ix,jx,kx,3),nby_int(ix,jx,kx,3),nbz_int(ix,jx,kx,3)
    double precision, intent(in) :: bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx),bb(ix,jx,kx)
    double precision :: bx_xs(ix,jx,kx),bx_ys(ix,jx,kx),bx_zs(ix,jx,kx)
    double precision :: by_xs(ix,jx,kx),by_ys(ix,jx,kx),by_zs(ix,jx,kx)
    double precision :: bz_xs(ix,jx,kx),bz_ys(ix,jx,kx),bz_zs(ix,jx,kx)
    double precision :: bb_xs(ix,jx,kx),bb_ys(ix,jx,kx),bb_zs(ix,jx,kx)
    double precision :: bbi_xs(ix,jx,kx),bbi_ys(ix,jx,kx),bbi_zs(ix,jx,kx)
    double precision :: rbbi_xs(ix,jx,kx),rbbi_ys(ix,jx,kx),rbbi_zs(ix,jx,kx)    

    double precision, parameter :: bb_lim = 1.d-10
    integer :: i,j,k
    
    !! initialization
    coef_hc = 0.d0 ; nbvec = 0.d0

    !! coefficients for x-component of the conduction flux
    do k=1,kx; do j=1,jx; do i=1,ix-1
       bx_xs(i,j,k) = 0.5d0*(bx(i,j,k)+bx(i+1,j,k))
       by_xs(i,j,k) = 0.5d0*(by(i,j,k)+by(i+1,j,k))
       bz_xs(i,j,k) = 0.5d0*(bz(i,j,k)+bz(i+1,j,k))       
    enddo;enddo;enddo
    bb_xs  = bx_xs*bx_xs+by_xs*by_xs+bz_xs*bz_xs
    bbi_xs = 1.d0/max(bb_xs,tiny)
    rbbi_xs = sqrt(bbi_xs)
    
    coef_hc(:,:,:,1,1) = bx_xs*bx_xs * bbi_xs
    where(bb_xs.lt.bb_lim) coef_hc(:,:,:,1,1) = 1.d0 !! for isotropic conduction
    nbvec(:,:,:,1) = bx_xs * rbbi_xs
    nbx_int(:,:,:,1) = bx_xs * rbbi_xs
    nby_int(:,:,:,1) = by_xs * rbbi_xs
    nbz_int(:,:,:,1) = bz_xs * rbbi_xs
    
    if(ndim.ge.2) then
       do k=1,kx; do j=1,jx-1; do i=1,ix
          bx_ys(i,j,k) = 0.5d0*(bx(i,j,k)+bx(i,j+1,k))
          by_ys(i,j,k) = 0.5d0*(by(i,j,k)+by(i,j+1,k))
          bz_ys(i,j,k) = 0.5d0*(bz(i,j,k)+bz(i,j+1,k))
       enddo;enddo;enddo
       bb_ys = bx_ys*bx_ys + by_ys*by_ys + bz_ys*bz_ys
       bbi_ys= 1.d0/max(bb_ys,tiny)
       rbbi_ys = sqrt(bbi_ys)
       
       coef_hc(:,:,:,2,1) = bx_xs*by_xs * bbi_xs
       coef_hc(:,:,:,1,2) = bx_ys*by_ys * bbi_ys
       coef_hc(:,:,:,2,2) = by_ys*by_ys * bbi_ys
       where(bb_xs.lt.bb_lim) !! for isotropic conduction
          coef_hc(:,:,:,2,1) = 0.d0
       end where
       where(bb_ys.lt.bb_lim) !! for isotropic conduction
          coef_hc(:,:,:,1,2) = 0.d0
          coef_hc(:,:,:,2,2) = 1.d0
       end where
       nbvec(:,:,:,2) = by_ys * rbbi_ys
       nbx_int(:,:,:,2) = bx_ys * rbbi_ys
       nby_int(:,:,:,2) = by_ys * rbbi_ys
       nbz_int(:,:,:,2) = bz_ys * rbbi_ys

       if(ndim.ge.3) then
          do k=1,kx-1
             bx_zs(:,:,k) = 0.5d0*(bx(:,:,k)+bx(:,:,k+1))
             by_zs(:,:,k) = 0.5d0*(by(:,:,k)+by(:,:,k+1))
             bz_zs(:,:,k) = 0.5d0*(bz(:,:,k)+bz(:,:,k+1))
          enddo
          bb_zs = bx_zs*bx_zs + by_zs*by_zs + bz_zs*bz_zs
          bbi_zs= 1.d0/max(bb_zs,tiny)
          rbbi_zs = sqrt(bbi_zs)
          
          coef_hc(:,:,:,3,1) = bx_xs*bz_xs * bbi_xs
          coef_hc(:,:,:,3,2) = by_ys*bz_ys * bbi_ys
          coef_hc(:,:,:,1,3) = bx_zs*bz_zs * bbi_zs
          coef_hc(:,:,:,2,3) = by_zs*bz_zs * bbi_zs
          coef_hc(:,:,:,3,3) = bz_zs*bz_zs * bbi_zs
          where(bb_xs.lt.bb_lim)
             coef_hc(:,:,:,3,1) = 0.d0
          end where
          where(bb_ys.lt.bb_lim)
             coef_hc(:,:,:,3,2) = 0.d0
          end where
          where(bb_zs.lt.bb_lim)
             coef_hc(:,:,:,1,3) = 0.d0
             coef_hc(:,:,:,2,3) = 0.d0
             coef_hc(:,:,:,3,3) = 1.d0
          end where
          nbvec(:,:,:,3) = bz_zs * rbbi_zs
          nbx_int(:,:,:,3) = bx_zs * rbbi_zs
          nby_int(:,:,:,3) = by_zs * rbbi_zs
          nbz_int(:,:,:,3) = bz_zs * rbbi_zs

       endif
    endif
          
  end subroutine cal_coef_hc

  function minmod_two(a,b,ix1,jx1,kx1) result(res)
    double precision,intent(in) :: a(ix1,jx1,kx1),b(ix1,jx1,kx1)
    integer,intent(in) :: ix1,jx1,kx1
    double precision :: res(ix1,jx1,kx1)

    res = sign(1.d0,a)*max(0.d0,min(abs(a),sign(1.d0,a)*b))

  end function minmod_two

  
end module HC_rot
