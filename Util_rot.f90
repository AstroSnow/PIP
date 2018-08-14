module Util_rot
  implicit none

contains  
  subroutine get_word(A,key,ind_s)
    character*200,intent(in)::A
    character*15,intent(out)::key
    integer ind_s,ind_e
    ind_s=index(A,':')
    ind_e=index(A,';')
    key=trim(adjustl(A(ind_s+1:ind_e-1)))
  end subroutine get_word

  subroutine get_value_integer(line,key,output)
    character*200,intent(in)::line
    character*15,intent(in)::key
    integer,intent(out)::output
    integer ind_s
    ind_s=index(line,":")
    if(key.eq.trim(adjustl(line(1:ind_s-1)))) then
       read(line(ind_s+1:),*)output
    else
       output=0
    endif    
  end subroutine get_value_integer


  subroutine get_hd_flux(u,f,gm)
    double precision,intent(in)::u(5),gm
    double precision,intent(out)::f(5)
    double precision pr
    pr=(gm-1.0d0)*(u(5)-0.5d0*(u(2)**2+u(3)**2+u(4)**2)/u(1))
    f(1)=u(2) 
    f(2)=pr +u(2)**2/u(1)
    f(3)=u(2)*u(3)/u(1)
    f(4)=u(2)*u(4)/u(1)
    f(5)=u(2)*(u(5)+pr)/u(1)
  end subroutine get_hd_flux
  subroutine get_mhd_flux(u,f,gm)
    double precision,intent(in)::u(8),gm
    double precision,intent(out)::f(8)
    double precision pr
    pr=(gm-1.0d0)*(u(5) &
         -0.5d0*(u(2)**2+u(3)**2+u(4)**2)/u(1)  &
         -0.5d0*(u(6)**2+u(7)**2+u(8)**2))

    f(1)=u(2)
    f(2)=u(2)**2/u(1)+pr-u(6)**2
    f(3)=u(2)*u(3)/u(1)-u(6)*u(7)
    f(4)=u(2)*u(4)/u(1)-u(6)*u(8)
  !  f(5)=v_l(1)*(e_l+pt_l)-bn*(v_l(1)*bn_l+v_l(2)*bt1_l+v_l(3)*bt2_l)
    f(5)=(u(2)*(u(5)+pr+0.5d0*(u(6)**2+u(7)**2+u(8)**2)) -&
         u(6)*(u(2)*u(6)+u(3)*u(7)+u(4)*u(8)))/u(1)
    f(6)=0.0d0
    f(7)=(u(2)*u(7)-u(3)*u(6))/u(1)
    f(8)=(u(2)*u(8)-u(4)*u(6))/u(1)
  end subroutine get_mhd_flux

  function minmod_func(U1,U2,U3,ix,jx,kx,nvar)
    integer,intent(in)::ix,jx,kx,nvar
    double precision::minmod_func(ix,jx,kx,nvar)
    double precision,intent(in)::U1(ix,jx,kx,nvar)
    double precision,intent(in)::U2(ix,jx,kx,nvar)
    double precision,intent(in)::U3(ix,jx,kx,nvar)
    double precision ::X(ix,jx,kx,nvar),Y(ix,jx,kx,nvar),SX(ix,jx,kx,nvar)
!    double precision ::Z(ix,jx,kx,nvar)
!    integer i,j,k
    
    X=U3-U2 ; Y=U2-U1
    SX=sign(1.0d0,X)
    minmod_func=SX*max(0.0d0,min(abs(X),SX*Y))
!    X=2.0d0*(U3-U2);Y=2.0d0*(U2-U1);Z=0.5d0*(U3-U1)
!    SX=sign(1.0d0,X)
!    minmod_func=SX*max(0.0d0,min(abs(X),SX*Y,SX*Z))
    return    
  end function minmod_func


  subroutine get_divergence(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U,D)
    integer,intent(in)::ix,jx,kx,ndim,s_order
    double precision,intent(in)::U(ix,jx,kx,3)
    double precision,intent(in)::dxc(ix),dyc(jx),dzc(kx)
    double precision,intent(out)::D(ix,jx,kx)
    integer i,j,k
    D=0.0d0
    if(s_order.eq.4) then              
       if(ndim.eq.1) then 
          do i=3,ix-2
             D(i,:,:)=(-U(i+2,:,:,1) &
                  +8.0d0*(U(i+1,:,:,1)-U(i-1,:,:,1))&
                  +U(i-2,:,:,1))/(12.0d0*dxc(i))
          enddo
       elseif(ndim.eq.2)then
          do j=3,jx-2
             do i=3,ix-2
                D(i,j,:)=(-U(i+2,j,:,1) &
                  +8.0d0*(U(i+1,j,:,1)-U(i-1,j,:,1))&
                  +U(i-2,j,:,1))/(12.0d0*dxc(i)) &
                  +(-U(i,j+2,:,2) &
                  +8.0d0*(U(i,j+1,:,2)-U(i,j-1,:,2))&
                  +U(i,j-2,:,2))/(12.0d0*dyc(j))
             enddo
          enddo
       elseif(ndim.eq.3)then
          do k=3,kx-2
             do j=3,jx-2
                do i=3,ix-2
                   D(i,j,k)=(-U(i+2,j,k,1) &
                        +8.0d0*(U(i+1,j,k,1)-U(i-1,j,k,1))&
                        +U(i-2,j,k,1))/(12.0d0*dxc(i)) &
                        +(-U(i,j+2,k,2) &
                        +8.0d0*(U(i,j+1,k,2)-U(i,j-1,k,2))&
                        +U(i,j-2,k,2))/(12.0d0*dyc(j)) &
                        +(-U(i,j,k+2,3) &
                        +8.0d0*(U(i,j,k+1,3)-U(i,j,k-1,3))&
                        +U(i,j,k-2,3))/(12.0d0*dzc(k))
                enddo
             enddo
          enddo
       endif
    else
       if(ndim.eq.1) then 
          do i=2,ix-1
             D(i,:,:)=(U(i+1,:,:,1)-U(i-1,:,:,1))/(2.0d0*dxc(i))
          enddo
       elseif(ndim.eq.2)then
          do j=3,jx-2
             do i=3,ix-2
                D(i,j,:)=(U(i+1,j,:,1)-U(i-1,j,:,1))/(2.0d0*dxc(i)) &
                  +(U(i,j+1,:,2)-U(i,j-1,:,2))/(2.0d0*dyc(j))
             enddo
          enddo
       elseif(ndim.eq.3)then
          do k=3,kx-2
             do j=3,jx-2
                do i=3,ix-2
                   D(i,j,k)=(U(i+1,j,k,1)-U(i-1,j,k,1))/(2.0d0*dxc(i)) &
                        +(U(i,j+1,k,2)-U(i,j-1,k,2))/(2.0d0*dyc(j)) &
                        +(U(i,j,k+1,3)-U(i,j,k-1,3))/(2.0d0*dzc(k))
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine get_divergence


  subroutine get_rotation(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U,R)
    integer,intent(in)::ix,jx,kx,ndim,s_order
    double precision,intent(in)::U(ix,jx,kx,3)
    double precision,intent(in)::dxc(ix),dyc(jx),dzc(kx)
    double precision,intent(out)::R(ix,jx,kx,3)
    integer i,j,k
    R(:,:,:,:)=0.0d0
    if(s_order.eq.4) then              
       if(ndim.eq.1) then 
          do i=3,ix-2
             R(i,:,:,2)=-(-U(i+2,:,:,3)+8.0d0*U(i+1,:,:,3) &
                  -8.0d0*U(i-1,:,:,3)+U(i-2,:,:,3)) &           
                  /(12.0d0*dxc(i))
             R(i,:,:,3)=(-U(i+2,:,:,2)+8.0d0*U(i+1,:,:,2) &
                  -8.0d0*U(i-1,:,:,2)+U(i-2,:,:,2)) &           
                  /(12.0d0*dxc(i))
          enddo
       else if (ndim.eq.2) then
          do j=3,jx-2
             do i=3,ix-2
                R(i,j,:,1)=(-U(i,j+2,:,3)+8.0d0*U(i,j+1,:,3) &
                     -8.0d0*U(i,j-1,:,3)+U(i,j-2,:,3))/(12.0d0*dyc(j))
                R(i,j,:,2)=-(-U(i+2,j,:,3)+8.0d0*U(i+1,j,:,3) &
                     -8.0d0*U(i-1,j,:,3)+U(i-2,j,:,3))/(12.0d0*dxc(i))
                R(i,j,:,3)=( &
                     (-U(i+2,j,:,2)+8.0d0*U(i+1,j,:,2) &
                     -8.0d0*U(i-1,j,:,2)+U(i-2,j,:,2)) &
                     /(12.0d0*dxc(i)) &
                     -(-U(i,j+2,:,1)+8.0d0*U(i,j+1,:,1) &
                     -8.0d0*U(i,j-1,:,1)+U(i,j-2,:,1))/(12.0d0*dyc(j)))
             enddo
          enddo
       else if (ndim.eq.3) then
          do k=3,kx-2
             do j=3,jx-2
                do i=3,ix-2
                   R(i,j,k,1)=( &
                        (-U(i,j+2,k,3)+8.0d0*U(i,j+1,k,3) &
                        -8.0d0*U(i,j-1,k,3)+U(i,j-2,k,3)) &
                        /(12.0d0*dyc(j)) &
                        -(-U(i,j,k+2,2)+8.0d0*U(i,j,k+1,2) &
                        -8.0d0*U(i,j,k-1,2)+U(i,j,k-2,2))/(12.0d0*dzc(k)))
                   R(i,j,k,2)=( &
                        (-U(i,j,k+2,1)+8.0d0*U(i,j,k+1,1) &
                        -8.0d0*U(i,j,k-1,1)+U(i,j,k-2,1)) &
                        /(12.0d0*dzc(k)) &
                        -(-U(i+2,j,k,3)+8.0d0*U(i+1,j,k,3) &
                        -8.0d0*U(i-1,j,k,3)+U(i-2,j,k,3))/(12.0d0*dxc(i)))
                   R(i,j,k,3)=( &
                        (-U(i+2,j,k,2)+8.0d0*U(i+1,j,k,2) &
                        -8.0d0*U(i-1,j,k,2)+U(i-2,j,k,2)) &
                        /(12.0d0*dxc(i)) &
                        -(-U(i,j+2,k,1)+8.0d0*U(i,j+1,k,1) &
                        -8.0d0*U(i,j-1,k,1)+U(i,j-2,k,1))/(12.0d0*dyc(j)))
                enddo
             enddo
          enddo
       endif
    else
       if(ndim.eq.1) then 
          do i=2,ix-1
             R(i,:,:,2)=-(U(i+1,:,:,3)-U(i-1,:,:,3))/(2.0d0*dxc(i))
             R(i,:,:,3)=(U(i+1,:,:,2)-U(i-1,:,:,2))/(2.0d0*dxc(i))
          enddo
       else if (ndim.eq.2) then
          do j=2,jx-1
             do i=2,ix-1
                R(i,j,:,1)=(U(i,j+1,:,3)-U(i,j-1,:,3))/(2.0d0*dyc(j))
                R(i,j,:,2)=-(U(i+1,j,:,3)-U(i-1,j,:,3))/(2.0d0*dxc(i))     
                R(i,j,:,3)=((U(i+1,j,:,2)-U(i-1,j,:,2))/(2.0d0*dxc(i)) &      
                     -(U(i,j+1,:,1)-U(i,j-1,:,1))/(2.0d0*dyc(j)))
             enddo
          enddo
       else if (ndim.eq.3) then
          do k=2,kx-1
             do j=2,jx-1
                do i=2,ix-1
                   R(i,j,k,1)=((U(i,j+1,k,3)-U(i,j-1,k,3))/(2.0d0*dyc(j)) &
                        -(U(i,j,k+1,2)-U(i,j,k-1,2))/(2.0d0*dzc(k)))
                   R(i,j,k,2)=((U(i,j,k+1,1)-U(i,j,k-1,1))/(2.0d0*dzc(k)) &
                        -(U(i+1,j,k,3)-U(i-1,j,k,3))/(2.0d0*dxc(i)))
                   R(i,j,k,3)=((U(i+1,j,k,2)-U(i-1,j,k,2))/(2.0d0*dxc(i)) &
                        -(U(i,j+1,k,1)-U(i,j-1,k,1))/(2.0d0*dyc(j)))
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine get_rotation


end module Util_rot
