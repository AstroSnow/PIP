subroutine set_fold_grid(ix,margin,x,dx,dxc,s_order)
  integer,intent(in)::ix,margin,s_order
  double precision,intent(inout)::x(ix),dx(ix),dxc(ix)
  integer ix0,factor,ix1,ori_x
  double precision dx0,L_x

  factor=2
  L_x=2.0d0

  ix0=ix-2*margin
  ix1=ix0/(2*(1+factor))*factor

  dx0=L_x/4.0d0/ix1

  ori_x=ix/2+1
  x(ori_x)=dx0*0.5d0
  x(ori_x-1)=-dx0*0.5d0
  do i=1,ix1-1
     x(i+ori_x)=x(i+ori_x-1)+dx0
     x(ori_x-i-1)=x(ori_x-i)-dx0
  enddo
  x(ori_x+ix1)=x(ori_x+ix1-1)+0.5*(1.0+factor)*dx0
  x(ori_x-ix1-1)=x(ori_x-ix1)-0.5*(1.0+factor)*dx0
  
  do i=ix1+1,ix/2-1
     x(ori_x+i)=x(ori_x+i-1)+factor*dx0
     x(ori_x-i-1)=x(ori_x-i)-factor*dx0
  enddo

  
  dx(1:ix-1)=x(2:ix)-x(1:ix-1)
  dx(ix)=dx(ix-1)
  
  if(s_order.eq.4) then 
     do i=3,ix-1
        dxc(i)=0.25d0*(dx(i+1)+dx(i)+dx(i-1)+dx(i-2))
     enddo
     dxc(1:2)=dxc(3)
     dxc(ix)=dxc(ix-1)
  else
     do i=2,ix
        dxc(i)=0.5d0*(dx(i)+dx(i-1))
     enddo
     dxc(1)=dxc(2)
  endif

end subroutine set_fold_grid

