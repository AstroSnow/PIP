module matrix_rot
  implicit none
contains
  subroutine inverse_tridiagonal(A,x,b,ix)
    ! solve linear simultaneous equations 
    
    !    A(:,1) :left comonent of tridiagonal matrix
    !    A(:,2) :diagonal comonent of tridiagonal matrix
    !    A(:,3) :right comonent of tridiagonal matrix
    
    integer,intent(in)::ix
    double precision,intent(in)::A(ix,3),b(ix)
    double precision,intent(inout)::x(ix)
    double precision d(ix),bb(ix),fact
    !    double precision c(ix),d(ix),e(ix),fact
    integer i
    !    c(:)=A(:,1);d(:)=A(:,2);e(:)=A(:,3)    
    d(:)=A(:,2);bb(:)=b
    do i=2,ix
       fact=-A(i,1)/d(i-1)
       d(i)=d(i)+fact*A(i-1,3)       
       bb(i)=bb(i)+fact*bb(i-1)
    enddo
    x(ix)=bb(ix)/d(ix)
    do i=ix-1,1,-1
       x(i)=(bb(i)-A(i,3)*x(i+1))/d(i)
    enddo    
  end subroutine inverse_tridiagonal
  
end module matrix_rot
