module bmp_rot
  implicit none
  integer,save::width,height,bit_count,infoheader_length
  integer,save::compression,n_color,file_size
  integer,parameter::unit=11,fortran_header=4
  integer,allocatable,save::image(:,:,:) 
 
contains
  subroutine read_BMP(filename)
    character*100,intent(in)::filename
    compression=0
    open(unit,file=trim(filename),access="direct",recl=1)
    call read_fileheader
    if(infoheader_length.eq.12) then
       call read_infoheader12
    else if(infoheader_length.eq.40) then
       call read_infoheader40
    else
       print *,"NOT PROPER HEADER LENGTH"
    endif    

    if(bit_count.eq.24) then
       
    else
       print *,"COMPRESSION BMP IS NOT ALLOWED"
    endif
    close(unit)
    n_color=3    
    allocate(image(n_color,width,height))    
    file_size=14+infoheader_length+file_size
    call read_image(filename)    
!    open(11,file="tmp.dat",form="unformatted",status="replace")
!    write(11)image
!    close(11)


  end subroutine read_BMP

  subroutine get_BMP(out_image)
    double precision,intent(out)::out_image(n_color,width,height)
    out_image=image
  end subroutine get_BMP

  subroutine read_image(filename)
    character*100,intent(in)::filename
    integer n,i,j,start,count
    integer*1 buf(file_size)
    integer rec_size
    start=14+infoheader_length+1
!    rec_size=start+file_size-1
    open(unit,file=trim(filename),access="direct",recl=file_size)
    read(unit,rec=1)buf
    close(unit)
!    count=1
    image=reshape(buf(start:),(/n_color,width,height/))
    image=(mod(image+256,256))
!    print *,size(buf(start:)),size(image)
!    do j=1,height
!       do i=1,width
!         do n=1,n_color
!             call read2buf(start,height*width*n_color,buf)
!             image(n,i,j)=get_integer(4,buf(count:count+3))
!             image(n,i,j)=buf(count)
!             count=count+1
!          enddo
!       enddo
!    enddo
  end subroutine read_image


  subroutine read_fileheader
    integer,parameter::header_size=18
    integer*1 buf(header_size)    
    call read2buf(1,header_size,buf)
    infoheader_length=get_integer(4,buf(15:18))    
  end subroutine read_fileheader

  subroutine read_infoheader12
    integer,parameter::header_size=12
    integer*1 buf(header_size)    
    call read2buf(15,header_size,buf)
    width=get_integer(2,buf(5:6))    
    height=get_integer(2,buf(7:8))   
    bit_count=get_integer(2,buf(11:12)) 
  end subroutine read_infoheader12

  subroutine read_infoheader40
    integer,parameter::header_size=40
    integer*1 buf(header_size)    
    
    call read2buf(15,header_size,buf)
    width=get_integer(4,buf(5:8))    
    height=get_integer(4,buf(9:12))    
    bit_count=get_integer(2,buf(15:16))
    compression=get_integer(4,buf(17:20))
    file_size=get_integer(4,buf(21:24))
    if(file_size.eq.0) then
       file_size=width*height*3
    endif
!    print *,"WIDTH,HEIGHT,BIT_COUNT,COMPRESSION,FILE_SIZE"
!    print *,WIDTH,HEIGHT,BIT_COUNT,COMPRESSION,FILE_SIZE
  end subroutine read_infoheader40

  subroutine read2buf(start,size,buf)
    integer,intent(in)::start,size
    integer*1,intent(out)::buf(size)
    integer i
    do i=start,start+size-1
       read(unit,rec=i)buf(i-start+1)
    enddo
  end subroutine read2buf

  !get integer from bits (little endian)
  function get_integer(N,bits)
    integer,intent(in)::N
    integer*1,intent(in):: bits(N)
    integer i
    integer get_integer
    get_integer=0
    do i=1,N
       get_integer=get_integer+(mod(bits(i)+256,256))*256**(i-1)
    enddo
  end function get_integer
  
end module bmp_rot
