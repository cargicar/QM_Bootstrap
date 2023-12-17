program main
   implicit none 
   include 'parameters.inc'
   integer :: i, j, p, q
   complex :: xs, xc
   real :: start, finish, tmp, xss, ess, xinter, einter
   complex, dimension (:,:), allocatable :: hankel
   complex, dimension(:,:),  allocatable :: L
   complex, dimension(:),  allocatable :: diagL
   
  ! open(newunit=io, file="log.txt", status="replace", action="write")
   call cpu_time(start)
   xinter=xupper-xlower
   einter=eupper-elower
   allocate ( hankel(s+1,s+1) )
   allocate ( L(s+1,s+1) )
   allocate ( diagL(s+1) )
   ! ! write(6,'(3f10.5)') (L(i,:), i=1,n)
   write(1,*) s, xsize, g
   do i = 1, xsize
      xss=xlower+i*xinter/xsize
      xc= complex(xss, 0)
      do j= 1, esize
         ess=elower+j*einter/esize
          call hankel_matrix(s, ess, g, w2, xc , hankel)
          call cholesky(s+1, hankel, L)
          !print*, "for xss = ",xss," and  ess =  ",ess," L is ", size(L, dim=1)
          ploop: do p= 1,  size(L, dim=1)
             diagL(p)=L(p,p)
          end do ploop
          if (.not.any(imagpart(diagL).ne.0)) then
               write(2,*) ess, xss
          end if
          ! if (ALL( (imagpart(L))==0 )) then
          !        write(2,*) ess, xss            
          ! end if
          !  end do qloop
       end do
   end do    

!   close(io)
   ! if (isnan()) stop '"x" is a NaN'
     
   call cpu_time(finish)
   print '("Time = ",f6.3," seconds.")',finish-start
 stop

end program main

recursive function xs(s, E, g, w2, x2) result(res)
!  Equation 6 Hartnoll.  Input: x2 = x**2  Returns x**s
implicit none
integer::s
complex :: x2, res
real ::  E, g, w2
if (MOD(s,2)==1) then
   res=0;
else if (s < 0) then
   res =0;
else if(s==0) then
   res=1
else if(s==2) then
   res = x2
else
   res = (4*(s-3)*E*xs(s-4,E,g, w2, x2 )+(s-3)*(s-4)*(s-5)*xs(s-6,E,g, w2, x2)-4*(s-2)*w2*xs(s-2,E,g, w2, x2))/(4*g*(s-1))
return
end if
end function 

! recursive function xs(s, E, g, x2) result(res)
! implicit none
! integer::s
! complex :: x2, res
! real ::  E, g 
! if (MOD(s,2)==1) then
!    res=0;
! else if (s < 0) then
!    res =0;
! else if(s==0) then
!    res=1
! else if(s==2) then
!    res = x2
! else
!    res = (4*(s-3)*E*xs(s-4,E,g,x2)+(s-3)*(s-4)*(s-5)*xs(s-6,E,g,x2)-4*(s-2)*xs(s-2,E,g,x2))/(4*g*(s-1))
! !     write(3,*) "x^",s,"=",res
! return
! end if
! end function 

subroutine hankel_matrix(s, E, g, w2, x2, hankel)
  integer :: i, j
  integer, intent(in) :: s
  complex :: xs, xsc
  real :: start, finish, tmp
  real, intent(in) :: E, g, w2
  complex, intent(in) :: x2
  complex, dimension (s+1,s+1), intent(out) :: hankel
 
 ! allocate ( hankel(s+1,s+1) )
  do i = 1, s+1           
     do j = 1, s+1
       ! xsc=complex(xs(i+j-2,E,g,x2),0)
        hankel(i,j) = xs(i+j-2,E,g, w2, x2)               
      end do      
   end do

end subroutine hankel_matrix


subroutine cholesky(n, A, L)
   implicit none
   integer,    intent(in)    :: n
   complex,     intent(in)    :: A(n,n)
   complex,    intent(out)   :: L(n,n)
   integer    :: i,j
   do i = 1, size(A,1)
    L(i,i) = sqrt(A(i,i) - dot_product(L(1:i-1,i), L(1:i-1,i)))
    L(i,i+1:) = (A(i,i+1:) - matmul(conjg(L(1:i-1,i)), L(1:i-1,i+1:))) / L(i,i)
   end do
   ! do j = 1, n
   !     L(j,j) = sqrt( A(j,j) - dot_product(L(j,1:j-1),L(j,1:j-1)) )  
   !      do i = j+1, n
   !         L(i,j)  = ( A(i,j) - dot_product(L(i,1:j-1),L(j,1:j-1)) ) / L(j,j)
   !      end do
   !   end do
  end subroutine cholesky


! subroutine cholesky(n, A, L)
!    implicit none
!    integer,    intent(in)    :: n
!    complex,     intent(in)    :: A(n,n)
!    complex,    intent(out)   :: L(n,n)
!    integer    :: i,j
!    do j = 1, n
!        L(j,j) = sqrt( A(j,j) - dot_product(L(j,1:j-1),L(j,1:j-1)) )  
!         do i = j+1, n
!            L(i,j)  = ( A(i,j) - dot_product(L(i,1:j-1),L(j,1:j-1)) ) / L(j,j)
!         end do
!      end do
!   end subroutine cholesky
