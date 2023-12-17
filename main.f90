program main
   implicit none 
   include 'parameters.inc'
  ! integer, parameter :: dp = real64
   integer :: i, j, p, q
   real(dp) :: xs, start, finish, tmp, xss, ess, xinter, einter
   real(dp), dimension (:,:), allocatable :: hankel
   real(dp), dimension(:,:),  allocatable :: L
  ! open(newunit=io, file="log.txt", status="replace", action="write")
   call cpu_time(start)
   xinter=xupper-xlower
   einter=eupper-elower
   allocate ( hankel(k+1,k+1) )
   allocate ( L(k+1,k+1) )
   ! ! write(6,'(3f10.5)') (L(i,:), i=1,n)
   write(1,*) k, xsize, g
   do i = 1, xsize
      xss=xlower+i*xinter/xsize
      do j= 1, esize
         ess=elower+j*einter/esize
          call hankel_matrix(k, ess, g, w2, xss, hankel)
          call cholesky(k+1, hankel, L)
         ! print*, "for xss = ",xss," and  ess =  ",ess," L is ", L
          if (.not.ANY( isnan(L) )) then
             write(2,*) ess, xss            
          end if 
       end do
   end do    

!   close(io)
   ! if (isnan()) stop '"x" is a NaN'
     
   call cpu_time(finish)
   print '("Time = ",f6.3," seconds.")',finish-start
 stop

end program main

recursive function xs(s, E, ge, w22, x22) result(res)
  implicit none
  include 'parameters.inc'
  integer::s
  real(dp) ::  res, E, ge, x22, w22
  if (MOD(s,2)==1) then
     res=0;
  else if (s < 0) then
     res =0;
  else if(s==0) then
     res=1
  else if(s==2) then
     res = x22
  else
     res=(4*(s-3)*E*xs(s-4,E,ge, w22, x22)+(s-3)*(s-4)*(s-5)*xs(s-6,E,ge, w22, x22)-4*(s-2)*w22*xs(s-2,E,ge, w22, x22))/(4*ge*(s-1))
   !   write(3,*), "x^",s,"=",res
     return
  end if
end function 

subroutine hankel_matrix(s, E, ge, w22, x2, hankel)
  include 'parameters.inc'
  integer :: i, j
  integer, intent(in) :: s
  real(dp) :: xs, start, finish, tmp
  real(dp), intent(in) :: E, ge, x2, w22
  real(dp), dimension (s+1,s+1), intent(out) :: hankel
 
 ! allocate ( hankel(s+1,s+1) )
  do i = 1, s+1           
     do j = 1, s+1                
        hankel(i,j) = xs(i+j-2,E,ge, w2, x2)               
      end do      
   end do

end subroutine hankel_matrix

subroutine cholesky(n, A, L)
   implicit none
   include 'parameters.inc'
   integer,    intent(in)    :: n
   real(kind=dp),     intent(in)    :: A(n,n)
   real(kind=dp),    intent(out)   :: L(n,n)
   integer    :: i,j
   do j = 1, n
       L(j,j) = sqrt( A(j,j) - dot_product(L(j,1:j-1),L(j,1:j-1)) )  
        do i = j+1, n
           L(i,j)  = ( A(i,j) - dot_product(L(i,1:j-1),L(j,1:j-1)) ) / L(j,j)
        end do
     end do
  end subroutine cholesky
