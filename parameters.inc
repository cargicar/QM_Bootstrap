! s= matrix size; g= coupling;
! xlower=lower bound x^2; xupper=upper bound x^2;
! elower=lower bound enery; eupper=upper bound energy;
!     xsize = grid size in x^2; esize = grid size in Energy;
integer, parameter :: dp = selected_real_kind(30, 307)
integer, parameter :: k=7, xsize = 3000, esize = 3000
real(dp), parameter :: w2=1, g=-0.07, xlower = 0.0, xupper = 1,  elower =0.0,  eupper = 2
