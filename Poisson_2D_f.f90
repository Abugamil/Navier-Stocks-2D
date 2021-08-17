subroutine Poisson_2D

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson
use variables
use constants
integer nx,ny
parameter(nx=MaxN-1, ny=MaxN-1)


integer i,j, stat
integer ipar(128)
real ax, bx, ay, by, lx, ly, hx, hy, xi, yi, cx, cy, c1
real spar(13*nx/2+7)
! Note that proper packing of data in right-hand side array f is
! automatically provided by the following declarations of the array
real f(nx+1,ny+1)
real bd_ax(ny+1), bd_bx(ny+1), bd_ay(nx+1), bd_by(nx+1)
real q
type(DFTI_DESCRIPTOR), pointer :: xhandle
character(4) BCtype

! Defining the rectangular domain 0<x<1, 0<y<1 for 2D Poisson Solver
ax=0.0E0
bx=1.0E0
ay=0.0E0
by=1.0E0

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the way to use Helmholtz Solver to solve Poisson problem!
!*******************************************************************************
q=0.0E0

! Computing the mesh size hx in x-direction
lx=bx-ax
hx=lx/nx
! Computing the mesh size hy in y-direction
ly=by-ay
hy=ly/ny

! Filling in the values of the TRUE solution u(x,y)=sin(2*pi*x)*sin(2*pi*y)+1
! in the mesh points into the array u
! Filling in the right-hand side f(x,y)=(8*pi*pi+q)*sin(2*pi*x)*sin(2*pi*y)+q
! in the mesh points into the array f.
! We choose the right-hand side to correspond to the TRUE solution of
! Poisson equation.
! Here we are using the mesh sizes hx and hy computed before to compute
! the coordinates (xi,yi) of the mesh points
do j=1,ny
   do i=1,nx
      f(i,j)=((Uon(i+1,j)-Uon(i,j))/dx+(Von(i,j+1)-Von(i,j))/dy)/dt
   enddo
enddo

! Setting the type of the boundary conditions on each side of the rectangular domain:
! On the boundary laying on the line x=0(=ax) Dirichlet boundary condition
! will be used
! On the boundary laying on the line x=1(=bx) Dirichlet boundary condition
! will be used
! On the boundary laying on the line y=0(=ay) Neumann boundary condition will be used
! On the boundary laying on the line y=1(=by) Neumann boundary condition will be used
BCtype = 'NNNN'

! Setting the values of the boundary function G(x,y) that is equal to
! the TRUE solution in the mesh points laying on Dirichlet boundaries
do i = 1,ny
   bd_ax(i) = 0.0
   bd_bx(i) = 0.0
enddo
! Setting the values of the boundary function g(x,y) that is equal to
! the normal derivative of the TRUE solution in the mesh points laying on
! Neumann boundaries
do i = 1,nx
   bd_ay(i) = 0.0
   bd_by(i) = 0.0
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for 2D Poisson Solver
call s_init_Helmholtz_2D(ax, bx, ay, by, nx, ny, BCtype, q, ipar, spar, stat)
call s_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, stat)
call s_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, stat)
call free_Helmholtz_2D(xhandle, ipar, stat)
call mkl_freebuffers
do j=1,ny
   do i=1,nx
      p(i,j)=f(i,j)
   enddo
enddo

end subroutine Poisson_2D
