!-----------------------------------------------!
!     A  df/dx ->  Thomas algorithm             !
!-----------------------------------------------!
function addx(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=0.5/dx
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=1.0/dx
 else
     tempout=0.0
 endif
 addx=tempout
return
end function addx

!-----------------------------------------------!
!     B  df/dx ->  Thomas algorithm             !
!-----------------------------------------------!
function bddx(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=0.0
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=-1.0/dx
 elseif (Cscheme.eq.BackwardDifference) then
     tempout=1.0/dx
 else
     tempout=0.0
 endif
 bddx=tempout
return
end function bddx


!-----------------------------------------------!
!     C  df/dx ->  Thomas algorithm             !
!-----------------------------------------------!
function cddx(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=-0.5/dx
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=0.0
 elseif (Cscheme.eq.BackwardDifference) then
     tempout=-1.0/dx
 else
     tempout=0.0
 endif
 cddx=tempout
return
end function cddx

!-----------------------------------------------!
!     A  df/dy ->  Thomas algorithm             !
!-----------------------------------------------!
function addy(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=0.5/dy
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=1.0/dy
 else
     tempout=0.0
 endif
 addy=tempout
return
end function addy

!-----------------------------------------------!
!     B  df/dy ->  Thomas algorithm             !
!-----------------------------------------------!
function bddy(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=0.0
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=-1.0/dy
 elseif (Cscheme.eq.BackwardDifference) then
     tempout=1.0/dy
 else
     tempout=0.0
 endif
 bddy=tempout
return
end function bddy


!-----------------------------------------------!
!     C  df/dy ->  Thomas algorithm             !
!-----------------------------------------------!
function cddy(i,j)
use variables
use constants
real :: tempout
 if (Cscheme.eq.CentralDifference) then
     tempout=-0.5/dy
 elseif (Cscheme.eq.ForwardDifference) then
     tempout=0.0
 elseif (Cscheme.eq.BackwardDifference) then
     tempout=-1.0/dy
 else
     tempout=0.0
 endif
 cddy=tempout
return
end function cddy
!-----------------------------------------------!
!     A  d2f/dx2 ->  Thomas algorithm           !
!-----------------------------------------------!
function ad2dx2(i,j)
use variables
use constants
real :: tempout
 tempout=1.0/(dx*dx)
 ad2dx2=tempout
return
end function ad2dx2

!-----------------------------------------------!
!     B  d2f/dx2 ->  Thomas algorithm           !
!-----------------------------------------------!
function bd2dx2(i,j)
use variables
use constants
real :: tempout
 tempout=-2.0/(dx*dx)
 bd2dx2=tempout
return
end function bd2dx2

!-----------------------------------------------!
!    C  d2f/dx2 ->  Thomas algorithm            !
!-----------------------------------------------!
function cd2dx2(i,j)
use variables
use constants
real :: tempout
 tempout=1.0/(dx*dx)
 cd2dx2=tempout
return
end function cd2dx2

!-----------------------------------------------!
!     A  d2f/dy2 ->  Thomas algorithm           !
!-----------------------------------------------!
function ad2dy2(i,j)
use variables
use constants
real :: tempout
 tempout=1.0/(dy*dy)
 ad2dy2=tempout
return
end function ad2dy2

!-----------------------------------------------!
!     B  d2f/dy2 ->  Thomas algorithm           !
!-----------------------------------------------!
function bd2dy2(i,j)
use variables
use constants
real :: tempout
 tempout=-2.0/(dy*dy)
 bd2dy2=tempout
return
end function bd2dy2

!-----------------------------------------------!
!    C  d2f/dy2 ->  Thomas algorithm            !
!-----------------------------------------------!
function cd2dy2(i,j)
use variables
use constants
real :: tempout
 tempout=1.0/(dy*dy)
 cd2dy2=tempout
return
end function cd2dy2

!-----------------------------------------------!
!    du/dx in i,j                               !
!-----------------------------------------------!
function dudx(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   integer :: sloi
   real :: tempout=0.0
   if (sloi.eq.0) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(u(i+1,j)-u(i-1,j))/(2.0*dx)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(u(i+1,j)-u(i,j))/dx
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(u(i,j)-u(i-1,j))/dx
      else
         tempout=0.0
      endif
   elseif (sloi.eq.1) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(Uo12(i+1,j)-Uo12(i-1,j))/(2.0*dx)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(Uo12(i+1,j)-Uo12(i,j))/dx
      elseif (Cscheme.eq.BackwardDifference) then
          tempout=(Uo12(i,j)-Uo12(i-1,j))/dx
      else
          tempout=0.0
      endif
   endif

   dudx=tempout
   return
end function dudx

!-----------------------------------------------!
!    du/dy in i,j                               !
!-----------------------------------------------!
function dudy(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   integer :: sloi
   real :: tempout=0.0
   if (sloi.eq.0) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(u(i,j+1)-u(i,j-1))/(2.0*dy)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(u(i,j+1)-u(i,j))/dy
      elseif (Cscheme.eq.BackwardDifference) then
          tempout=(u(i,j)-u(i,j-1))/dy
      else
          tempout=0.0
      endif
   elseif (sloi.eq.1) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(Uo12(i,j+1)-Uo12(i,j-1))/(2.0*dy)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(Uo12(i,j+1)-Uo12(i,j))/dy
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(Uo12(i,j)-Uo12(i,j-1))/dy
      else
         tempout=0.0
      endif
   endif

   dudy=tempout
   return
end function dudy

!-----------------------------------------------!
!    dv/dx in i,j                               !
!-----------------------------------------------!
function dvdx(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   integer :: sloi
   real :: tempout=0.0
   if (sloi.eq.0) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(v(i+1,j)-v(i-1,j))/(2.0*dx)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(v(i+1,j)-v(i,j))/dx
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(v(i,j)-v(i-1,j))/dx
      else
         tempout=0.0
      endif
   elseif (sloi.eq.1) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(Vo12(i+1,j)-Vo12(i-1,j))/(2.0*dx)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(Vo12(i+1,j)-Vo12(i,j))/dx
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(Vo12(i,j)-Vo12(i-1,j))/dx
      else
         tempout=0.0
      endif
   endif

   dvdx=tempout
   return
end function dvdx

!-----------------------------------------------!
!    dv/dy in i,j                               !
!-----------------------------------------------!
function dvdy(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   integer :: sloi
   real :: tempout=0.0
   if (sloi.eq.0) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(v(i,j+1)-v(i,j-1))/(2.0*dy)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(v(i,j+1)-v(i,j))/dy
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(v(i,j)-v(i,j-1))/dy
      else
         tempout=0.0
      endif
   elseif (sloi.eq.1) then
      if (Cscheme.eq.CentralDifference) then
          tempout=(Vo12(i,j+1)-Vo12(i,j-1))/(2.0*dy)
      elseif (Cscheme.eq.ForwardDifference) then
          tempout=(Vo12(i,j+1)-Vo12(i,j))/dy
      elseif (Cscheme.eq.BackwardDifference) then
         tempout=(Vo12(i,j)-Vo12(i,j-1))/dy
      else
         tempout=0.0
      endif
   endif

   dvdy=tempout
   return
end function dvdy

!-----------------------------------------------!
!    d2u/dx2 in i,j                             !
!-----------------------------------------------!
function d2udx2(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   integer :: sloi
   real :: tempout=0.0
   if (sloi.eq.0) then
       tempout = (u(i+1,j)-2.0*u(i,j)+u(i-1,j))/(dx*dx)
   elseif (sloi.eq.1) then
       tempout = (Uo12(i+1,j)-2.0*Uo12(i,j)+Uo12(i-1,j))/(dx*dx)
   endif    
   d2udx2 = tempout
   return
end function d2udx2

!-----------------------------------------------!
!    d2v/dx2 in i,j                             !
!-----------------------------------------------!
function d2vdx2(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   real :: tempout=0.0
   integer :: sloi
   if (sloi.eq.0) then
       tempout = (v(i+1,j)-2.0*v(i,j)+v(i-1,j))/(dx*dx)
   elseif (sloi.eq.1) then
       tempout = (Vo12(i+1,j)-2.0*Vo12(i,j)+Vo12(i-1,j))/(dx*dx)
   endif    
   d2vdx2 = tempout
   return
end function d2vdx2

!-----------------------------------------------!
!    d2u/dy2 in i,j                             !
!-----------------------------------------------!
function d2udy2(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   real :: tempout
   integer :: sloi
   if (sloi.eq.0) then
       tempout = (u(i,j+1)-2.0*u(i,j)+u(i,j-1))/(dy*dy)
   elseif (sloi.eq.1) then
       tempout = (Uo12(i,j+1)-2.0*Uo12(i,j)+Uo12(i,j-1))/(dy*dy)
   endif    
   d2udy2 = tempout
   return
end function d2udy2

!-----------------------------------------------!
!    d2v/dy2 in i,j                             !
!-----------------------------------------------!
function d2vdy2(i,j,sloi) 
   use variables
   use constants
   integer :: i
   integer :: j
   real :: tempout
   integer :: sloi
   if (sloi.eq.0) then
       tempout = (v(i,j+1)-2.0*v(i,j)+v(i,j-1))/(dy*dy)
   elseif (sloi.eq.1) then
       tempout = (Vo12(i,j+1)-2.0*Vo12(i,j)+Vo12(i,j-1))/(dy*dy)
   endif    
   d2vdy2 = tempout
   return
end function d2vdy2

function S1(Ax)
   real :: Ax
   if (Ax.ge.0.0) then
       S1=0.0
   else
       S1=Ax
   endif
   return 
end function S1
function S2(Ax)
   real :: Ax
   if (Ax.ge.0.0) then
       S2=Ax
   else
       S2=0.0
   endif
   return 
end function S2
function S3(Bx)
   real :: Bx
   if (Bx.ge.0.0) then
       S3=0.0
   else
       S3=Bx
   endif
   return 
end function S3
function S4(Bx)
   real :: Bx
   if (Bx.ge.0.0) then
       S4=Bx
   else
       S4=0.0
   endif
   return 
end function S4