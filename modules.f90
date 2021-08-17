module constants
   integer, parameter :: MaxN = 180
   integer, parameter :: np = MaxN-1
   integer, parameter :: n = np*np
   integer, parameter :: nnz = n*3-2*np
   
   integer, parameter :: MaxTime=90000
   integer, parameter :: PrintEvery=100
   integer, parameter :: CentralDifference = 1
   integer, parameter :: ForwardDifference = 2
   integer, parameter :: BackwardDifference = 3
   integer, parameter :: CrossWind = 2

   double precision, parameter :: Cs = 0.1
   
   integer, parameter :: VarU = 1
   integer, parameter :: VarV = 2
   integer, parameter :: VarP = 4

   double precision, parameter :: dx = 1.0/MaxN
   double precision, parameter :: dy = 1.0/MaxN
   double precision, parameter :: dz = 1.0/MaxN
   double precision, parameter :: dt = dx-dx*dx*dx*dx
   
 
   double precision, parameter :: Re = 1000.0
end module constants
module variables
   integer :: Iter
   integer :: Cscheme = 1
   double precision, allocatable :: aUx1(:),bUx1(:),aUx2(:),bUx2(:)
   double precision, allocatable :: aUy1(:),bUy1(:),aUy2(:),bUy2(:)
   double precision, allocatable :: aVx1(:),bVx1(:),aVx2(:),bVx2(:)
   double precision, allocatable :: aVy1(:),bVy1(:),aVy2(:),bVy2(:)
   
   double precision, allocatable :: U(:,:),V(:,:),P(:,:)
   double precision, allocatable :: Uon(:,:),Von(:,:)
   double precision, allocatable :: Uo12(:,:),Vo12(:,:)
end module variables
