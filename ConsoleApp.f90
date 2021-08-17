program prog
   use variables
   use constants

   call InitiliseVariables
   call InitFields
  
   
   
   call UVSolver  
   
   call DestroyVariables 
    
end program prog
!-------------------------------------!
! Setting boundary condition for UV   !
!-------------------------------------!
subroutine InitFields
   use variables
   use constants
   
   do i=1,MaxN
      aUx1(i)=0.0
      bUx1(i)=0.0
      aUx2(i)=0.0
      bUx2(i)=0.0
      
      aUy1(i)=0.0
      bUy1(i)=0.0
      aUy2(i)=0.0
      bUy2(i)=1.0
      
      aVx1(i)=0.0
      bVx1(i)=0.0
      aVx2(i)=0.0
      bVx2(i)=0.0
      
      aVy1(i)=0.0
      bVy1(i)=0.0
      aVy2(i)=0.0
      bVy2(i)=0.0
      
   enddo
   
   
   
   do i=0,MaxN+1
       do j=0,MaxN+1
           u(i,j)=0.0
           v(i,j)=0.0
           p(i,j)=0.0
           Uo12(i,j)=0.0
           Uon(i,j)=0.0
           Vo12(i,j)=0.0
           Von(i,j)=0.0
       enddo
   enddo
 !  forall (j=1:MaxN) u(0,j)=1.0
end subroutine InitFields
!-------------------------------------!
! Run solver for velocity             !
! Rouch 140 stranica                  !
!-------------------------------------!
subroutine UVSolver
   use variables
   use constants
   real, dimension(0:MaxN) :: A 
   real, dimension(0:MaxN) :: B
   real, dimension(0:MaxN) :: C
   real, dimension(0:MaxN) :: F
   real, dimension(1:MaxN) :: X
   integer :: i,j,k
   integer :: printtime
   do iter=1,MaxTime
	!***************  U po x       
	   do j=1,MaxN
		 do i=1,MaxN
            A(i) = -(-u(i,j)*addx(i,j) + ad2dx2(i,j)/Re)
            B(i) = -(-u(i,j)*bddx(i,j) + bd2dx2(i,j)/Re) + 2.0/dt
            C(i) = -(-u(i,j)*cddx(i,j) + cd2dx2(i,j)/Re)
            F(i) =(-v(i,j)*dudy(i,j,0) + d2udy2(i,j,0)/Re) + 2.0*u(i,j)/dt
         enddo
         A(0)=aUx1(j)
         B(0)=bUx1(j)
         C(0)=aUx2(j)
         F(0)=bUx2(j)       
		 call sweepF(A,B,C,F,X)
         forall (i=1:MaxN)	Uo12(i,j)=X(i)  
       enddo
!       forall (i=1:MaxN, j=1:MaxN) u(i,j)=Uo12(i,j)
!       call WriteFile(VarU,Iter+MaxTime)
	!***************  U po y       
	   do i=1,MaxN
		 do j=1,MaxN
            A(j) = -(-v(i,j)*addy(i,j) + ad2dy2(i,j)/Re)
            B(j) = -(-v(i,j)*bddy(i,j) + bd2dy2(i,j)/Re) + 2.0/dt
            C(j) = -(-v(i,j)*cddy(i,j) + cd2dy2(i,j)/Re)
            F(j) =(-u(i,j)*dudx(i,j,1) + d2udx2(i,j,1)/Re) + 2.0*Uo12(i,j)/dt
         enddo
         A(0)=aUy1(i)
         B(0)=bUy1(i)
         C(0)=aUy2(i)
         F(0)=bUy2(i)
    	 call sweepF(A,B,C,F,X)
         forall (j=1:MaxN) Uon(i,j)=X(j)
       enddo
!	forall (i=1:MaxN, j=1:MaxN) u(i,j)=Uo23(i,j)
	!***************  V po x       
	   do j=1,MaxN
		 do i=1,MaxN
            A(i) = -(-u(i,j)*addx(i,j) + ad2dx2(i,j)/Re)
            B(i) = -(-u(i,j)*bddx(i,j) + bd2dx2(i,j)/Re) + 2.0/dt
            C(i) = -(-u(i,j)*cddx(i,j) + cd2dx2(i,j)/Re)
            F(i) =(-v(i,j)*dvdy(i,j,0) + d2vdy2(i,j,0)/Re) + 2.0*v(i,j)/dt
         enddo
         A(0)=aVx1(j)
         B(0)=bVx1(j)
         C(0)=aVx2(j)
         F(0)=bVx2(j)
		 call sweepF(A,B,C,F,X)
         forall (i=1:MaxN)	Vo12(i,j)=X(i)  
       enddo
	forall (i=1:MaxN, j=1:MaxN) v(i,j)=Vo12(i,j)
	!**************  V po y
	   do i=1,MaxN
		 do j=1,MaxN
            A(j) = -(-v(i,j)*addy(i,j) + ad2dy2(i,j)/Re)
            B(j) = -(-v(i,j)*bddy(i,j) + bd2dy2(i,j)/Re) + 2.0/dt
            C(j) = -(-v(i,j)*cddy(i,j) + cd2dy2(i,j)/Re)
            F(j) =(-u(i,j)*dvdx(i,j,1) + d2vdx2(i,j,1)/Re) + 2.0*Vo12(i,j)/dt
         enddo
         A(0)=aVy1(i)
         B(0)=bVy1(i)
         C(0)=aVy2(i)
         F(0)=bVy2(i)
		 call sweepF(A,B,C,F,X)
         forall (j=1:MaxN) Von(i,j)=X(j)
       enddo
!    forall (i=1:MaxN, j=1:MaxN) v(i,j)=Vo23(i,j)
!    forall (i=1:MaxN, j=1:MaxN) u(i,j)=Uo23(i,j)
    call Poisson_2D
    do j=1,MaxN
        do i=2,MaxN-1
            u(i,j)=Uon(i,j)-dt*(p(i,j)-p(i-1,j))/dx
            
        enddo
    enddo
    do i=1,MaxN
        do j=2,MaxN-1
            v(i,j)=Von(i,j)-dt*(p(i,j)-p(i,j-1))/dy
        enddo
    enddo
           
!    call WriteFile(VarU,Iter)
!    call WriteFile(VarV,Iter)
     if (Mod(Iter,PrintEvery).eq.0) then
         WRITE(*,*) 'Time: ', Iter
         printtime=Iter/PrintEvery
         call WriteIter(printtime)
     endif
   enddo
end subroutine UVSolver
!-------------------------------------!
! Thomas method                       !
!-------------------------------------!
subroutine sweepF(A,B,C,F,X) 
   use variables
   use constants
   real :: A(0:MaxN),B(0:MaxN),C(0:MaxN),F(0:MaxN)
   real :: X(1:MaxN)
   integer :: I
   real :: alfa(MaxN+1),beta(MaxN+1)
   real :: temp
   alfa(MaxN)=C(0)
   beta(MaxN)=F(0)
   do I=MaxN-1,1,-1
      temp =(B(I)+A(I)*alfa(I+1))
      
      alfa(I)=-C(I)/temp
      beta(I)=(F(I)-A(I)*beta(I+1))/temp
   enddo
     X(1)=(alfa(1)*B(0)+beta(1))/(1.0-alfa(1)*A(0))
   do I=1,  MaxN-1
      X(I+1)=alfa(I+1)*X(I)+beta(I+1)
   enddo
end subroutine sweepF

subroutine sweepB(A,B,C,F,X) 
   use variables
   use constants
   real :: A(0:MaxN),B(0:MaxN),C(0:MaxN),F(0:MaxN)
   real :: X(1:MaxN)
   integer :: I,M
   real :: alfa(MaxN+1),beta(MaxN+1)
   alfa(1)=A(0)
   beta(1)=B(0)
   do I=1,MaxN-1
      alfa(I+1)=-C(I)/(B(I)+A(I)*alfa(I) )
      beta(I+1)=(F(I)-A(I)*beta(I))/(B(I)+A(I)*alfa(I))
   enddo
     X(MaxN)=(F(0)+C(0)*beta(MaxN))/(1-C(0)*alfa(MaxN))
   do I=MaxN-1,1,-1
      X(I)=alfa(I+1)*X(I+1)+beta(I+1)
   enddo
end subroutine sweepB

!   elseif (ProblemDim.eq.3D) then
!
!      dudx(i,j,k) = (u(i+1,j,k)-u(i-1,j,k))/(2.0*dx)
!      dvdx(i,j,k) = (v(i+1,j,k)-v(i-1,j,k))/(2.0*dx)
!      dwdx(i,j,k) = (w(i+1,j,k)-w(i-1,j,k))/(2.0*dx)
!
!      dudy(i,j,k) = (u(i,j+1,k)-u(i,j-1,k))/(2.0*dy)
!      dvdy(i,j,k) = (v(i,j+1,k)-v(i,j-1,k))/(2.0*dy)
!      dwdy(i,j,k) = (w(i,j+1,k)-w(i,j-1,k))/(2.0*dy)
!
!      dudz(i,j,k) = (u(i,j,k+1)-u(i,j,k-1))/(2.0*dz)
!      dvdz(i,j,k) = (v(i,j,k+1)-v(i,j,k-1))/(2.0*dz)
!      dwdz(i,j,k) = (w(i,j,k+1)-w(i,j,k-1))/(2.0*dz)
!
!      d2udx2(i,j,k) = (u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))/(2.0*dx*dx)
!      d2vdx2(i,j,k) = (v(i+1,j,k)-2.0*v(i,j,k)+v(i-1,j,k))/(2.0*dx*dx)
!      d2wdx2(i,j,k) = (w(i+1,j,k)-2.0*w(i,j,k)+w(i-1,j,k))/(2.0*dx*dx)
!
!      d2udy2(i,j,k) = (u(i,j+1,k)-2.0*u(i,j,k)+u(i,j-1,k))/(2.0*dy*dy)
!      d2vdy2(i,j,k) = (v(i,j+1,k)-2.0*v(i,j,k)+v(i,j-1,k))/(2.0*dy*dy)
!      d2wdy2(i,j,k) = (w(i,j+1,k)-2.0*w(i,j,k)+w(i,j-1,k))/(2.0*dy*dy)
!
!      d2udz2(i,j,k) = (u(i,j,k+1)-2.0*u(i,j,k)+u(i,j,k-1))/(2.0*dz*dz)
!      d2vdz2(i,j,k) = (v(i,j,k+1)-2.0*v(i,j,k)+v(i,j,k-1))/(2.0*dz*dz)
!      d2wdz2(i,j,k) = (w(i,j,k+1)-2.0*w(i,j,k)+w(i,j,k-1))/(2.0*dz*dz)
!   endif
!--------------------------------------------!
! Koefficienti of raznosti protiv potoka 2   !
!--------------------------------------------!



subroutine WriteFile(TVar,Num)
   use variables
   use constants
   Character(LEN=15) :: FName
   Integer :: TVar
   Integer :: Num
   If (TVar.eq.VarU) then
    Write(FName,"(A,I5.5,A)") 'U',Num,'.dat'
   else
    Write(FName,"(A,I5.5,A)") 'V',Num,'.dat'
   endif
    
   open(TVar,file=trim(adjustl(FName)))
   write(TVar,*) 'TITLE = "Example: Simple 2D-Plane Data"'
   write(TVar,*) 'VARIABLES = "X", "Y", "Var"'
   write(TVar,*) 'ZONE I=',MaxN,'  J=',MaxN,' F=POINT'
   do i=1,MaxN
     do j=1,MaxN
       if (TVar.eq.VarU) then
          write(TVar,*) i,j,u(i,j)
       elseif (TVar.eq.VarV) then
          write(TVar,*) i,j,v(i,j)
!       else
!          write(1,*) i,j,p(i,j)
       endif
     enddo
   enddo
   close(TVar)  
end subroutine WriteFile
subroutine WriteIter(Num)
   use variables
   use constants
   Character(LEN=15) :: FName
   Integer :: TVar=10
   Integer :: Num
   Write(FName,"(A,I5.5,A)") 'Out',Num,'.dat'
    
   open(TVar,file=trim(adjustl(FName)))
   write(TVar,*) 'TITLE = "Example: Simple 2D-Plane Data"'
   write(TVar,*) 'VARIABLES = "X", "Y", "U", "V"'
   write(TVar,*) 'ZONE I=',MaxN,'  J=',MaxN,' F=POINT'
   do i=1,MaxN
     do j=1,MaxN
          write(TVar,*) i,j,u(i,j),v(i,j)
     enddo
   enddo
   close(TVar)  
end subroutine WriteIter

!-------------------------------------!
! Set Variables                       !
!-------------------------------------!
subroutine InitiliseVariables
   use variables
   use constants
   
   allocate(P(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate P ****" 
   
   allocate(aUx1(1:MaxN),STAT=istat)
   allocate(bUx1(1:MaxN),STAT=istat)
   allocate(aUx2(1:MaxN),STAT=istat)
   allocate(bUx2(1:MaxN),STAT=istat)
   allocate(aUy1(1:MaxN),STAT=istat)
   allocate(bUy1(1:MaxN),STAT=istat)
   allocate(aUy2(1:MaxN),STAT=istat)
   allocate(bUy2(1:MaxN),STAT=istat)
   
   allocate(aVx1(1:MaxN),STAT=istat)
   allocate(bVx1(1:MaxN),STAT=istat)
   allocate(aVx2(1:MaxN),STAT=istat)
   allocate(bVx2(1:MaxN),STAT=istat)
   allocate(aVy1(1:MaxN),STAT=istat)
   allocate(bVy1(1:MaxN),STAT=istat)
   allocate(aVy2(1:MaxN),STAT=istat)
   allocate(bVy2(1:MaxN),STAT=istat)
   
   allocate(U(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate U ****"
   allocate(V(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate V ****"
   allocate(Uon(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate Uon1 ****"
   allocate(Von(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate Von1 ****"

   allocate(Uo12(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate Uo13 ****"
   allocate(Vo12(0:MaxN+1,0:MaxN+1),STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Not enough memory to allocate Vo13 ****"
end subroutine InitiliseVariables
!-------------------------------------!
! Destroy Variables                   !
!-------------------------------------!
subroutine DestroyVariables
   use variables
   use constants
   
   deallocate(P,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate P ****"
   
   deallocate(aUx1,STAT=istat)
   deallocate(bUx1,STAT=istat)
   deallocate(aUx2,STAT=istat)
   deallocate(bUx2,STAT=istat)
   deallocate(aUy1,STAT=istat)
   deallocate(bUy1,STAT=istat)
   deallocate(aUy2,STAT=istat)
   deallocate(bUy2,STAT=istat)
   
   deallocate(aVx1,STAT=istat)
   deallocate(bVx1,STAT=istat)
   deallocate(aVx2,STAT=istat)
   deallocate(bVx2,STAT=istat)
   deallocate(aVy1,STAT=istat)
   deallocate(bVy1,STAT=istat)
   deallocate(aVy2,STAT=istat)
   deallocate(bVy2,STAT=istat)
   
   deallocate(U,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate U ****"
   deallocate(V,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate V ****"
   deallocate(Uon,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate Uon1 ****"
   deallocate(Von,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate Von1 ****"

   deallocate(Uo12,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate Uo13 ****"
   deallocate(Vo12,STAT=istat)   
   if (istat.ne.0) write(*,*) "**** Can't deallocate Vo13 ****"
end subroutine DestroyVariables
