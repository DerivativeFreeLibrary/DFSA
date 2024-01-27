!     =================================================================
!     =================================================================
!     Module: Subroutines that define the problem
!     =================================================================
!     Define objective function, constaints and bounds
!     ******************************************************************
!     ******************************************************************

SUBROUTINE SETDIM(N,NI,NE)
	IMPLICIT NONE
	
	INTEGER				:: N, NI, NE

	n = 5
	ni = 5
	ne = 1

	RETURN
END SUBROUTINE SETDIM 

SUBROUTINE FUNCTINIT(N,X,LB,UB,nomefun,fglob)
	IMPLICIT NONE
	
	INTEGER				:: N
	DOUBLE PRECISION	:: X(N), LB(N), UB(N)
	DOUBLE PRECISION	:: fglob
	CHARACTER(LEN=40)	:: nomefun

	fglob = -30665.54d0
	nomefun = 'ex3_1_2 of the COCONUT benchmark'

	x(1) = 78.d0
	x(2) = 33.d0
	x(3) = 29.9953d0
	x(4) = 45.d0
	x(5) = 36.7758d0

	lb(1) = 78.d0; ub(1) = 102.d0 
	lb(2) = 33.d0; ub(2) = 45.d0
	lb(3) = 27.d0; ub(3) = 45.d0
	lb(4) = 27.d0; ub(4) = 45.d0
	lb(5) = 27.d0; ub(5) = 45.d0

	RETURN
END SUBROUTINE FUNCTINIT 

subroutine funob(n,x,f)

      implicit none

!     SCALAR ARGUMENTS
      integer n, i
      double precision f

!     ARRAY ARGUMENTS
      double precision x(n)

	f = 0.8356891d0*x(1)*x(5) + 37.293239d0*x(1) + 5.3578547d0*x(3)*x(3) - 40792.141d0

     return
end subroutine funob
!     ******************************************************************

subroutine fconstreq(n,m,x,c)
      
      implicit none

!     SCALAR ARGUMENTS
      integer n,m,i
!     ARRAY ARGUMENTS
      double precision c(m)
      double precision x(n)
      
	c(1) = 0.0056858d0*x(2)*x(5) - 0.0022053d0*x(3)*x(5) + 0.0006262d0*x(1)*x(4) - 6.665593d0

end subroutine fconstreq


!     ******************************************************************

subroutine fconstrin(n,m,x,c)
      
      implicit none

!     SCALAR ARGUMENTS
      integer n,m
!     ARRAY ARGUMENTS
      double precision x(n)
      double precision c(m)

! inequality constraints must be given as g(x) <= 0

	c(1) = 0.0022053d0*x(3)*x(5) - 0.0056858d0*x(2)*x(5) - 0.0006262d0*x(1)*x(4) - 85.334407d0
	c(2) = 0.0071317d0*x(2)*x(5) + 0.0021813d0*x(3)*x(3) + 0.0029955d0*x(1)*x(2) - 29.48751d0
	c(3) =-0.0071317d0*x(2)*x(5) - 0.0021813d0*x(3)*x(3) - 0.0029955d0*x(1)*x(2) + 9.48751d0
	c(4) = 0.0047026d0*x(3)*x(5) + 0.0019085d0*x(3)*x(4) + 0.0012547d0*x(1)*x(3) - 15.599039d0
	c(5) =-0.0047026d0*x(3)*x(5) - 0.0019085d0*x(3)*x(4) - 0.0012547d0*x(1)*x(3) + 10.699039d0
      
end subroutine fconstrin

!     ******************************************************************
