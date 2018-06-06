! Description:
! random_normal: returns a normally distributed pseudo-random number with zero mean and unit variance
! 

FUNCTION random_normal() RESULT (ran_norm)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.   This version uses the default
!  uniform random number generator which is in your fortran library.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

!  Fortran 90 version by Alan Miller (alan @ mel.dms.csiro.au)

IMPLICIT NONE
REAL :: ran_norm

!     Local variables
REAL, PARAMETER :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,   &
                   half = 0.5, r1 = 0.27597, r2 = 0.27846
REAL            :: u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
ran_norm = v/u
RETURN

END FUNCTION random_normal

    
! choose k subjects from n randomly
! omega(i) = 1 if i'th subject is chosen, 0 if not
SUBROUTINE RanSample (omega, n, k)
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER :: n
    INTEGER :: k
    INTEGER :: omega (n)
    ! - - - local declarations - - -
    integer :: i,j
    integer :: index(n)
    double precision :: u    ! uniform random variable
    integer :: temp
    
    omega = 0
    index = (/ (i,i=1,n) /)
    do i=1,k
        call random_number(u)
        j = int(float(n-i+1)*u) + i
        temp = index(i)
        index(i) = index(j)
        index(j) = temp
    end do
    omega(index(1:k)) = 1
END SUBROUTINE RanSample


subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
    
    
subroutine vectorsquare(a,n,c)
!============================================================
! Calculate a*a', where a is a vector
! ----------------------------------------------------------
! input ...
! a(n) - vector a
! n    - dimension
! output ...
! c(n,n) - a*a'
!!!!!  next step: modify c to storage compression format
!!!!!  c is a positive definite matrix
!!!!!    stored by rows in lower triangular form 
!!!!!    as a one dimensional array, in the sequence
!===========================================================
implicit none 
integer n
double precision a(n), c(n,n)

integer i, j

do i=1,n-1
    c(i,i) = a(i)*a(i)
    do j=i+1,n
        c(i,j) = a(i)*a(j)
        c(j,i) = c(i,j)
    end do
end do
c(n,n) = a(n)*a(n)
end subroutine vectorsquare

subroutine vectorproduct(a,b,n,c)
!============================================================
! Calculate a*b', where a,b are vectors
! ----------------------------------------------------------
! input ...
! a(n) - vector a
! b(n) - vector b
! n    - dimension
! output ...
! c(n,n) - a*b'
!===========================================================
implicit none 
integer n
double precision a(n), b(n), c(n,n)

integer i, j

do i=1,n
    do j=1,n
        c(i,j) = a(i)*b(j)
    end do
end do
    end subroutine vectorproduct

subroutine synvectorproduct(a,b,n,c)
!============================================================
! Calculate a*b'+b*a', where a,b are vectors
! ----------------------------------------------------------
! input ...
! a(n) - vector a
! b(n) - vector b
! n    - dimension
! output ...
! c(n,n) - a*b'+b*a'
!===========================================================
implicit none 
integer n
double precision a(n), b(n), c(n,n)

integer i, j

do i=1,n
    do j=1,n
        c(i,j) = a(i)*b(j) + a(j)*b(i)
    end do
end do
end subroutine synvectorproduct
    
subroutine vectorproduct2(a,b,n,m,c)
!============================================================
! Calculate a*b', where a,b are vectors
! ----------------------------------------------------------
! input ...
! a(n) - vector a
! b(m) - vector b
! n,m    - dimension
! output ...
! c(n,m) - a*b'
!===========================================================
implicit none 
integer n, m
double precision a(n), b(m), c(n,m)

integer i, j

do i=1,n
    do j=1,m
        c(i,j) = a(i)*b(j)
    end do
end do
end subroutine vectorproduct2
    
subroutine syminverse(a,c,n)
!============================================================
! Inverse symmetric matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! Both a and c are symmetric
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)

double precision aa(n*(n+1)/2), cc(n*(n+1)/2)
integer i, j, k
integer nullty, ifault

k = 0
do i=1,n
	do j=1,i
		k = k + 1
		aa(k) = a(i,j)
	end do
end do
call syminv(aa, n, cc, nullty, ifault)
k = 0
do i = 1, n
	do j=1,i-1
		k = k + 1
		c(i,j) = cc(k)
		c(j,i) = cc(k)
	end do
	k = k + 1
	c(i,i) = cc(k)
end do
end subroutine syminverse
