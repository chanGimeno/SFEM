!======================================================================
MODULE hermite
!======================================================================

CONTAINS

  !======================================================================
  RECURSIVE FUNCTION factorial(n) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(in) :: n
    INTEGER :: res

    IF (n .EQ. 0) THEN
       res = 1
    ELSE
       res = n * factorial(n - 1)
    END IF

  END FUNCTION factorial
  !======================================================================

  ! !======================================================================
  ! RECURSIVE FUNCTION h_k(k, x, dimx) RESULT(hermite)

  !   IMPLICIT NONE 
  !   INTEGER, INTENT(in) :: k
  !   INTEGER, INTENT(in) :: dimx
  !   DOUBLE PRECISION, INTENT(in), DIMENSION(dimx) :: x
  !   DOUBLE PRECISION, DIMENSION(dimx) :: hermite

  !   IF (k == 0) THEN
  !      hermite(:) = 1.d0
  !   ELSE IF (k == 1) THEN
  !      hermite(:) = x(:)
  !   ELSE        
  !      hermite = x*h_k(k-1,x,dimx)-(k-1)*h_k(k-2,x,dimx)
  !   END IF

  ! END FUNCTION h_k
  ! !======================================================================

  ! !======================================================================
  ! FUNCTION h_k_norm(k, x, dimx) RESULT(hermiteNorm)

  !   IMPLICIT NONE 
  !   INTEGER, INTENT(in) :: k
  !   INTEGER, INTENT(in) :: dimx
  !   DOUBLE PRECISION, INTENT(in), DIMENSION(dimx) :: x
  !   DOUBLE PRECISION, DIMENSION(dimx) :: hermiteNorm
  !   INTEGER, external :: factorial
  !   DOUBLE PRECISION, DIMENSION(dimx), external :: h_k

  !   hermiteNorm = h_k(k,x,dimx)/(1.d0*factorial(k))

  ! END FUNCTION h_k_norm
  ! !======================================================================

  !======================================================================
  RECURSIVE SUBROUTINE h_k(k, x, dimx, hermite)

    IMPLICIT NONE 
    ! INOUT
    INTEGER, INTENT(in) :: k
    INTEGER, INTENT(in) :: dimx
    DOUBLE PRECISION, INTENT(in), DIMENSION(dimx) :: x
    DOUBLE PRECISION, DIMENSION(dimx), INTENT(inout) :: hermite
    ! LOCAL
    DOUBLE PRECISION, DIMENSION(dimx) :: h_k_m1, h_k_m2

    IF (k == 0) THEN
       hermite(:) = 1.d0
    ELSE IF (k == 1) THEN
       hermite(:) = x(:)
    ELSE        
       CALL h_k(k-1, x, dimx, h_k_m1)
       CALL h_k(k-2, x, dimx, h_k_m2)
       hermite = x*h_k_m1-(k-1)*h_k_m2
    END IF

  END SUBROUTINE h_k
  !======================================================================

  !======================================================================
  SUBROUTINE h_k_norm(k, x, dimx, hermiteNorm)

    IMPLICIT NONE 
    INTEGER, INTENT(in) :: k
    INTEGER, INTENT(in) :: dimx
    DOUBLE PRECISION, DIMENSION(dimx), INTENT(in) :: x
    DOUBLE PRECISION, DIMENSION(dimx), INTENT(out) :: hermiteNorm
    ! INTEGER, external :: factorial

    CALL h_k(k, x, dimx, hermiteNorm)
    hermiteNorm = hermiteNorm/SQRT(1.d0*factorial(k))

  END SUBROUTINE h_k_norm
  !======================================================================

  !======================================================================
  SUBROUTINE hermPol2(p, x, dimx, hermitePols)

    IMPLICIT NONE 
    ! inout
    INTEGER, INTENT(in) :: p
    INTEGER, INTENT(in) :: dimx
    DOUBLE PRECISION, INTENT(in), DIMENSION(dimx) :: x
    DOUBLE PRECISION, DIMENSION(p,dimx), INTENT(out) :: hermitePols
    ! local
    INTEGER :: k
    DOUBLE PRECISION, DIMENSION(dimx) :: hermitePol

    DO k=0, p-1
       CALL h_k_norm(k, x, dimx, hermitePol)
       hermitePols(k+1,:) = hermitePol
    END DO

  END SUBROUTINE hermPol2
  !======================================================================

  !======================================================================
  SUBROUTINE hermPol(p, x, dimx, hermitePols)

    IMPLICIT NONE 

    ! inout
    INTEGER, INTENT(in) :: p
    INTEGER, INTENT(in) :: dimx
    DOUBLE PRECISION, INTENT(in), DIMENSION(dimx) :: x
    DOUBLE PRECISION, DIMENSION(p,dimx), INTENT(out) :: hermitePols
    ! local
    INTEGER :: k
    DOUBLE PRECISION, DIMENSION(dimx) :: h_k_m1, h_k_m2
    DOUBLE PRECISION, DIMENSION(dimx) :: hermitePol

    DO k=0, p-1
       IF (k == 0) THEN
          hermitePol = 1.d0
          h_k_m2 = hermitePol
       ELSE IF (k == 1) THEN
          hermitePol = x
          h_k_m1 = hermitePol
       ELSE        
          hermitePol = x*h_k_m1 - (k-1)*h_k_m2
          h_k_m2 = h_k_m1
          h_k_m1 = hermitePol
       END IF
       hermitePols(k+1,:) = hermitePol/SQRT(1.d0*factorial(k))
    END DO

  END SUBROUTINE hermPol
  !======================================================================

!======================================================================
END MODULE hermite
!======================================================================
