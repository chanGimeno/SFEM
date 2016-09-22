PROGRAM test_hemPol

  USE hermite
  
  IMPLICIT NONE 

  INTEGER, PARAMETER :: p = 4
  INTEGER, PARAMETER :: dimx = 101
  DOUBLE PRECISION, DIMENSION(dimx) :: x
  DOUBLE PRECISION, DIMENSION(dimx) :: hermitePol
  DOUBLE PRECISION, DIMENSION(p,dimx) :: hermitePols

  CALL random_seed()
  CALL random_number(x)
  PRINT *, factorial(p)
  call h_k(p, x, dimx, hermitePol)
  call h_k_norm(p, x, dimx, hermitePol)
  CALL hermPol(p, x, dimx, hermitePols)

END PROGRAM test_hemPol
