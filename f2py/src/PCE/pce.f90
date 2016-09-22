!======================================================================
MODULE m_pce
!======================================================================

  USE hermite
    
CONTAINS

  !======================================================================
  SUBROUTINE pce(p, d, nump, Seq, U, Psi)
    
    IMPLICIT NONE 

    ! inout
    INTEGER, INTENT(in) :: p
    INTEGER, INTENT(in) :: d
    INTEGER, INTENT(in) :: nump
    DOUBLE PRECISION, DIMENSION(d), INTENT(in) :: U
    INTEGER, DIMENSION(nump,d), INTENT(in) :: Seq
    DOUBLE PRECISION, DIMENSION(nump), INTENT(out) :: Psi
    ! local
    INTEGER :: k, j
    DOUBLE PRECISION, DIMENSION(0:p-1,d) :: hermitePols

    ! Evaluate Hermite polynomial 
    CALL hermPol(p, U, d, hermitePols(0:,:))

    ! Evaluate PCEs
    Psi = 1.d0
    DO k = 1,nump
       ! PRINT *, Seq(k,:)
       DO j = 1,d
          Psi(k) = Psi(k) * hermitePols(Seq(k,j),j)
       END DO
    END DO
    
  END SUBROUTINE pce
  !======================================================================

!======================================================================
END MODULE m_pce
!======================================================================
