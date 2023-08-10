!!!  HIGH ORDER IN SPACE AND TIME DEFERRED CORRECTION (EXPLICIT) 
!!!     RESIDUAL DISTRIBUTION METHOD 
!!!  DESIGNED FOR THE SYSTEM GIVEN BY THE EULER EQUATIONS in 1D and 2D
!!!
!!!  Authors:
!!!  Remi Abgrall (University of Zurich),
!!!  Paola Bacigaluppi (University of Zurich),
!!!  Svetlana Tokareva (University of Zurich)
!!!  Institute of Mathematics and Institute of Computational Sciences
!!!  University of Zurich
!!!  July 10, 2018
!!!  Correspondance:	remi.abgrall@math.uzh.ch
!!!  ------------------------------------------


!*----------------------------------------------------------------------------------------
!*This module (only) provides the function inverse to invert a square matrix
!*----------------------------------------------------------------------------------------
MODULE algebra
  USE PRECISION
  IMPLICIT NONE

  LOGICAL, SAVE, PRIVATE :: initialise = .FALSE.
  LOGICAL, SAVE, PRIVATE :: module_debug = .FALSE.

  LOGICAL, SAVE, PUBLIC               :: optim1_plain = .FALSE.
  LOGICAL, SAVE, PUBLIC               :: optim2_plain = .FALSE.
  LOGICAL, SAVE, PUBLIC               :: optim3_plain = .FALSE.
  LOGICAL, SAVE, PUBLIC               :: optim4_plain = .FALSE.
  LOGICAL, SAVE, PUBLIC               :: optim5_plain = .FALSE.
CONTAINS

  !------------------------------------------
  !!--  Inversion d'une matrice carree par LU et pivot partiel par ligne
  !!--  Ref : Numerical recipes in C
  !------------------------------------------
  FUNCTION Inverse( Mat)  RESULT (Mat1)
    CHARACTER(LEN = *), PARAMETER :: mod_name = "InverseLu"
    REAL(DP), DIMENSION(:, : ), INTENT(IN) :: Mat
    REAL(DP), DIMENSION(SIZE(Mat, dim = 1), SIZE( Mat, dim = 1)) :: Mat1, mat2
    REAL(DP), DIMENSION(SIZE(Mat, dim = 1))   :: col
    REAL(DP)                        :: d
    INTEGER    :: j, Nsize, l, zz
    INTEGER, DIMENSION(SIZE(Mat, dim = 1)) :: indx !-- CCE R.B. 2008/10/21



    IF (.NOT. initialise) THEN
       optim1_plain = .FALSE.
       optim2_plain = .TRUE.
       optim3_plain = .TRUE. !-- celle là est particulièrement efficace
       optim4_plain = .TRUE.

       optim2_plain = .FALSE. !-- ces 2 optimisations là mettre MHD en l'air, et pourtant elles sont meilleures en précusion, ce qui ne saute pas aux yeux au niveau des résidus %%%%%%
       optim4_plain = .FALSE.

       optim3_plain = .FALSE. !-- cette optimisation met MHD en l'air uniquement en -O5, pas en -O2
    END IF

    Nsize = SIZE(Mat, dim=1)
    mat2  = Mat

    CALL ludcmp(mat2, Nsize, indx, d) !-- indx OUT

    DO j = 1, Nsize
       col = 0.0_DP
       col(j) = 1.0_DP
       CALL luksb(mat2, Nsize, indx, col) !-- col IN OUT; indx IN
       Mat1(:, j) = col(: )
    END DO


  CONTAINS
    SUBROUTINE ludcmp(mat3, Nsize, indx, d)
      CHARACTER(LEN = *), PARAMETER :: mod_name = "ludcmp"

      INTEGER, INTENT(IN) :: Nsize
      REAL(DP), DIMENSION(Nsize, Nsize ), INTENT(INOUT) :: mat3 !-- Passer en :,: pose des pbs avec le -O5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER, DIMENSION(Nsize ), INTENT(OUT) :: indx !-- Passer en :,: pose des pbs avec le -O5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL(DP), INTENT(OUT) :: d

      REAL(DP), DIMENSION(Nsize) :: vv
      INTEGER :: i, j, k, imax, zz
      REAL(DP) :: Big, dum, somme, temp

      d = 1.0_dp
      DO i = 1, Nsize
         IF (optim1_plain) THEN
            Big = MAXVAL(ABS(mat3(i, :)))
         ELSE
            Big = 0.0_dp
            DO j = 1, Nsize
               temp = ABS(mat3(i, j))
               IF (temp > Big) THEN
                  Big = temp
               END IF
            END DO
         END IF
         IF (Big == 0.0_dp) THEN
            WRITE( *, *) mod_name, " ERREUR : Matrice singuliere, i==", i
            STOP
         END IF
         vv(i) = 1.0_DP / Big
      END DO

      DO j = 1, Nsize

         DO i = 1, j -1
            IF (optim2_plain) THEN
               mat3(i, j) = mat3(i, j) - SUM(mat3(i, 1: i - 1) * mat3(1: i - 1, j))
            ELSE
               somme = mat3(i, j)
               DO k = 1, i -1
                  somme = somme - mat3(i, k) * mat3(k, j)
               END DO
               mat3(i, j) = somme
            END IF
         END DO !-- boucle sur i


         Big = 0.0_dp
         imax = -1
         DO i = j, Nsize
            somme = mat3(i, j)
            DO k = 1, j -1
               somme = somme - mat3(i, k) * mat3(k, j)

            END DO

            mat3(i, j) = somme
            dum = vv(i) * ABS(somme)
            IF (dum >= Big) THEN
               Big = dum
               imax = i
            END IF
         END DO !-- boucle sur i

         IF (j /= imax) THEN
            DO k = 1, Nsize
               dum = mat3(imax, k)
               mat3(imax, k) = mat3(j, k)
               mat3(j, k) = dum
            END DO !-- boucle sur k
            d = - d
            vv(imax ) = vv(j)
         END IF



         indx(j) = imax

         IF (ABS(mat3(j, j)) <= 1.0e-20_dp) THEN
            mat3(j, j) = SIGN(1.0e-20_dp, mat3(j, j)) !-- CCE 2007/04/24 !-- CCE 2008/12/17
         END IF

         IF (j /= Nsize) THEN

            DO i = j + 1, Nsize
               mat3(i, j) = mat3(i, j) / mat3(j, j)
            END DO !-- boucle sur i
         END IF



      END DO !-- boucle sur j

    END SUBROUTINE ludcmp


    SUBROUTINE luksb(mat2, Nsize, indx, col)
      CHARACTER(LEN = *), PARAMETER :: mod_name = "luksb"

      INTEGER, INTENT(IN) :: Nsize
      REAL(DP),  DIMENSION(Nsize, Nsize), INTENT(IN) :: mat2 !-- Passer en :,: pose des pbs avec le -O5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER,  DIMENSION(Nsize), INTENT(IN) :: indx !-- Passer en :,: pose des pbs avec le -O5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL(DP), DIMENSION(Nsize), INTENT(INOUT) :: col !-- Passer en :,: pose des pbs avec le -O5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER ::  i, ii, ip, j
      REAL(DP)  :: somme

      ii = 1
      DO i = 1, Nsize
         ip = indx(i)
         IF (ip <= 0 .OR. ip > SIZE(col)) THEN
            PRINT *, mod_name, " ERREUR : indx(", i, ") invalide", ip
            STOP
         END IF
         somme = col(ip)
         col(ip) = col(i)
         IF (ii > 0) THEN
            IF (optim3_plain) THEN
               somme = somme - SUM(mat2(i, ii: i -1) * col(ii: i -1))
            ELSE
               DO j = ii, i -1
                  somme = somme - mat2(i, j) * col(j)
               END DO
            END IF
         ELSE
            IF (somme == 0.0_dp) THEN
               ii = i
            END IF
         END IF

         col(i) = somme
      END DO

      DO i = Nsize, 1, -1
         IF (optim4_plain) THEN
            col(i) = (col(i) - SUM(mat2(i, i + 1: Nsize) * col(i + 1: Nsize))) / mat2(i, i)
         ELSE
            somme = col(i)
            DO j = i + 1, Nsize
               somme = somme - mat2(i, j) * col(j)
            END DO
            col(i) = somme / mat2(i, i)
         END IF
      END DO
    END SUBROUTINE  luksb

  END FUNCTION Inverse
END MODULE algebra
