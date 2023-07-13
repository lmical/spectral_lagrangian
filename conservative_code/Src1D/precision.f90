MODULE PRECISION !*It just defines what single, double and quadruple precisions will be (otherwise there is the risk that they may be compiler or machine dependent)
  USE, INTRINSIC :: iso_fortran_env
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = kind(1.0)!REAL32
  INTEGER, PARAMETER :: dp = kind(1.d0)!REAL64
  INTEGER, PARAMETER :: qp = 2*kind(1.d0)!REAL128
END MODULE PRECISION

