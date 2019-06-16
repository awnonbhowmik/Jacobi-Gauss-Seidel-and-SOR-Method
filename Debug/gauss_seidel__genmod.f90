        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 16 17:54:41 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSS_SEIDEL__genmod
          INTERFACE 
            SUBROUTINE GAUSS_SEIDEL(A,B,N,X0,TOL)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: A(N,N)
              REAL(KIND=4) :: B(N)
              REAL(KIND=4) :: X0(N)
              REAL(KIND=4) :: TOL
            END SUBROUTINE GAUSS_SEIDEL
          END INTERFACE 
        END MODULE GAUSS_SEIDEL__genmod
