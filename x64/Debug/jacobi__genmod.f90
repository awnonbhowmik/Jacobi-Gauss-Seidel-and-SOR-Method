        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 16 17:53:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE JACOBI__genmod
          INTERFACE 
            SUBROUTINE JACOBI(A,B,N,X0,TOL)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: A(N,N)
              REAL(KIND=4) :: B(N)
              REAL(KIND=4) :: X0(N)
              REAL(KIND=4) :: TOL
            END SUBROUTINE JACOBI
          END INTERFACE 
        END MODULE JACOBI__genmod
