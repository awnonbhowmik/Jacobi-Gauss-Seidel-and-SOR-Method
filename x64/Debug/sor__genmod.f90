        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 16 17:53:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOR__genmod
          INTERFACE 
            SUBROUTINE SOR(A,B,N,X0,TOL,W)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: A(N,N)
              REAL(KIND=4) :: B(N)
              REAL(KIND=4) :: X0(N)
              REAL(KIND=4) :: TOL
              REAL(KIND=4) :: W
            END SUBROUTINE SOR
          END INTERFACE 
        END MODULE SOR__genmod
