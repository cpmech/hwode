C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT HW-EQ1 EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
C link dr1_radau5 radau5 decsol dc_decsol   or
C link dr1_radau5 radau5 lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=1,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FCN,JAC,SOLOUT
        LOGICAL DEBUG
C --- PARAMETER IN THE DIFFERENTIAL EQUATION
        LAMB=-50.0D0
C --- DIMENSION OF THE SYSTEM
        N=1
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=1.5D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-4
        ATOL=1.0D0*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-4 
C --- SET DEFAULT VALUES 
        DO I=1,20
           IWORK(I)=0
           WORK(I)=0.D0
        END DO
C --- CALL OF THE SUBROUTINE RADAU5
        write(*,'(A)')''
        write(*,'(A)')'running radau5.f test'
        DEBUG=.TRUE.
        CALL RADAU5(N,FCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC,IJAC,MLJAC,MUJAC,
     &                  FCN,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,LAMB,IPAR,IDID)
C --- PRINT FINAL SOLUTION
      write(*,'(A,I0)')'Number of function evaluations   = ',1+IWORK(14)
      write(*,'(A,I0)')'Number of Jacobian evaluations   = ',IWORK(15)
      write(*,'(A,I0)')'Number of factorizations         = ',IWORK(19)
      write(*,'(A,I0)')'Number of lin sys solutions      = ',IWORK(20)
      write(*,'(A,I0)')'Number of performed steps        = ',IWORK(16)
      write(*,'(A,I0)')'Number of accepted steps         = ',IWORK(17)
      write(*,'(A,I0)')'Number of rejected steps         = ',IWORK(18)
      write(*,'(A,I0)')'Number of iterations (maximum)   = ',IWORK(21)
      write(*,'(A,ES23.15)')'y =',Y(1)
      write(*,'(A)')'.'
      write(*,'(A)',advance='no')'test result: ok. 1 passed; 0 failed;'
      write(*,'(A)',advance='no')' 0 ignored; 0 measured;'
      write(*,'(A)')' 0 filtered out; finished in 0.00s'
      write(*,'(A)')''
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,LAMB,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        RETURN
        END
C
C
        SUBROUTINE FCN(N,X,Y,F,LAMB,IPAR)
C --- RIGHT-HAND SIDE OF HW-EQ1
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=LAMB*(Y(1)-DCOS(X))
        RETURN
        END 
C
C
        SUBROUTINE JAC(N,X,Y,DFY,LDFY,LAMB,IPAR)
C --- JACOBIAN OF HW-EQ1
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=LAMB
        RETURN
        END


