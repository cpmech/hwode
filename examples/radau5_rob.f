C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT ROBERTSON'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
C link dr1_radau5 radau5 decsol dc_decsol   or
C link dr1_radau5 radau5 lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=3,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FROBER,JROBER,SOLOUT
        LOGICAL DEBUG
C --- DIMENSION OF THE SYSTEM
        N=3
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
        Y(1)=1.0D0
        Y(2)=0.0D0
        Y(3)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=0.3D0
C --- REQUIRED TOLERANCE (ATOL must be must smaller)
        RTOL=1.0D-2
        ATOL=(1.0D-6)*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- SET DEFAULT VALUES 
        DO I=1,20
           IWORK(I)=0
           WORK(I)=0.D0
        END DO
C --- CALL OF THE SUBROUTINE RADAU5
        write(*,'(A)')''
        write(*,'(A)')'running radau5.f test'
        DEBUG=.FALSE.
        CALL RADAU5(N,FROBER,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JROBER,IJAC,MLJAC,MUJAC,
     &                  FROBER,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
      write(*,'(A)',advance='no')'Radau5: Radau method (Radau IIA) '
      write(*,'(A)')'(implicit, order 5, embedded)'
      write(*,'(A,I0)')'Number of function evaluations   = ',1+IWORK(14)
      write(*,'(A,I0)')'Number of Jacobian evaluations   = ',IWORK(15)
      write(*,'(A,I0)')'Number of factorizations         = ',IWORK(19)
      write(*,'(A,I0)')'Number of lin sys solutions      = ',IWORK(20)
      write(*,'(A,I0)')'Number of performed steps        = ',IWORK(16)
      write(*,'(A,I0)')'Number of accepted steps         = ',IWORK(17)
      write(*,'(A,I0)')'Number of rejected steps         = ',IWORK(18)
      write(*,'(A,I0)')'Number of iterations (maximum)   = ',IWORK(21)
      write(*,'(A,3ES23.15)')'y =',Y(1),Y(2),Y(3)
      write(*,'(A,ES23.15)')'h =',H
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        WRITE (6,99) NR-1,X,Y(1),Y(2),Y(3)
 99     FORMAT('step =',I4,', x =',F5.2,', y =',3ES23.15)
        RETURN
        END
C
C
        SUBROUTINE FROBER(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF ROBERTSON EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=-0.04D0*Y(1)+1.D4*Y(2)*Y(3)
        F(3)=3.D7*Y(2)*Y(2)
        F(2)=-F(1)-F(3)
        RETURN
        END 
C
        SUBROUTINE JROBER(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF ROBERTSON EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        PROD1=1.0D4*Y(2)
        PROD2=1.0D4*Y(3)
        PROD3=6.0D7*Y(2)
        DFY(1,1)=-0.04D0
        DFY(1,2)=PROD2
        DFY(1,3)=PROD1
        DFY(2,1)=0.04D0
        DFY(2,2)=-PROD2-PROD3
        DFY(2,3)=-PROD1
        DFY(3,1)=0.D0
        DFY(3,2)=PROD3
        DFY(3,3)=0.D0
        RETURN
        END
