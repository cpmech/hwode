C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
C link dr1_radau5 radau5 decsol dc_decsol   or
C link dr1_radau5 radau5 lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=2,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,JVPOL,SOLOUT
        LOGICAL DEBUG
C --- PARAMETER IN THE DIFFERENTIAL EQUATION
        RPAR=1.0D-6
C --- DIMENSION OF THE SYSTEM
        N=2
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
        Y(1)=2.0D0
        Y(2)=-0.6D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-4
        ATOL=1.0D0*RTOL
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
        CALL RADAU5(N,FVPOL,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVPOL,IJAC,MLJAC,MUJAC,
     &                  FVPOL,IMAS,MLMAS,MUMAS,
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
      write(*,'(A,ES23.15,ES23.15)')'y =',Y(1),Y(2)
      write(*,'(A,ES23.15)')'h =',H
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) NR-1,X,Y(1),Y(2)
           XOUT=0.2D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C --- CONTINUOUS OUTPUT FOR RADAU5
              WRITE (6,99) NR-1,XOUT,CONTR5(1,XOUT,CONT,LRC),
     &                     CONTR5(2,XOUT,CONT,LRC)
              XOUT=XOUT+0.2D0
              GOTO 10
           END IF
        END IF
 99     FORMAT('step =',I4,', x =',F5.2,', y =',2ES23.15)
        RETURN
        END
C
C
        SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/RPAR
        RETURN
        END 
C
C
        SUBROUTINE JVPOL(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=0.0D0
        DFY(1,2)=1.0D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/RPAR
        DFY(2,2)=(1.0D0-Y(1)**2)/RPAR
        RETURN
        END


