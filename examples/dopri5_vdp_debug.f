C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON VAN DER POL EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=2,NRDENS=2)
        PARAMETER (LWORK=8*NDGL+5*NRDENS+21,LIWORK=NRDENS+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,SOLOUT
        LOGICAL DEBUG
C --- DIMENSION OF THE SYSTEM
        N=NDGL
        RPAR=0.003D0
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES AND ENDPOINT OF INTEGRATION
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=0.0D0
        XEND=2.0D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-3
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0  
C --- INITIAL STEPSIZE
        WORK(7)=1.0D-4
C --- Stiffness Detection all the time
        IWORK(4)=1
C --- CALL OF THE SUBROUTINE DOPRI5   
        write(*,'(A)')''
        write(*,'(A)')'running dopri5.f test'
        DEBUG=.TRUE.
        CALL DOPRI5(N,FVPOL,X,Y,XEND,
     &              RTOL,ATOL,ITOL,
     &              SOLOUT,IOUT,DEBUG,
     &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
      write(*,'(A)',advance='no')'DoPri5: Dormand-Prince method '
      write(*,'(A)')'(explicit, order 5(4), embedded)'
      write(*,'(A,I0)')'Number of function evaluations   = ',IWORK(17)
      write(*,'(A,I0)')'Number of performed steps        = ',IWORK(18)
      write(*,'(A,I0)')'Number of accepted steps         = ',IWORK(19)
      write(*,'(A,I0)')'Number of rejected steps         = ',IWORK(20)
      write(*,'(A,2ES23.15)')'y =',Y(1),Y(2)
      write(*,'(A,ES23.15)')'h =',WORK(7)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
C        WRITE (6,99) NR-1,X,Y(1),Y(2)
 99     FORMAT('step =',I4,', x =',F5.2,', y =',2ES23.15)
        RETURN
        END
C
        SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EPS=RPAR
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/EPS
        RETURN
        END 
