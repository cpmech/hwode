C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON HW-EQ1 EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=1,NRDENS=1)
        PARAMETER (LWORK=8*NDGL+5*NRDENS+21,LIWORK=NRDENS+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK),RPAR(1)
        EXTERNAL FCN,SOLOUT
        LOGICAL DEBUG
C --- DIMENSION OF THE SYSTEM
        N=NDGL
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=2
C --- INITIAL VALUES AND ENDPOINT OF INTEGRATION
        RPAR(1)=-50.0D0
        X=0.0D0
        Y(1)=0.0D0
        XEND=1.5D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-4
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0  
C --- DORIVAL: INITIAL STEPSIZE
        WORK(7)=1.0D-4
C --- DENSE OUTPUT IS USED
        IWORK(5)=NRDENS
        IWORK(21)=1
C --- CALL OF THE SUBROUTINE DOPRI5   
        write(*,'(A)')''
        write(*,'(A)')'running dopri5.f test'
        DEBUG=.FALSE.
        CALL DOPRI5(N,FCN,X,Y,XEND,
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
      write(*,'(A,4ES23.15)')'y =',Y(1)
      write(*,'(A,ES23.15)')'h =',WORK(7)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTD5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CON(5*ND),ICOMP(ND),RPAR(1)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) NR-1,X,Y(1)
           XOUT=X+0.1D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) NR-1,XOUT,CONTD5(1,XOUT,CON,ICOMP,ND)
              XOUT=XOUT+0.1D0
              GOTO 10
           END IF
        END IF
 99     FORMAT('step =',I4,', x =',F6.2,', y =',ES23.15)
        RETURN
        END
C
        SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C --- Hairer-Wanner Eq1
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N),RPAR(1)
        F(1)=RPAR(1)*(Y(1)-DCOS(X))
        RETURN
        END 

