C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON ROBERTSON'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=3,NRDENS=3)
        PARAMETER (LWORK=8*NDGL+5*NRDENS+21,LIWORK=NRDENS+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK),RPAR(2)
        EXTERNAL FROBER,SOLOUT
        LOGICAL DEBUG
C --- DIMENSION OF THE SYSTEM
        N=NDGL
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES AND ENDPOINT OF INTEGRATION
        X=0.0D0
        Y(1)=1.0D0
        Y(2)=0.0D0
        Y(3)=0.0D0
        XEND=0.3D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-2
        ATOL=(1.0D-6)*RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0  
C --- INITIAL STEPSIZE
        WORK(7)=1.0D-6
C --- CALL OF THE SUBROUTINE DOPRI5   
        write(*,'(A)')''
        write(*,'(A)')'running dopri5.f test'
        DEBUG=.FALSE.
        CALL DOPRI5(N,FROBER,X,Y,XEND,
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
      write(*,'(A,3ES23.15)')'y =',Y(1),Y(2),Y(3)
      write(*,'(A,ES23.15)')'h =',WORK(7)
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
