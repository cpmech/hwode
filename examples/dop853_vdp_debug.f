C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
cfeh dr_dop853 dop853
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=2,NRD=2)
        PARAMETER (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,SOLOUT
        LOGICAL DEBUG
C --- DIMENSION OF THE SYSTEM
        N=2
        RPAR=1.0D-3
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=2
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED (RELATIVE) TOLERANCE
        TOL=1.0D-9
        ITOL=0
        RTOL=TOL
        ATOL=TOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,10
        IWORK(I)=0
  10    WORK(I)=0.D0   
        IWORK(5)=N
        IWORK(4)=1
C --- DORIVAL: INITIAL STEPSIZE
        WORK(7)=1.0D-6
C --- CALL OF THE SUBROUTINE DOPRI8   
        write(*,'(A)')''
        write(*,'(A)')'running dop853.f test'
        DEBUG=.TRUE.
        CALL DOP853(N,FVPOL,X,Y,XEND,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
      write(*,'(A)',advance='no')'DoPri8: Dormand-Prince method '
      write(*,'(A)')'(explicit, order 8(5,3), embedded)'
      write(*,'(A,I0)')'Number of function evaluations   = ',IWORK(17)
      write(*,'(A,I0)')'Number of performed steps        = ',IWORK(18)
      write(*,'(A,I0)')'Number of accepted steps         = ',IWORK(19)
      write(*,'(A,I0)')'Number of rejected steps         = ',IWORK(20)
      write(*,'(A,ES23.15,ES23.15)')'y =',Y(1),Y(2)
      write(*,'(A,ES23.15)')'h =',WORK(7)
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CON(8*ND),ICOMP(ND)
        COMMON /INTERN/XOUT
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

