C * * * * * * * * * * * * * * * * * * * * * * * * *
C                 DRIVER FOR ODEX2 
C * * * * * * * * * * * * * * * * * * * * * * * * *
compile odex2
cfeh dr_odex2 odex2
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=2,KM=9,NRDENS=2,
     &     LWORK=NDGL*(2*KM+6)+5*KM+20+(KM*(2*KM+5)+6)*NRDENS,
     &     LIWORK=KM+20+NRDENS)
        DIMENSION Y(NDGL),YP(NDGL),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL TWOB,SOLOUT
C --- DIMENSION OF THE SYSTEM
        N=2
C --- OUTPUT ROUTINE AND DENSE OUTPUT IS USED DURING INTEGRATION
        IOUT=2
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=0.5D0
        Y(2)=0.0D0
        YP(1)=0.0D0
        YP(2)=SQRT(3.D0)
C --- ENDPOINT OF INTEGRATION
        XEND=20.0D0
C --- REQUIRED (RELATIVE) TOLERANCE
        TOL=1.0D-10
        ITOL=0
        RTOL=TOL
        ATOL=TOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,10
        IWORK(I)=0
  10    WORK(I)=0.D0  
        H=0.01D0   
C --- IF DENSE OUTPUT IS REQUIRED
        IWORK(8)=N
C --- CALL OF THE SUBROUTINE DOPRI5   
        CALL ODEX2(N,TWOB,X,Y,YP,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) TOL
 90     FORMAT('       tol=',D8.2)
        WRITE (6,91) (IWORK(J),J=17,20)
 91     FORMAT(' fcn=',I5,' step=',I4,' accpt=',I4,' rejct=',I3)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,YP,N,CON,NCON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTD5", THE CONTINUOUS COLLOCATION SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),YP(N),CON(NCON),ICOMP(ND)
        COMMON /INTERN/XOUT  
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=X+1.0D0
        ELSE
C           WRITE (6,99) X,Y(1),Y(2),NR-1
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              SOL1=CONTX2(1,XOUT,CON,NCON,ICOMP,ND)
              SOL2=CONTX2(2,XOUT,CON,NCON,ICOMP,ND)
              WRITE (6,99) XOUT,SOL1,SOL2,NR-1
              XOUT=XOUT+1.0D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE TWOB(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        RAD=Y(1)**2+Y(2)**2
        RAD=RAD*SQRT(RAD)
        F(1)=-Y(1)/RAD
        F(2)=-Y(2)/RAD
        RETURN
        END 

