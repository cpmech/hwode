C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RETARD 
C * * * * * * * * * * * * * * *
compile retard
cfeh dr_retard retard
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=3,NGRID=11,LWORK=8*NDGL+21+NGRID,LIWORK=20)
        PARAMETER (NRDENS=1,LRCONT=600,LICONT=NRDENS+1)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
        COMMON /CORER/RCONT(LRCONT) 
        COMMON /COREI/ICONT(LICONT)
        EXTERNAL FCN,SOLOUT
C --- DIMENSION OF THE SYSTEM
        N=NDGL
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES AND ENDPOINT OF INTEGRATION
        RPAR=0.1D0
        X=0.0D0
        Y(1)=5.0D0
        Y(2)=0.1D0
        Y(3)=1.0D0
        XEND=40.D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-5
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0  
C --- SECOND COMPONENT USES RETARDED ARGUMENT
        IWORK(5)=NRDENS
        ICONT(2)=2
C ---  USE AS GRID-POINTS
        IWORK(6)=NGRID
        DO 12 I=1,NGRID-1
  12    WORK(20+I)=I
        WORK(20+NGRID)=20.D0
C --- CALL OF THE SUBROUTINE RETARD   
        CALL RETARD(N,FCN,X,Y,XEND,
     &              RTOL,ATOL,ITOL,
     &              SOLOUT,IOUT,
     &              WORK,LWORK,IWORK,LIWORK,LRCONT,LICONT,
     &              RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) Y(1),Y(2),Y(3)
 99     FORMAT(1X,'X = XEND     Y =',3E18.10)
C --- PRINT STATISTICS
        WRITE (6,91) RTOL,(IWORK(J),J=17,20)
 91     FORMAT('     tol=',D8.2,'   fcn=',I5,' step=',I4,
     &             ' accpt=',I4,' rejct=',I3)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N)
        EXTERNAL PHI
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),NR-1
           XOUT=X+5.D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) X,Y(1),NR-1
              XOUT=XOUT+5.D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F6.2,'    Y =',E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EXTERNAL PHI
        Y2L1=YLAG(2,X-1.D0,PHI,RPAR,IPAR)
        Y2L10=YLAG(2,X-10.D0,PHI,RPAR,IPAR)
        F(1)=-Y(1)*Y2L1+Y2L10
        F(2)=Y(1)*Y2L1-Y(2)
        F(3)=Y(2)-Y2L10
        RETURN
        END
C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        IF (I.EQ.2) PHI=RPAR
        RETURN
        END
