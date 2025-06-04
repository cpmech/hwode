C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SDIRK4 AT VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
clink dr_sdirk4 sdirk4 decsol
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SDIRK4 (FULL JACOBIAN)
        PARAMETER (ND=2,LWORK=2*ND*ND+12*ND+7,LIWORK=2*ND+7)
        PARAMETER (LRCONT=5*ND+2)
        COMMON /CONT/ICONT(4),RCONT(LRCONT)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
        EXTERNAL FVPOL,JVPOL,SOLOUT
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
        Y(2)=-0.66D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-5
        ATOL=1.0D0*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES 
        DO 10 I=1,7
        IWORK(I)=0
  10    WORK(I)=0.D0
          IWORK(4)=0
C --- CALL OF THE SUBROUTINE SDIRK4
        CALL SDIRK4(N,FVPOL,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVPOL,IJAC,MLJAC,MUJAC,
     &                  FVPOL,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,LRCONT,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) RTOL
 90     FORMAT('       rtol=',D8.2)
        WRITE (6,91) NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTS4"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.1D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT-1.D-10) THEN
C --- CONTINUOUS OUTPUT FOR SDIRK4
              WRITE (6,99) XOUT,CONTS4(1,XOUT),CONTS4(2,XOUT),NR-1
              XOUT=XOUT+0.1D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I5)
        RETURN
        END
C
C
        SUBROUTINE FVPOL(N,X,Y,F)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EPS=1.0D-6
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/EPS
        RETURN
        END 
C
C
        SUBROUTINE JVPOL(N,X,Y,DFY,LDFY)
C --- JACOBIAN OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        EPS=1.0D-6
        DFY(1,1)=0.0D0
        DFY(1,2)=1.0D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/EPS
        DFY(2,2)=(1.0D0-Y(1)**2)/EPS
        RETURN
        END


