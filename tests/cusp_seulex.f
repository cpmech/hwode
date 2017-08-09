C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX ON CUSP PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (ND=96,KM=12,LWORK=2*ND*ND+(KM+9)*ND+4*KM+20)
        PARAMETER (LIWORK=2*ND+KM+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        COMMON/NERVES/NNERV
        COMMON/DIFFCOEF/DIFFUS
        EXTERNAL FCUSP,JCUSP,SOLOUT
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_seulex')
        REWIND 8
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
c -------- INITIAL CONSTANTES --------
        NNERV=32
        DIFFUS=1.D0*NNERV*NNERV/144.D0
        N=3*NNERV
C --- PROBLEM IS AUTONOMOUS
        IFCN=0
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=0
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=3
        MUJAC=3
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C ---------- VAL INIT ------------
        X=0.D0
        XEND=1.1D0
        ANERV=NNERV
        DEL=2.D0*3.14159265358979324D0/ANERV
        DO 25 INERV=1,NNERV
            Y(3*INERV-2)=0.D0
            Y(3*INERV-1)=-2.D0*COS(INERV*DEL)
  25        Y(3*INERV)=2.D0*SIN(INERV*DEL)
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES 
        DO 10 I=1,13
  10    WORK(I)=0.D0
        DO 12 I=1,10
  12    IWORK(I)=0
        IWORK(4)=3
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE SEULEX 
        CALL SEULEX(N,FCUSP,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  FCUSP,IJAC,MLJAC,MUJAC,
     &                  FCUSP,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,LRCONT,LICONT,IDID)
C --- PRINT SOLUTION
      DO 20 I=1,N
        WRITE (8,*) Y(I)
  20   CONTINUE
C --- PRINT STATISTICS
        CALL DTIME(TARRAY,TRESULT)
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(IWORK(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I6,' jac=',I5,' step=',I5,
     &        ' accpt=',I5,' rejct=',I4,' dec=',I5,
     &        ' sol=',I6)
C -------- NEW TOLERANCE ---
        IF (TARRAY(1).GT.500.) STOP
        TOLST=TOLST*TOLFC
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
           WRITE (6,99) X,Y(1),Y(2),Y(3),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
