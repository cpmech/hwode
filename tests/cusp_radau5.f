C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON CUSP PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=96,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        COMMON/NERVES/NNERV
        COMMON/DIFFCOEF/DIFFUS
        EXTERNAL FCUSP,JCUSP,SOLOUT
        LOGICAL DEBUG
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
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
C --- COMPUTE THE JACOBIAN SWITCH
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
        DO 10 I=1,20
  10    WORK(I)=0.D0
        DO 12 I=1,20
  12    IWORK(I)=0
C --- ENDPOINT OF INTEGRATION
        XEND=1.1D0
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE RADAU5
        DEBUG=.FALSE.
        CALL RADAU5(N,FCUSP,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  FCUSP,IJAC,MLJAC,MUJAC,
     &                  FCUSP,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        CALL DTIME(TARRAY,TRESULT)
      DO 20 I=1,N
        WRITE (8,*) Y(I)
  20   CONTINUE
C --- PRINT STATISTICS
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(IWORK(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
           write (6,*) y(1),y(n/2),y(n)
C -------- NEW TOLERANCE ---
        IF (TARRAY(1).GT.500.) STOP
        TOLST=TOLST*TOLFC
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F6.2,'    Y =',2F18.7,'    NSTEP =',I4)
        RETURN
        END
