C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON PLATE PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
      PARAMETER(MX=8,MY=5,MACHS1=2,MACHS2=4)
        PARAMETER (ND=2*MX*MY,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FPLATE,JPLATE,JPLATS,JPLATSB,SOLOUT
        LOGICAL DEBUG
C -------- CONSTANTES ------------------
      COMMON/TRANS/NX,NXM1,NY,NYM1,NDEMI,NACHS1,NACHS2,NDUMMY,
     &       OMEGA,STIFFN,DELX,USH4,FAC,WEIGHT
      NX=MX
      NY=MY
      NACHS1=MACHS1
      NACHS2=MACHS2
      NXM1=NX-1
      NYM1=NY-1
      NDEMI=NX*NY
      OMEGA=1000.D0
      STIFFN=100.D0
      WEIGHT=200.D0
      DENOM=NX+1
      DELX=2.D0/DENOM
      USH4=1.D0/(DELX**4)
      FAC=STIFFN*USH4
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
        REWIND 8
C -------- TWO OPTIONS -------
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
      DO 1 I=1,N
  1   Y(I)=0.D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL*1.0D-3
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-2
C --- SET DEFAULT VALUES 
        DO I=1,20
           WORK(I)=0.D0
           IWORK(I)=0
           ISTAT(I)=0
        END DO
C --- SECOND ORDER OPTION AND BANDED
           IJAC=1
           IWORK(9)=N/2
           MLJAC=2*MX
           MUJAC=2*MX
c          MLJAC=N
C --- ENDPOINT OF INTEGRATION
        XEND=7.D0
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE RADAU5 
        DEBUG=.FALSE.
        CALL RADAU5(N,FPLATE,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
c     &                  JPLATE,IJAC,MLJAC,MUJAC,
     &                  JPLATSB,IJAC,MLJAC,MUJAC,
     &                  FPLATE,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
      DO 15 K=1,N
 15   WRITE (8,*)Y(K)
      CALL DTIME(TARRAY,TRESULT)
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(ISTAT(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
C -------- NEW TOLERANCE ---
        TOLST=TOLST*TOLFC
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
      WRITE (6,991) X,NR-1
 991  FORMAT(1X,'X =',F9.5,'    NSTEP =',I4)
C      WRITE(6,992)(Y(I),I=1,N,N/10)
 992  FORMAT(10X,10F14.10)
        RETURN
        END
