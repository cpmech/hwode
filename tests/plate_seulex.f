C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX ON PLATE PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
      PARAMETER(MX=8,MY=5,MACHS1=2,MACHS2=4)
        PARAMETER (ND=2*MX*MY,KM=12,LWORK=2*ND*ND+(KM+9)*ND+4*KM+20)
        PARAMETER (LIWORK=2*ND+KM+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FPLATE,JPLATE,JPLATS,JPLATSB,XPLATE,SOLOUT
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
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- PROBLEM IS non-AUTONOMOUS
        IFCN=1
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
        WORK(11)=10.D0
        WORK(12)=100.D0 
        WORK(13)=5.D0
C --- SECOND ORDER OPTION AND BANDED
           IJAC=1
C           IJAC=0
           IWORK(9)=N/2
           MLJAC=2*MX
           MUJAC=2*MX
C --- ENDPOINT OF INTEGRATION
        XEND=7.0D0
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE SEULEX 
        CALL SEULEX(N,FPLATE,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JPLATSB,IJAC,MLJAC,MUJAC,
     &                  FPLATE,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
      DO 12 K=1,N
 12   WRITE (8,*)Y(K)
        CALL DTIME(TARRAY,TRESULT)
C --- PRINT STATISTICS
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
        IF (TARRAY(1).GT.100.) GOTO 40
        TOLST=TOLST*TOLFC
 30     CONTINUE
 40     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
           WRITE (6,99) X,Y(1),Y(40),Y(80),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
