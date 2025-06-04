C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON KS PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (BANDED JACOBIAN)
        PARAMETER (MMM=9)
        PARAMETER (NH=2**MMM,N=2*NH,ND=2*NH-2)
        PARAMETER (IJAC=1,MLJAC=0,MUJAC=0,IMAS=0)
        PARAMETER (LJAC=MLJAC+MUJAC+1,LE=2*MLJAC+MUJAC+1)
        PARAMETER (LWORK=ND*(LJAC+3*LE+12)+20,LIWORK=3*ND+20)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),U(N)
        COMMON/TRANS/QQ,UZERO
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FKS,JKS,SOLOUT
        LOGICAL DEBUG
C --- DATA FOR THE PROBLEM
        QQ=0.025D0
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
C --- OUTPUT DURING INTEGRATION ---
        IOUT=0
C --- INITIAL VALUES
        T=0.0D0
C --- INITIAL POSITIONS IN Y-SPACE ----
        AN=N
        DO I=1,N
         DELX=1.D0/AN
         X=DELX*(I-1)
         U1=MIN(X-0.0D0,0.1D0-X)
         U2=20.D0*(X-0.2D0)*(0.3D0-X)
         U3=MIN(X-0.6D0,0.7D0-X)
         U4=MIN(X-0.9D0,1.0D0-X)
         U(I)=16.D0*MAX(0.D0,U1,U2,U3,U4)
        END DO
C ---  FOURIER TRANSFORM ---
        CALL REALFT(U,NH,+1)
        DO I=1,N
         U(I)=U(I)/AN
        END DO
C --- INITIAL POSITION IN FOURIER MODES ---
        UZERO=U(1)
        DO I=1,ND
         Y(I)=U(I+2)
        END DO
        U(1)=UZERO
        U(2)=0.D0
        DO I=3,N
         U(I)=Y(I-2)
        END DO
        CALL REALFT(U,NH,-1)
        DO I=1,N
         U(I)=2.D0*U(I)
        END DO
C --- ENDPOINT OF INTEGRATION
         TEND=100.0D0
C --- REQUIRED TOLERANCE
         RTOL=TOLST
         ATOL=RTOL
         ITOL=0
C --- INITIAL STEP SIZE
         H=1.0D-6 
C --- SET DEFAULT VALUES
         DO I=1,20
          IWORK(I)=0
          WORK(I)=0.D0
         END DO
C --- CALL OF THE SUBROUTINE RADAU5
        CALL DTIME(TARRAY,TRESULT)
        DEBUG=.FALSE.
         CALL RADAU5(ND,FKS,T,Y,TEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JKS,IJAC,MLJAC,MUJAC,
     &                  JKS,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        CALL DTIME(TARRAY,TRESULT)
        WRITE(8,9921)(Y(I),I=1,10)
        WRITE(8,9921)(Y(I),I=11,100,5)
        WRITE(8,9921)(Y(I),I=101,1024,50)
 9921   FORMAT(1X,F22.16)
C --- PRINT STATISTICS
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(IWORK(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
C -------- NEW TOLERANCE ---
        TOLST=TOLST*TOLFC
        IF (TARRAY(1).GT.1000.) STOP
 30     CONTINUE
        STOP
        END

        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
          IF(MOD(NR,10).EQ.1) WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F7.4,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
