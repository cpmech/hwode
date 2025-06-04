C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON BEAM PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=80,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        COMMON/NNNN/NCOM,NNCOM,NSQ,NQUATR,DELTAS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FTIGE,SOLOUT
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
C ---------- CONSTANTS --------------
        N=40
        NN=2*N
        NCOM=N
        NSQ=N*N
        NQUATR=NSQ*NSQ
        NNCOM=NN
        AN=N
        DELTAS=1.D0/AN
C --- SET DEFAULT VALUES 
        DO 12 I=1,20
        WORK(I)=0.D0
  12    IWORK(I)=0
C --------- TUNED VALUES OF WORK PARAMETERS 
        WORK(3)=0.1
        WORK(4)=0.3
        WORK(5)=0.99
        WORK(6)=2.0
C --------- INITIAL VALUES -------------
        T=0.D0
        DO 1 I=1,NN
    1   Y(I)=0.D0
        TEND=5.D0
        H=1.D-3
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=0
C --------- SWITCH FOR OUTPUT --------
        IOUT=0
C --- SECOND ORDER OPTION
           IJAC=0
           IWORK(9)=N
           MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE RADAU5
        DEBUG=.FALSE.
        CALL RADAU5(NN,FTIGE,T,Y,TEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  FTIGE,IJAC,MLJAC,MUJAC,
     &                  FTIGE,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
        CALL DTIME(TARRAY,TRESULT)
C --- PRINT SOLUTION
        WRITE(8,9921)(Y(I),I=1,NN)
 9921   FORMAT(1X,F22.16)
        WRITE(8,*)TARRAY(1)
        WRITE (8,*) (IWORK(I),I=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (IWORK(I),I=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
C -------- NEW TOLERANCE ---
        IF(TARRAY(1).GT.1000.)STOP
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
