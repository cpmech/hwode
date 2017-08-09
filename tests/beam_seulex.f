C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX ON BEAM PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (ND=80,KM=12,KM2=2+KM*(KM+3)/2,NRDENS=0)
        PARAMETER (LWORK=2*ND*ND+(KM+8)*ND+4*KM+20+KM2*NRDENS)
        PARAMETER (LIWORK=2*ND+KM+20+NRDENS) 
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        COMMON/NNNN/NCOM,NNCOM,NSQ,NQUATR,DELTAS
        EXTERNAL FTIGE,SOLOUT
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
C --- PROBLEM IS NON-AUTONOMOUS
        IFCN=1
C --- SECOND ORDER OPTION
           IJAC=0
           IWORK(9)=N
           MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
        CALL DTIME(TARRAY,TRESULT)
C --- CALL OF THE SUBROUTINE SEULEX 
        CALL SEULEX(NN,FTIGE,IFCN,T,Y,TEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  FTIGE,IJAC,MLJAC,MUJAC,
     &                  FTIGE,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
        CALL DTIME(TARRAY,TRESULT)
C --- PRINT SOLUTION
        WRITE(8,9921)(Y(I),I=1,NN)
 9921   FORMAT(1X,F22.16)
        WRITE(8,*)TARRAY(1)
        WRITE (8,*) (IWORK(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
C -------- NEW TOLERANCE ---
        IF(TARRAY(1).GT.500.)STOP
        TOLST=TOLST*TOLFC
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,
     &                     RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
           WRITE (6,99) X,Y(1),Y(2),Y(3),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
