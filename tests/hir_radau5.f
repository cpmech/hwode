C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON HIRES PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=8,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FHIRES,JHIRES,SOLOUT
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
C --- DIMENSION OF THE SYSTEM
        N=8
C --- COMPUTE THE JACOBIAN
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
      Y (1) = 1.D0
      Y (2) = 0.D0
      Y (3) = 0.D0
      Y (4) = 0.D0
      Y (5) = 0.D0
      Y (6) = 0.D0
      Y (7) = 0.D0
      Y (8) = 0.0057D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL*1.D-4
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES 
        DO I=1,20
           WORK(I)=0.D0
           IWORK(I)=0
           ISTAT(I)=0
        END DO
C --- ENDPOINT OF INTEGRATION
        XEND=321.8122D0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,2
C --- CALL OF THE SUBROUTINE RADAU5
        DEBUG=.FALSE.
        CALL RADAU5(N,FHIRES,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JHIRES,IJAC,MLJAC,MUJAC,
     &                  FHIRES,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
      DO 15 K=1,8
 15   WRITE (8,*)Y(K)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND+100.D0
        CALL DTIME(TARRAY,TRESULT)
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
           WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F6.2,'    Y =',2F18.7,'    NSTEP =',I4)
        RETURN
        END
