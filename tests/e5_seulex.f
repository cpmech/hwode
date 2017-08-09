C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX ON E5 PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (ND=4,KM=12,LWORK=2*ND*ND+(KM+9)*ND+4*KM+20)
        PARAMETER (LIWORK=2*ND+KM+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FE5,JE5,SOLOUT
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
        N=4
C --- PROBLEM IS AUTONOMOUS
        IFCN=0
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=1.76D-3
        Y(2)=0.0D0
        Y(3)=0.0D0
        Y(4)=0.0D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=1.7D-24
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
        XEND=10.0D0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,7
C --- CALL OF THE SUBROUTINE SEULEX 
        CALL SEULEX(N,FE5,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JE5,IJAC,MLJAC,MUJAC,
     &                  FE5,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        WRITE (8,*) Y(1)
        WRITE (8,*) Y(2)
        WRITE (8,*) Y(3)
        WRITE (8,*) Y(4)
c         write (6,*) y(2), y(3)+y(4)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND*100.D0
        CALL DTIME(TARRAY,TRESULT)
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(ISTAT(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
 91     FORMAT(' fcn=',I6,' jac=',I5,' step=',I5,
     &        ' accpt=',I5,' rejct=',I4,' dec=',I5,
     &        ' sol=',I6)
C -------- NEW TOLERANCE ---
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
