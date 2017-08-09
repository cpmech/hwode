C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX AT OREGO PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (ND=3,KM=12,LWORK=2*ND*ND+(KM+9)*ND+4*KM+20)
        PARAMETER (LIWORK=2*ND+KM+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FOREGON,JOREGON,SOLOUT
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
        N=3
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
        Y(1)=1.0D0
        Y(2)=2.0D0
        Y(3)=3.0D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=1.0D-6*RTOL
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
        XEND=30.0D0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,12
C --- CALL OF THE SUBROUTINE SEULEX 
        CALL SEULEX(N,FOREGON,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JOREGON,IJAC,MLJAC,MUJAC,
     &                  FOREGON,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        WRITE (8,*) Y(1)
        WRITE (8,*) Y(2)
        WRITE (8,*) Y(3)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND+30.D0
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
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTR5" (OR "CONTS4")
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
        COMMON /INTERN/XOUT
           WRITE (6,99) X,Y(1),Y(2),Y(3),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
