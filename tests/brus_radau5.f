C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON BRUSSELATOR 1D PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 
        PARAMETER (ND=1000,NL=2,NU=2)
        PARAMETER (LWORK=(7*NL+4*NU+16)*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
        COMMON/PARAM/N,N2,GAMMA,GAMMA2
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FBRUS,JBRUS,SOLOUT
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
        REWIND 8
c -------- INITIAL CONSTANTES --------
        PI=3.14159265358979324D0
        N=500
        N2=2*N
        USDELQ=(DBLE(N+1))**2
        GAMMA=0.02D0*USDELQ
        GAMMA2=2.D0*GAMMA
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
C ---------- VAL INIT ------------
        X=0.D0
        XEND=10.D0
        DO 1 I=1,N
        ANP1=N+1
        XI=I/ANP1
        Y(2*I)=3.D0
   1    Y(2*I-1)=1.D0+0.5D0*DSIN(2.D0*PI*XI)
C --- COMPUTE THE JACOBIAN NUMERICALLY
        IJAC=1
C --- JACOBIAN IS NOT FULL
        MLJAC=NL
        MUJAC=NU
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- CALL OF THE SUBROUTINE RADAU5 
C --- SET DEFAULT VALUES 
        DO I=1,20
           WORK(I)=0.D0
           IWORK(I)=0
           ISTAT(I)=0
        END DO
        CALL DTIME(TARRAY,TRESULT)
        CALL RADAU5(N2,FBRUS,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBRUS,IJAC,MLJAC,MUJAC,
     &                  FBRUS,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        CALL DTIME(TARRAY,TRESULT)
        WRITE(8,9921)(Y(I),I=1,995,7)
 9921   FORMAT(1X,F22.16)
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
        TOLST=TOLST*TOLFC
        IF (TARRAY(1).GT.500.) STOP
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),Y(3),NR-1
 99     FORMAT(1X,'X =',F11.4,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
