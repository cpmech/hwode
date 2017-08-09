C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON BECKER DOERING PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=5000,LJAC=5,LMAS=1,LE=8)
        PARAMETER (LWORK=ND*(LJAC+LMAS+3*LE+12)+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FBECK,JBECK,FBECKD,JBECKD,MBECK,SOLOUT,JBECKF
        COMMON /PRJECTION/RHO,FACT
        COMMON /COEFF/B(5000),RHO1
        R0=0.0D0
        DO I=1,ND-1
           R1=I**(2.D0/3.D0)
           B(I+1)=EXP((2*I-1.0D0)/(R1*R1+R1*R0+R0*R0))
           R0=R1
        END DO
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
        N=ND
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- THE FOLLOWING VALUES HAVE TO BE USED ONLY IF
C --- DC_DECSOL_BD.F  IS USED
        MLJAC=3
        MUJAC=1
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
        MLMAS=0
        MUMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        RHO=7.5D0
        RHO1=RHO
        X=0.0D0
        Y(1)=RHO
        DO I=2,N
           Y(I)=0.0D0
        END DO
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=1.0D-3*RTOL
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
        XEND=1.0D0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,15
C --- CALL OF THE SUBROUTINE RADAU5
        CALL RADAU5(N,FBECK,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBECK,IJAC,MLJAC,MUJAC,
     &                  MBECK,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        WRITE (8,*) Y(1)
        DO J=2,N-1,1000
           WRITE (8,*) Y(J)
        END DO
        WRITE (8,*) Y(N)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND*10.D0
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
           SUM=Y(1)
           DO I=2,N
              SUM=SUM+I*Y(I)
           END DO
           write (6,*) '  rho ',sum,rho
        IF (TARRAY(1).GT.500.) STOP
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        LOGICAL PROJECT
        COMMON /INTERN/XOUT
        COMMON /PRJECTION/RHO,FACT
        COMMON /COEFF/B(5000),RHO1
        PROJECT=.TRUE.
        IF (PROJECT) THEN
           IF (NR.EQ.1) FACT=N*(N+1.D0)*(1+2*N)/6.0D0
           SUM=0.D0
           DO I=N,1,-1
              SUM=SUM+I*Y(I)
           END DO
           ALAM=(SUM-RHO)/FACT
           DO I=1,N
              Y(I)=Y(I)-I*ALAM
           END DO
        END IF
        RETURN
        END
