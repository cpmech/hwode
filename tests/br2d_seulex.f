C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR SEULEX ON BRUSS-2D
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (NSD=128,ND=NSD*NSD*2,KM=12,LJAC=1,LMAS=0,LE=1)
        PARAMETER (LWORK=ND*(LJAC+LMAS+LE+KM+8)+4*KM+20)
        PARAMETER (LIWORK=2*ND+KM+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FBRUS,JBRUS,SOLOUT
        COMMON/TRANS/ALF,NS,NSSQ,NSNSM1,NSM1SQ
        COMMON /FOURIER/TCOS(512),NF(2),ALPH,NDIM,NSF,NSSQF
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_seulex')
        REWIND 8
C ----- DIMENSIONS --------
      NS=128
      NSSQ=NS*NS
      NSNSM1=NS*(NS-1)
      N=2*NSSQ 
      ALF=1.0D-1
      ALPH=ALF
      WRITE(6,*)'  NS=',NS,'   N=',N,'   ALF=',ALF
      PI=4.0D0*ATAN(1.0D0)
      NDIM=2
      NF(1)=NS
      NF(2)=NS
      DO I=1,NS
         TCOS(I)=2*COS(PI*(I-1)*2.0D0/NS)*ALF*NSSQ
      END DO
      NSF=NS
      NSSQF=NSSQ
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
C --- FOR THE USE OF  DC_DECSOL_2D.F
        IJAC=1
        MLJAC=0
        MUJAC=0
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
        MLMAS=0
        MUMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
        ANS=NS
        DO J=1,NS
            YY=(J-1)/ANS
            DO I=1,NS
               Y((J-1)*NS+I)=22.D0*YY*(1.D0-YY)**(1.5D0)
            END DO
        END DO
        DO I=1,NS
            XX=(I-1)/ANS
            DO J=1,NS
               Y((J-1)*NS+I+NSSQ)=27.D0*XX*(1.D0-XX)**(1.5D0)
            END DO
        END DO
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-4 
C --- SET DEFAULT VALUES 
        DO I=1,20
           WORK(I)=0.D0
           IWORK(I)=0
           ISTAT(I)=0
        END DO
C --- ENDPOINT OF INTEGRATION
        XEND=1.5D0
        IFCN=0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,2
C --- CALL OF THE SUBROUTINE RADAU5
        CALL SEULEX(N,FBRUS,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBRUS,IJAC,MLJAC,MUJAC,
     &                  FBRUS,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        DO J=1,N,267
           WRITE (8,*) Y(J)
        END DO
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND+10.D0
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
        IF (TARRAY(1).GT.5000.) STOP
c           write (6,*) y(1),y(n/2),y(n)
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
           WRITE (6,99) X,Y(1),Y(N/2),Y(N),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
