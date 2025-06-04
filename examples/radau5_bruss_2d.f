C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON BRUSS-2D
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (NSD=36,ND=NSD*NSD*2,LJAC=1,LMAS=0,LE=1)
        PARAMETER (LWORK=ND*(LJAC+LMAS+3*LE+12)+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2),TRESULT
        EXTERNAL FBRUS,JBRUSF,SOLOUT
        LOGICAL DEBUG
        COMMON/TRANS/ALF,NS,NSSQ,NSNSM1,NSM1SQ
        COMMON /FOURIER/TCOS(512),NF(2),ALPH,NDIM,NSF,NSSQF
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
        REWIND 8
C ----- DIMENSIONS --------
      NS=36
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
        NRLOOP=1
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
        XEND=0.01D0
        CALL DTIME(TARRAY,TRESULT)
        DO 20 I=1,2
C --- CALL OF THE SUBROUTINE RADAU5
        DEBUG=.FALSE.
        CALL RADAU5(N,FBRUS,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBRUSF,IJAC,MLJAC,MUJAC,
     &                  FBRUS,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
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
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(N/2),Y(N),NR-1
 99     FORMAT(1X,'X =',F5.2,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE FBRUS(N,X,Y,F,RPAR,IPAR)
C ------ BRUSSELATOR WITH DIFFUSION IN D2 ----
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        COMMON /COUNT / NFCN 
        COMMON/TRANS/ALF,NS,NSSQ,NSNSM1,NSM1SQ
        NFCN=NFCN+1
C ---- CONSTANTS FOR INHOMOGENITY ---
        ANS=NS
        RADSQ=0.1D0**2
        IF(X.GE.1.1D0)THEN
         BET=5.00D0
        ELSE
         BET=0.00D0
        END IF
C ------- GRANDE BOUCLE -----
      DO I=1,NSSQ
C ------- LEFT NEIGHBOUR ---
         IF(MOD(I,NS).EQ.1)THEN
            ULEFT=Y(I+NS-1)
            VLEFT=Y(NSSQ+I+NS-1)
         ELSE
            ULEFT=Y(I-1)
            VLEFT=Y(NSSQ+I-1)
         END IF
C ------- RIGHT NEIGHBOUR ---
         IF(MOD(I,NS).EQ.0)THEN
            URIGHT=Y(I-NS+1)
            VRIGHT=Y(NSSQ+I-NS+1)
         ELSE
            URIGHT=Y(I+1)
            VRIGHT=Y(NSSQ+I+1)
         END IF
C ------- LOWER NEIGHBOUR ---
         IF(I.LE.NS)THEN
            ULOW=Y(I+NSNSM1)
            VLOW=Y(NSSQ+I+NSNSM1)
         ELSE
            ULOW=Y(I-NS)
            VLOW=Y(NSSQ+I-NS)
         END IF
C ------- UPPER NEIGHBOUR ---
         IF(I.GT.NSNSM1)THEN
            UUP=Y(I-NSNSM1)
            VUP=Y(NSSQ+I-NSNSM1)
         ELSE
            UUP=Y(I+NS)
            VUP=Y(NSSQ+I+NS)
         END IF
C ------ LA DERIVEE --------
         UIJ=Y(I)
         VIJ=Y(I+NSSQ)
         F(I)=1.D0+UIJ*UIJ*VIJ-4.4D0*UIJ
     &        +ALF*NSSQ*(ULEFT+URIGHT+ULOW+UUP-4.D0*UIJ)
         F(I+NSSQ)=3.4D0*UIJ - UIJ*UIJ*VIJ
     &        +ALF*NSSQ*(VLEFT+VRIGHT+VLOW+VUP-4.D0*VIJ)
C ----- INHOMOGENITE ----
         IY=(I-1)/NS+1
         IX=I-(IY-1)*NS
         YY=IY/ANS
         XX=IX/ANS
         IF(((XX-0.3D0)**2+(YY-0.6D0)**2).LE.RADSQ)THEN
           F(I)=F(I)+BET
         END IF
      END DO
      RETURN
      END
C
        SUBROUTINE JBRUSF(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF BRUSS-2D (FULL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        COMMON/TRANS/ALF,NS,NSSQ,NSNSM1,NSM1SQ
        FAC=ALF*NSSQ
         DO I=1,N
            DO J=1,N
               DFY(I,J)=0.0D0
            END DO
            DFY(I,I)=-FAC*4.0D0
         END DO
         DO I=1,N
            IL=I-1
            IF (MOD(I,NS).EQ.1) IL=I-1+NS
            DFY(I,IL)=FAC
            IR=I+1
            IF (MOD(I,NS).EQ.0) IR=I+1-NS
            DFY(I,IR)=FAC
            IU=I+NS
            IF (I.GT.NS*(NS-1).AND.I.LE.NS*NS) IU=I-NS*(NS-1)
            IF (I.GT.NS*(2*NS-1)) IU=I-NS*(NS-1)
            DFY(I,IU)=FAC
            ID=I-NS
            IF (I.LE.NS) ID=I+NS*(NS-1)
            IF (I.LE.NS*(NS+1).AND.I.GT.NS*NS) ID=I+NS*(NS-1)
            DFY(I,ID)=FAC
         END DO
      RETURN
      END

