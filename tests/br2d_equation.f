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
        SUBROUTINE JBRUS(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- IN THE CASE THAT  DC_DECSOL_2D.F  IS USED, NO JACOBIAN
C --- HAS TO BE SPECIFIED. THEREFORE THIS DUMMY SUBROUTINE
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
c
          return
c
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

