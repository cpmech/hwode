C
        SUBROUTINE FTIGE(NN,T,TH,F,RPAR,IPAR)
C ---RIGHT-HAND SIDE OF THE BEAM EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION F(NN),TH(150),U(150),V(150),W(150)
        DIMENSION ALPHA(150),BETA(150),STH(150),CTH(150)
        COMMON/NNNN/N,NNCOM,NSQ,NQUATR,DELTAS
C ----- CALCUL DES TH(I) ET DES SIN ET COS -------------
        DO 22 I=2,N
        THDIFF=TH(I)-TH(I-1)
        STH(I)=DSIN(THDIFF)
   22   CTH(I)=DCOS(THDIFF)
C -------- CALCUL DU COTE DROIT DU SYSTEME LINEAIRE -----
        IF(T.GT.3.14159265358979324D0)THEN
C --------- CASE T GREATER PI ------------
C ---------- I=1 ------------
            TERM1=(-3.D0*TH(1)+TH(2))*NQUATR
            V(1)=TERM1
C -------- I=2,..,N-1 -----------
            DO 32 I=2,N-1
            TERM1=(TH(I-1)-2.D0*TH(I)+TH(I+1))*NQUATR
  32        V(I)=TERM1
C ----------- I=N -------------
            TERM1=(TH(N-1)-TH(N))*NQUATR
            V(N)=TERM1
        ELSE
C --------- CASE T LESS EQUAL PI ------------
            FABS=1.5D0*DSIN(T)*DSIN(T)
            FX=-FABS
            FY= FABS
C ---------- I=1 ------------
            TERM1=(-3.D0*TH(1)+TH(2))*NQUATR
            TERM2=NSQ*(FY*DCOS(TH(1))-FX*DSIN(TH(1)))
            V(1)=TERM1+TERM2
C -------- I=2,..,N-1 -----------
            DO 34 I=2,N-1
            TERM1=(TH(I-1)-2.D0*TH(I)+TH(I+1))*NQUATR
            TERM2=NSQ*(FY*DCOS(TH(I))-FX*DSIN(TH(I)))
  34        V(I)=TERM1+TERM2
C ----------- I=N -------------
            TERM1=(TH(N-1)-TH(N))*NQUATR
            TERM2=NSQ*(FY*DCOS(TH(N))-FX*DSIN(TH(N)))
            V(N)=TERM1+TERM2
        END IF
C -------- COMPUTE PRODUCT DV=W ------------
        W(1)=STH(2)*V(2)
        DO 43 I=2,N-1
   43   W(I)=-STH(I)*V(I-1)+STH(I+1)*V(I+1)
        W(N)=-STH(N)*V(N-1)
C -------- TERME 3 -----------------
        DO 435 I=1,N
  435   W(I)=W(I)+TH(N+I)*TH(N+I) 
C ------- SOLVE SYSTEM CW=W ---------
        ALPHA(1)=1.D0
        DO 44 I=2,N
        ALPHA(I)=2.D0
   44   BETA(I-1)=-CTH(I)
        ALPHA(N)=3.D0
        DO 45 I=N-1,1,-1
        Q=BETA(I)/ALPHA(I+1)
        W(I)=W(I)-W(I+1)*Q
   45   ALPHA(I)=ALPHA(I)-BETA(I)*Q
        W(1)=W(1)/ALPHA(1)
        DO 46 I=2,N
   46   W(I)=(W(I)-BETA(I-1)*W(I-1))/ALPHA(I)
C -------- COMPUTE U=CV+DW ---------
        U(1)=V(1)-CTH(2)*V(2)+STH(2)*W(2)
        DO 47 I=2,N-1
   47   U(I)=2.D0*V(I)-CTH(I)*V(I-1)-CTH(I+1)*V(I+1)
     &                -STH(I)*W(I-1)+STH(I+1)*W(I+1)
        U(N)=3.D0*V(N)-CTH(N)*V(N-1)-STH(N)*W(N-1)
C -------- PUT  DERIVATIVES IN RIGHT PLACE -------------
        DO 54 I=1,N
        F(I)=TH(N+I)
  54    F(N+I)=U(I)
        RETURN
        END
