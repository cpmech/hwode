      FUNCTION HINIT(N,FCN,X,Y,XEND,POSNEG,F0,F1,Y1,IORD,
     &                 HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
C ----------------------------------------------------------
C ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)
C ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
C ----   H = 0.01 * NORM (Y0) / NORM (F0)
C ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
C ---- COMPARED TO THE SOLUTION
      DNF=0.0D0
      DNY=0.0D0 
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      IF (ITOL.EQ.0) THEN   
        DO 10 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  10    DNY=DNY+(Y(I)/SK)**2 
      ELSE
        DO 11 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  11    DNY=DNY+(Y(I)/SK)**2 
      END IF
      IF (DNF.LE.1.D-10.OR.DNY.LE.1.D-10) THEN
         H=1.0D-6
      ELSE
         H=SQRT(DNY/DNF)*0.01D0 
      END IF
      H=MIN(H,HMAX)
      H=SIGN(H,POSNEG) 
C ---- PERFORM AN EXPLICIT EULER STEP
      DO 12 I=1,N
  12  Y1(I)=Y(I)+H*F0(I)
      CALL FCN(N,X+H,Y1,F1,RPAR,IPAR) 
C ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
      DER2=0.0D0
      IF (ITOL.EQ.0) THEN   
        DO 15 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
  15    DER2=DER2+((F1(I)-F0(I))/SK)**2   
      ELSE
        DO 16 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
  16    DER2=DER2+((F1(I)-F0(I))/SK)**2   
      END IF
      DER2=SQRT(DER2)/H
C ---- STEP SIZE IS COMPUTED SUCH THAT
C ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
      DER12=MAX(ABS(DER2),SQRT(DNF))
      IF (DER12.LE.1.D-15) THEN
         H1=MAX(1.0D-6,ABS(H)*1.0D-3)
      ELSE
         H1=(0.01D0/DER12)**(1.D0/IORD) 
      END IF
      H=MIN(100*ABS(H),H1,HMAX)
      HINIT=SIGN(H,POSNEG)  
      RETURN
      END 
