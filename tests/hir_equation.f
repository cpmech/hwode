C
      SUBROUTINE FHIRES (N, X, Y, DY, RPAR, IPAR)
      INTEGER N
      DOUBLE PRECISION X, Y, DY
      DIMENSION Y (8), DY (8)
      DY (1) = -1.71D0*Y(1) + 0.43D0*Y(2) + 8.32D0*Y(3) + 0.0007D0
      DY (2) = 1.71D0*Y(1) - 8.75D0*Y(2)
      DY (3) = -10.03D0*Y(3) + 0.43D0*Y(4) + 0.035D0*Y(5)
      DY (4) = 8.32D0*Y(2) + 1.71D0*Y(3) - 1.12D0*Y(4)
      DY (5) = -1.745D0*Y(5) + 0.43D0*Y(6) + 0.43D0*Y(7)
      DY (6) = -280.0D0*Y(6)*Y(8) + 0.69D0*Y(4) + 1.71D0*Y(5) -
     &         0.43D0*Y(6) + 0.69D0*Y(7)
      DY (7) = 280.0D0*Y(6)*Y(8) - 1.81D0*Y(7)
      DY (8) = -DY (7)
      RETURN 
      END


      SUBROUTINE JHIRES(N,X,Y,DFY,LDFY,RPAR,IPAR)
C -------- JACOBIAN FOR HIRES PROBLEM --------
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
C------ METTRE A ZERO -------
      DO 1 I=1,N
      DO 1 J=1,N
  1   DFY(I,J)=0.D0
C
      DFY(1,1)=  -1.71D0
      DFY(1,2)=   0.43D0
      DFY(1,3)=   + 8.32D0
C
      DFY(2,1)=   1.71D0
      DFY(2,2)=   - 8.75D0
C
      DFY(3,3)=   -10.03D0
      DFY(3,4)=   0.43D0
      DFY(3,5)=   + 0.035D0
C
      DFY(4,2)=    8.32D0
      DFY(4,3)=   + 1.71D0
      DFY(4,4)=   - 1.12D0
C
      DFY(5,5)=   -1.745D0
      DFY(5,6)=   + 0.43D0
      DFY(5,7)=   + 0.43D0
C
      DFY(6,4)=    + 0.69D0
      DFY(6,5)=    + 1.71D0
      DFY(6,6)=    - 0.43D0  -280.0D0*Y(8)
      DFY(6,7)=    + 0.69D0
      DFY(6,8)=     -280.0D0*Y(6)
C 
      DFY(7,6)=     280.0D0*Y(8)
      DFY(7,7)=      - 1.81D0
      DFY(7,8)=     280.0D0*Y(6)
C
      DFY(8,6)=    - 280.0D0*Y(8)
      DFY(8,7)=       1.81D0
      DFY(8,8)=    - 280.0D0*Y(6)
      RETURN
      END


