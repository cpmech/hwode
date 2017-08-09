        SUBROUTINE FVANDER(N,X,Y,F,EPS,IPAR)
C --- RIGHT-HAND SIDE OF VANDERPOL EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=Y(2)
        PROD=1.D0-Y(1)*Y(1)
        F(2)=(PROD*Y(2)-Y(1))/EPS
        RETURN
        END 
C
        SUBROUTINE JVANDER(N,X,Y,DFY,LDFY,EPS,IPAR)
C --- JACOBIAN OF VANDERPOL EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=0.D0
        DFY(1,2)=1.D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/EPS
        DFY(2,2)=(1.0D0-Y(1)**2)/EPS
        RETURN
        END

