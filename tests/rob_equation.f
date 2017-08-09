        SUBROUTINE FROBER(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF ROBERTSON EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=-0.04D0*Y(1)+1.D4*Y(2)*Y(3)
        F(3)=3.D7*Y(2)*Y(2)
        F(2)=-F(1)-F(3)
        RETURN
        END 
C
        SUBROUTINE JROBER(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF ROBERTSON EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        PROD1=1.0D4*Y(2)
        PROD2=1.0D4*Y(3)
        PROD3=6.0D7*Y(2)
        DFY(1,1)=-0.04D0
        DFY(1,2)=PROD2
        DFY(1,3)=PROD1
        DFY(2,1)=0.04D0
        DFY(2,2)=-PROD2-PROD3
        DFY(2,3)=-PROD1
        DFY(3,1)=0.D0
        DFY(3,2)=PROD3
        DFY(3,3)=0.D0
        RETURN
        END

