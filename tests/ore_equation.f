        SUBROUTINE FOREGON(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF OREGO EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=77.27D0*(Y(2)+Y(1)*(1.D0-8.375D-6*Y(1)-Y(2)))
        F(2)=(Y(3)-(1.D0+Y(1))*Y(2))/77.27D0
        F(3)=0.161D0*(Y(1)-Y(3))
        RETURN
        END 
C
        SUBROUTINE JOREGON(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF OREGO EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=77.27D0*(1.D0-2.D0*8.375D-6*Y(1)-Y(2))
        DFY(1,2)=77.27D0*(1.D0-Y(1))
        DFY(1,3)=0.D0
        DFY(2,1)=-Y(2)/77.27D0
        DFY(2,2)=-(1.D0+Y(1))/77.27D0
        DFY(2,3)=1.D0/77.27D0
        DFY(3,1)=.161D0
        DFY(3,2)=.0D0
        DFY(3,3)=-.161D0
        RETURN
        END

