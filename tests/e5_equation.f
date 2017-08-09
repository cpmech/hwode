        SUBROUTINE FE5(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF E5 PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        PROD1=7.89D-10*Y(1)
        PROD2=1.1D7*Y(1)*Y(3)
        PROD3=1.13D9*Y(2)*Y(3)
        PROD4=1.13D3*Y(4)
        F(1)=-PROD1-PROD2
        F(2)=PROD1-PROD3
        F(4)=PROD2-PROD4
        F(3)=F(2)-F(4)
        RETURN
        END 
C
        SUBROUTINE JE5(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF E5 PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        A=7.89D-10
        B=1.1D7
        CM=1.13D9
        C=1.13D3
        DFY(1,1)=-A-B*Y(3)
        DFY(1,2)=0.D0
        DFY(1,3)=-B*Y(1)
        DFY(1,4)=0.D0
        DFY(2,1)=A
        DFY(2,2)=-CM*Y(3)
        DFY(2,3)=-CM*Y(2)
        DFY(2,4)=0.D0
        DFY(3,1)=A-B*Y(3)
        DFY(3,2)=-CM*Y(3)
        DFY(3,3)=-B*Y(1)-CM*Y(2)
        DFY(3,4)=C
        DFY(4,1)=B*Y(3)
        DFY(4,2)=0.D0
        DFY(4,3)=B*Y(1)
        DFY(4,4)=-C
        RETURN
        END

