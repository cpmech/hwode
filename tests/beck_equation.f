        SUBROUTINE FBECK(N,X,Y,F,RPAR,IPAR)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),F(N)
        COMMON /COEFF/B(5000),RHO
        Y1=Y(1)
        AJO=Y(N-1)*Y1-B(N)*Y(N)
        F(N)=AJO
        SUMJ=AJO
        DO I=N-1,2,-1
           AJ=Y(I-1)*Y1-B(I)*Y(I)
           F(I)=AJ-AJO
           SUMJ=SUMJ+AJ
           AJO=AJ
        END DO
        F(1)=-SUMJ-AJ
        RETURN
        END 
C
        SUBROUTINE JBECK(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF BECKDO IN THE CASE THAT THE LINEAR ALGEBRA
C --- OF  DC_DECSOL_BD.F  IS USED
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        COMMON /COEFF/B(5000),RHO
        Y1=Y(1)
        SUM=0.0D0
        DO I=N-1,2,-1
           SUM=SUM+Y(I)
        END DO
        DFY(2,1)=-4*Y1-SUM
        DO I=3,N
           DFY(1,I)=B(I)
        END DO
        DO I=2,N-1
           DFY(2,I)=-B(I)-Y1
           DFY(3,I)=Y1
        END DO
        DFY(2,N)=-B(N)
C
        DO I=3,N-1
           DFY(4,I)=B(I)-Y1
           DFY(5,I)=Y(I-1)-Y(I)
        END DO
        DFY(4,2)=2*B(2)-Y1
        DFY(4,N)=B(N)
        DFY(5,2)=2*Y1-Y(2)
        DFY(5,N)=Y(N-1)
        RETURN
        END
C
        SUBROUTINE JBECKF(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF BECKDO (AS FULL MATRIX)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        COMMON /COEFF/B(5000),RHO
        Y1=Y(1)
        SUM=0.0D0
        DO I=N-1,2,-1
           SUM=SUM+Y(I)
        END DO
        DFY(1,1)=-4*Y1-SUM
        DO I=3,N
           DFY(I-1,I)=B(I)
        END DO
        DO I=2,N-1
           DFY(I,I)=-B(I)-Y1
           DFY(I+1,I)=Y1
        END DO
        DFY(N,N)=-B(N)
C
        DO I=3,N-1
           DFY(1,I)=B(I)-Y1
           DFY(I,1)=Y(I-1)-Y(I)
        END DO
        DFY(1,2)=2*B(2)-Y1
        DFY(1,N)=B(N)
        DFY(2,1)=2*Y1-Y(2)
        DFY(N,1)=Y(N-1)
        RETURN
        END
C
C
        SUBROUTINE FBECKD(N,X,Y,F,RPAR,IPAR)
C --- BECKDO AS DAE OF INDEX 1
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),F(N)
        COMMON /COEFF/B(5000),RHO
        Y1=Y(1)
        AJO=Y(N-1)*Y1-B(N)*Y(N)
        F(N)=AJO
        DO I=N-1,2,-1
           AJ=Y(I-1)*Y1-B(I)*Y(I)
           F(I)=AJ-AJO
           AJO=AJ
        END DO
        SUM=0.D0
        DO I=N,1,-1
           SUM=SUM+I*Y(I)
        END DO
        F(1)=SUM-RHO
        RETURN
        END 
C
        SUBROUTINE JBECKD(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF BECKDO AS DAE OF INDEX 1
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        COMMON /COEFF/B(5000),RHO
        Y1=Y(1)
        DFY(2,1)=1.0D0
        DO I=3,N
           DFY(1,I)=B(I)
        END DO
        DO I=2,N-1
           DFY(2,I)=-B(I)-Y1
           DFY(3,I)=Y1
        END DO
        DFY(2,N)=-B(N)
C
        DO I=3,N-1
           DFY(4,I)=I
           DFY(5,I)=Y(I-1)-Y(I)
        END DO
        DFY(4,2)=2
        DFY(4,N)=N
        DFY(5,2)=2*Y1-Y(2)
        DFY(5,N)=Y(N-1)
        RETURN
        END
C
        SUBROUTINE MBECK(N,AM,LMAS,RPAR,IPAR)
        DOUBLE PRECISION AM(LMAS,N)
        AM(1,1)=0.0D0
        DO I=2,N
           AM(1,I)=1.0D0
        END DO
        RETURN
        END
