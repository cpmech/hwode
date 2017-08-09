        SUBROUTINE FCUSP(N,T,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF THE CUSP EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        COMMON/NERVES/NNERV
        COMMON/DIFFCOEF/DIFFUS
c-----------
        DO 25 INERV=1,NNERV
        X=Y(3*INERV-2)
        A=Y(3*INERV-1)
        B=Y(3*INERV)
        IF(INERV.EQ.1)THEN
            XRIGHT=Y(3*NNERV-2)
            ARIGHT=Y(3*NNERV-1)
            BRIGHT=Y(3*NNERV)
        ELSE
            XRIGHT=Y(3*INERV-5)
            ARIGHT=Y(3*INERV-4)
            BRIGHT=Y(3*INERV-3)
        END IF
        IF(INERV.EQ.NNERV)THEN
            XLEFT=Y(1)
            ALEFT=Y(2)
            BLEFT=Y(3)
        ELSE
            XLEFT=Y(3*INERV+1)
            ALEFT=Y(3*INERV+2)
            BLEFT=Y(3*INERV+3)
        END IF  
        XDOT=-10000.D0*(B+X*(A+X*X))
        U=(X-0.7D0)*(X-1.3D0)
        V=U/(U+0.1D0)
        ADOT=B+0.07D0*V
        BDOT=(1.D0*(1.D0-A*A)*B-A)-0.4D0*X+0.035D0*V
        F(3*INERV-2)=XDOT+DIFFUS*(XLEFT-2.D0*X+XRIGHT)
        F(3*INERV-1)=ADOT+DIFFUS*(ALEFT-2.D0*A+ARIGHT)
        F(3*INERV)  =BDOT+DIFFUS*(BLEFT-2.D0*B+BRIGHT)
  25    CONTINUE
        RETURN
        END

