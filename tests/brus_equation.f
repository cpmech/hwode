C
        SUBROUTINE FBRUS(NNN,X,Y,F,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(NNN),F(NNN)
        COMMON/PARAM/N,N2,GAMMA,GAMMA2
        I=1
        IU=2*I-1
        IV=2*I
        UI=Y(IU)
        VI=Y(IV)
            UIM=1.D0
            VIM=3.D0
            UIP=Y(IU+2)
            VIP=Y(IV+2)
        PROD=UI*UI*VI
        F(IU)=1.D0+PROD-4.D0*UI+GAMMA*(UIM-2.D0*UI+UIP)
        F(IV)=3.D0*UI-PROD+GAMMA*(VIM-2.D0*VI+VIP)
        DO 5 I=2,N-1
        IU=2*I-1
        IV=2*I
        UI=Y(IU)
        VI=Y(IV)
            UIM=Y(IU-2)
            VIM=Y(IV-2)
            UIP=Y(IU+2)
            VIP=Y(IV+2)
        PROD=UI*UI*VI
        F(IU)=1.D0+PROD-4.D0*UI+GAMMA*(UIM-2.D0*UI+UIP)
        F(IV)=3.D0*UI-PROD+GAMMA*(VIM-2.D0*VI+VIP)
    5   CONTINUE
        I=N
        IU=2*I-1
        IV=2*I
        UI=Y(IU)
        VI=Y(IV)
            UIM=Y(IU-2)
            VIM=Y(IV-2)
            UIP=1.D0
            VIP=3.D0
        PROD=UI*UI*VI
        F(IU)=1.D0+PROD-4.D0*UI+GAMMA*(UIM-2.D0*UI+UIP)
        F(IV)=3.D0*UI-PROD+GAMMA*(VIM-2.D0*VI+VIP)
        RETURN
        END
C
        SUBROUTINE JBRUS(NN,X,Y,DFY,LDFY,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(NN),DFY(LDFY,NN)
        COMMON/PARAM/N,N2,GAMMA,GAMMA2
        DO 1 I=1,N
        IU=2*I-1
        IV=2*I
        UI=Y(IU)
        VI=Y(IV)
        UIVI=UI*VI
        UI2=UI*UI
        DFY(3,IU)=2.D0*UIVI-4.D0-GAMMA2
        DFY(2,IV)=UI2
        DFY(4,IU)=3.D0-2.D0*UIVI
        DFY(3,IV)=-UI2-GAMMA2
        DFY(2,IU)=0.D0
        DFY(4,IV)=0.D0
    1   CONTINUE
        DO 2 I=1,N2-2
        DFY(1,I+2)=GAMMA
        DFY(5,I)=GAMMA
    2   CONTINUE
        RETURN
        END
