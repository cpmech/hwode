C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT HW-EQ1 EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
C link dr1_radau5 radau5 decsol dc_decsol   or
C link dr1_radau5 radau5 lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=1,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FCN,JAC,SOLOUT
C --- PARAMETER IN THE DIFFERENTIAL EQUATION
        RPAR=1.0D-6
C --- DIMENSION OF THE SYSTEM
        N=1
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=1.5D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-4
        ATOL=1.0D0*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-4 
C --- SET DEFAULT VALUES 
        DO I=1,20
           IWORK(I)=0
           WORK(I)=0.D0
        END DO
C --- CALL OF THE SUBROUTINE RADAU5
        CALL RADAU5(N,FCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC,IJAC,MLJAC,MUJAC,
     &                  FCN,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1)
 99     FORMAT(1X,'X =',F5.2,'    Y =',1E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) RTOL
 90     FORMAT('       rtol=',D8.2)
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,' accpt=',I4,
     &        ' rejct=',I3,' dec=',I4,' sol=',I5)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),NR-1
           XOUT=0.2D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C --- CONTINUOUS OUTPUT FOR RADAU5
              WRITE (6,99) XOUT,CONTR5(1,XOUT,CONT,LRC),NR-1
              XOUT=XOUT+0.2D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',1E18.10,'    NSTEP =',I4)
        RETURN
        END
C
C
        SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF HW-EQ1
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=RPAR*(Y(1)-DCOS(X))
        RETURN
        END 
C
C
        SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF HW-EQ1
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=RPAR
        RETURN
        END


