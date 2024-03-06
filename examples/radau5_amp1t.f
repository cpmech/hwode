C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT THE ONE-TRANSISTOR AMPLIFIER PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (BANDED JACOBIAN)
        PARAMETER (ND=5,LJAC=4,LMAS=3,LE=6)
        PARAMETER (LWORK=ND*(LJAC+LMAS+3*LE+12)+20,LIWORK=3*ND+20)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),RPAR(7)
        EXTERNAL FAMPL,JBAMPL,BBAMPL,SOLOUT
        LOGICAL DEBUG
C --- DATA FOR THE PROBLEM
        UE=0.4D0
          RPAR(1)=UE
        UB=6.0D0
          RPAR(2)=UB
        UF=0.026D0
          RPAR(3)=UF
        ALPHA=0.99D0
          RPAR(4)=ALPHA
        BETA=1.0D-6
          RPAR(5)=BETA
        RR=1000.0D0
          RPAR(6)=RR
        SS=9000.0D0
          RPAR(7)=SS
C --- DIMENSION OF THE SYSTEM
        N=5
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A BANDED MATRIX (LOWER AND UPPER BAND-WIDTHS ARE 2, 1)
        MLJAC=2
        MUJAC=1
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
        IMAS=1
        MLMAS=1
        MUMAS=1
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=0.0D0
        Y(2)=UB/2.0D0
        Y(3)=UB/2.0D0
        Y(4)=UB
        Y(5)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=0.05D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-4
        ATOL=1.0D-4
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
C --- CALL OF THE SUBROUTINE RADAU5
        write(*,'(A)')''
        write(*,'(A)')'running radau5.f test'
        DEBUG=.FALSE.
        CALL RADAU5(N,FAMPL,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBAMPL,IJAC,MLJAC,MUJAC,
     &                  BBAMPL,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,DEBUG,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
      write(*,'(A)',advance='no')'Radau5: Radau method (Radau IIA) '
      write(*,'(A)')'(implicit, order 5, embedded)'
      write(*,'(A,I0)')'Number of function evaluations   = ',1+IWORK(14)
      write(*,'(A,I0)')'Number of Jacobian evaluations   = ',IWORK(15)
      write(*,'(A,I0)')'Number of factorizations         = ',IWORK(19)
      write(*,'(A,I0)')'Number of lin sys solutions      = ',IWORK(20)
      write(*,'(A,I0)')'Number of performed steps        = ',IWORK(16)
      write(*,'(A,I0)')'Number of accepted steps         = ',IWORK(17)
      write(*,'(A,I0)')'Number of rejected steps         = ',IWORK(18)
      write(*,'(A,I0)')'Number of iterations (maximum)   = ',IWORK(21)
      write(*,'(A,3ES23.15)')'y1to3 =',Y(1),Y(2),Y(3)
      write(*,'(A,2ES23.15)')'y4to5 =',Y(4),Y(5)
      write(*,'(A,ES23.15)')'h =',H
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTR5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) NR-1,X,Y(1),Y(5)
           XOUT=0.001D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C --- CONTINUOUS OUTPUT FOR RADAU5
              WRITE (6,99) NR-1,XOUT,CONTR5(1,XOUT,CONT,LRC),
     &                     CONTR5(5,XOUT,CONT,LRC)
              XOUT=XOUT+0.001D0
              GOTO 10
           END IF
        END IF
 99     FORMAT('step =',I4,', x =',F7.4,', y1and5 =',2ES23.15)
        RETURN
        END
C
C
        SUBROUTINE FAMPL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 Y(N),F(N),RPAR(7)
        UE=RPAR(1)
        UB=RPAR(2)
        UF=RPAR(3)
        ALPHA=RPAR(4)
        BETA=RPAR(5)
        RR=RPAR(6)
        SS=RPAR(7)
        GAMMA=1.0D0-ALPHA
        PI=3.141592653589793238462643383279502884197169399375105
     &     82097494459230781640628620899862803482534211706798214D0;
        UET=UE*DSIN(200.D0*PI*X)
        F23=BETA*(DEXP((Y(2)-Y(3))/UF)-1.D0)
        F(1)=(Y(1)-UET)/RR
        F(2)=(2.0D0*Y(2)-UB)/SS+GAMMA*F23
        F(3)=Y(3)/SS-F23
        F(4)=(Y(4)-UB)/SS+ALPHA*F23
        F(5)=Y(5)/SS
        RETURN
        END
C
C
        SUBROUTINE JBAMPL(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 Y(N),DFY(LDFY,N),RPAR(7)
        UE=RPAR(1)
        UB=RPAR(2)
        UF=RPAR(3)
        ALPHA=RPAR(4)
        BETA=RPAR(5)
        RR=RPAR(6)
        SS=RPAR(7)
        GAMMA=1.0D0-ALPHA
        G23=BETA*DEXP((Y(2)-Y(3))/UF)/UF
        DFY(1,1)=0.0D0
        DFY(1,2)=0.0D0
        DFY(1,3)=-GAMMA*G23
        DFY(1,4)=0.0D0
        DFY(1,5)=0.0D0
        DFY(2,1)=1.0D0/RR
        DFY(2,2)=2.0D0/SS+GAMMA*G23
        DFY(2,3)=1.0D0/SS+G23
        DFY(2,4)=1.0D0/SS
        DFY(2,5)=1.0D0/SS
        DFY(3,1)=0.0D0
        DFY(3,2)=-G23
        DFY(3,3)=-ALPHA*G23
        DFY(3,4)=0.0D0
        DFY(3,5)=0.0D0
        DFY(4,1)=0.0D0
        DFY(4,2)=ALPHA*G23
        DFY(4,3)=0.0D0
        DFY(4,4)=0.0D0
        DFY(4,5)=0.0D0
        RETURN
        END 
C
        SUBROUTINE BBAMPL(N,B,LB,RPAR,IPAR)
C --- MATRIX "M" FOR THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 B(LB,N),RPAR(7)
        C1=1.0D-6
        C2=2.0D-6
        C3=3.0D-6
C
        B(1,1)=0.0D0
        B(1,2)=C1
        B(1,3)=0.0D0
        B(1,4)=0.0D0
        B(1,5)=C3
        B(2,1)=-C1
        B(2,2)=-C1
        B(2,3)=-C2
        B(2,4)=-C3
        B(2,5)=-C3
        B(3,1)=C1
        B(3,2)=0.0D0
        B(3,3)=0.0D0
        B(3,4)=C3
        B(3,5)=0.0D0
        RETURN
        END



