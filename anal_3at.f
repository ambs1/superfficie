      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=3,NT=3*NA,NTES=5000,PREC=1.0D-4)
      PARAMETER(MA=112,ECM1=219474.62492D0,ATMAS=1.82288853D3)
      DIMENSION X(NT),R(NA*(NA-1)/2)
      DIMENSION SOL(NA*(NA-1)/2,NTES),ENER(NTES),
     &  N1(NTES),N2(NTES),VAL(NA*(NA-1)/2),XSOL(NT,NTES),
     &  VALP(NT),VALPT(NT,NTES),CINT(NA*(NA-1)/2),VECP(NT,NT),
     &  VECPT(NT,NT,NTES),RMM(NA),FREQ(NA*(NA-1)/2+1)
      CHARACTER*1 IFR(4)
      CHARACTER*2 ATOM(NA)
      OPEN(5,FILE='min_3at.dat',STATUS='old')
      OPEN(43,FILE='min_3at.res',STATUS='unknown')
      OPEN(44,FILE='min_3at.tex',STATUS='unknown')
      OPEN(45,FILE='min_3at.tab',STATUS='unknown')
      OPEN(98,FILE='min_3at.molden',STATUS='unknown')
      RNA=REAL(NA)
      KTES=1
      DO I=1,NA
      READ(5,*)ATOM(I),RMM(I)
      ENDDO
      READ(5,*)ECUT,DRMAX,NTEST,IREAD,IWRIT

      IF(IREAD.EQ.0)THEN
      NTEST=1
      READ(5,*)(X(I),I=1,NT)
         CALL TRANS(X,R)          
         E=POTEN(R)
         IF(IWRIT.EQ.1)THEN
         WRITE(6,'(12F12.4)')R,E
         ENDIF
         CALL NR(X,IND,VALP,VECP,RMM,ECUT)
         IF(IND.EQ.100)STOP ' USE A NEW ESTIMATE...'  
         CALL TRANS(X,R)          
         E=POTEN(R)
               DO KTR=1,NA*(NA-1)/2
                   SOL(KTR,1)=R(KTR)
               ENDDO
                NN=0
                NP=0
               DO KTR=1,NT
                 XSOL(KTR,1)=X(KTR)
                 VALPT(KTR,1)=VALP(KTR)
                 DO LTR=1,NT
                  VECPT(KTR,LTR,1)=VECP(KTR,LTR)
                 ENDDO
                 IF(ABS(VALP(KTR)).GT.PREC)THEN
                    IF(VALP(KTR).LT.0.0)THEN
                      NN=NN+1
                    ELSE
                     NP=NP+1
                    ENDIF
                  ENDIF
                 ENDDO
           ENER(KTES)=E
         N1(1)=NN
         N2(1)=NP

      ELSE IF(IREAD.EQ.1)THEN

      READ(5,*)(X(I),I=1,NT)
         CALL TRANS(X,R)          
         E=POTEN(R)
         IF(IWRIT.EQ.1)THEN
         WRITE(6,'(12F12.4)')R,E
         ENDIF
         CALL NR(X,IND,VALP,VECP,RMM,ECUT)
         IF(IND.EQ.100)STOP ' USE A NEW ESTIMATE...'  
         CALL TRANS(X,R)          
         E=POTEN(R)
               DO KTR=1,NA*(NA-1)/2
                   SOL(KTR,1)=R(KTR)
               ENDDO
                NN=0
                NP=0
               DO KTR=1,NT
                 XSOL(KTR,1)=X(KTR)
                 VALPT(KTR,1)=VALP(KTR)
                 DO LTR=1,NT
                  VECPT(KTR,LTR,1)=VECP(KTR,LTR)
                 ENDDO
                 IF(ABS(VALP(KTR)).GT.PREC)THEN
                    IF(VALP(KTR).LT.0.0)THEN
                      NN=NN+1
                    ELSE
                     NP=NP+1
                    ENDIF
                  ENDIF
                 ENDDO
           ENER(KTES)=E
         N1(1)=NN
         N2(1)=NP

      DO K=1,NTEST
         DO I=1,NT
              X(I)=DRMAX*(1.0-2.0*RANF(DUM))
         ENDDO
         CALL TRANS(X,R)
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(I5,12F12.4)')K,R,E
         ENDIF
         CALL NR(X,IND,VALP,VECP,RMM,ECUT)
         IF(IND.EQ.100)GOTO 1213
         CALL TRANS(X,R)
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(I5,12F12.4)')K,R,E
         ENDIF
        DO KT=1,KTES
              NNT=0
           DO KTR=1,NA*(NA-1)/2
             IF(ABS(R(KTR)-SOL(KTR,KT)).LT.2*PREC)NNT=NNT+1
           ENDDO
           IF(NNT.EQ.(NA*(NA-1)/2))GOTO 11
        ENDDO
        KTES=KTES+1
             DO KTR=1,NA*(NA-1)/2
              SOL(KTR,KTES)=R(KTR)
            ENDDO
            NN=0
            NP=0
             DO KTR=1,NT
              VALPT(KTR,KTES)=VALP(KTR)
               DO LTR=1,NT
                 VECPT(KTR,LTR,KTES)=VECP(KTR,LTR)
               ENDDO
              XSOL(KTR,KTES)=X(KTR)
              IF(ABS(VALP(KTR)).GT.PREC)THEN
                 IF(VALP(KTR).LT.0.0)THEN
                    NN=NN+1
                  ELSE
                   NP=NP+1
                  ENDIF
               ENDIF
             ENDDO
       ENER(KTES)=E
       N1(KTES)=NN
       N2(KTES)=NP
 11    CONTINUE
 1213  CONTINUE
       ENDDO

      ELSE

       IF(IWRIT.GT.0)THEN
       WRITE(6,*)' FIRST STATONARY POINT'
       ENDIF
      NT0=0
 1211 CONTINUE
      NT0=NT0+1
      IF(NT0.GT.NTES)STOP 'PLEASE TRY TO CHANGE DRMAX'
      DO I=1,NT
              X(I)=DRMAX*(1.0-2.0*RANF(DUM))
            ENDDO
         CALL TRANS(X,R)          
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(I5,12F12.4)')NT0,R,E
         ENDIF
         CALL NR(X,IND,VALP,VECP,RMM,ECUT)
         IF(IND.EQ.100)GOTO 1211 
         CALL TRANS(X,R)          
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(A5,12F12.4)')'SOLUC',R,E
         ENDIF
               DO KTR=1,NA*(NA-1)/2
                   SOL(KTR,1)=R(KTR)
               ENDDO
                NN=0
                NP=0
               DO KTR=1,NT
                 XSOL(KTR,1)=X(KTR)
                 VALPT(KTR,1)=VALP(KTR)
                 DO LTR=1,NT
                  VECPT(KTR,LTR,1)=VECP(KTR,LTR)
                 ENDDO
                 IF(ABS(VALP(KTR)).GT.PREC)THEN
                    IF(VALP(KTR).LT.0.0)THEN
                      NN=NN+1
                    ELSE
                     NP=NP+1
                    ENDIF
                  ENDIF
                 ENDDO
            ENER(KTES)=E
            N1(1)=NN
            N2(1)=NP

      DO K=1,NTEST
         DO I=1,NT
              X(I)=DRMAX*(1.0-2.0*RANF(DUM))
         ENDDO
         CALL TRANS(X,R)
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(I5,12F12.4)')K,R,E
         ENDIF
         CALL NR(X,IND,VALP,VECP,RMM,ECUT)
         IF(IND.EQ.100)GOTO 1212
         CALL TRANS(X,R)
         E=POTEN(R)
         IF(IWRIT.GT.0)THEN
         WRITE(6,'(I5,12F12.4)')K,R,E
         ENDIF
        DO KT=1,KTES
              NNT=0
           DO KTR=1,NA*(NA-1)/2
             IF(ABS(R(KTR)-SOL(KTR,KT)).LT.2*PREC)NNT=NNT+1
           ENDDO
           IF(NNT.EQ.(NA*(NA-1)/2))GOTO 10
        ENDDO
        KTES=KTES+1
             DO KTR=1,NA*(NA-1)/2
              SOL(KTR,KTES)=R(KTR)
            ENDDO
            NN=0
            NP=0
             DO KTR=1,NT
              VALPT(KTR,KTES)=VALP(KTR)
               DO LTR=1,NT
                 VECPT(KTR,LTR,KTES)=VECP(KTR,LTR)
               ENDDO
              XSOL(KTR,KTES)=X(KTR)
              IF(ABS(VALP(KTR)).GT.PREC)THEN
                 IF(VALP(KTR).LT.0.0)THEN
                    NN=NN+1
                  ELSE
                   NP=NP+1
                  ENDIF
               ENDIF
             ENDDO
       ENER(KTES)=E
       N1(KTES)=NN
       N2(KTES)=NP
 10    CONTINUE
1212   CONTINUE
       ENDDO

      ENDIF

      DO I=1,KTES
      WRITE(98,*)'[Molden Format]'
      WRITE(98,*)'               '
      WRITE(98,*)'[FREQ]'
      DO K=1,NT
      WRITE(98,*)SQRT(ABS(VALPT(K,I)/ATMAS))*ECM1*
     &VALPT(K,I)/sqrt(VALPT(K,I)**2)
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'[FR-COORD]'
      DO IN=1,NA
      NK=3*(IN-1)
      WRITE(98,*)ATOM(IN),(XSOL(K,I),K=NK+1,NK+3)
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'[FR-NORM-COORD]'
      DO J=1,NT
      WRITE(98,*)'vibration',J
      WRITE(98,'(3F20.7)')(VECPT(K,J,I),K=1,NT)
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'               '
      WRITE(98,*)'               '
      WRITE(98,*)'               '
        WRITE(43,*)'  ******************************** '
        WRITE(43,1005)ENER(I),N1(I),N2(I)
        WRITE(43,*)' INTERNUCLEAR DISTANCES:'        
        WRITE(43,1000)(SOL(K,I),K=1,NA*(NA-1)/2)
        WRITE(43,*)' HESSIAN EIGENVALUES:'        
        WRITE(43,1000)(VALPT(K,I),K=1,NT)
        WRITE(43,*)' FREQUENCIES '
        NVIB=0
        PRECM=0.1
        DO K=1,NT
        FTEST=SQRT(ABS(VALPT(K,I)/ATMAS))*ECM1
        DO L=2,K
        IF(ABS(FTEST-FREQ(L)).LT.PRECM)GOTO 134
        ENDDO
        IF(FTEST.GT.PREC.AND.VALPT(K,I).GT.0.0)THEN
        NVIB=NVIB+1
        FREQ(NVIB)=FTEST
        IFR(NVIB)=' '
        ELSE IF(FTEST.GT.PRECM.AND.VALPT(K,I).LT.0.0)THEN
        NVIB=NVIB+1
        FREQ(NVIB)=FTEST
        IFR(NVIB)='i'
        ENDIF
 134    CONTINUE
        ENDDO       
C **
        WRITE(44,2005)(SOL(K,I),K=1,NA*(NA-1)/2),ENER(I),
     &       (NINT(FREQ(K)),IFR(K),K=1,NT-6)
        WRITE(45,3005)(SOL(K,I),K=1,NA*(NA-1)/2),ENER(I),
     &       (NINT(FREQ(K)),IFR(K),K=1,NT-6)
 2005   FORMAT(2X,3(F7.3,'&'),F10.4,4(I6,A1,'&'),'\\')
 3005   FORMAT(2X,3(F7.3),F10.4,4(I6,A1))
C **
        WRITE(43,1011)(SQRT(ABS(VALPT(K,I)/ATMAS))*ECM1,K=1,NT)
        WRITE(43,*)' HESSIAN EIGENVECTORS:'
        DO K=1,NT        
        WRITE(43,1000)(VECPT(K,J,I),J=1,NT)
        ENDDO
        DO KTR=1,NT
         X(KTR)=XSOL(KTR,I)
        ENDDO
        CALL DERR(X,VAL,CINT)
        WRITE(43,*)' INTERNAL COORDINATES:'        
        WRITE(43,1000)(CINT(II),II=1,NT-6)
        WRITE(43,*)' CARTESIAN COORDINATES:'        
        WRITE(43,1000)(X(K),K=1,NT)
        WRITE(43,*)' FIRST DERIVATIVE INTERNAL COORDINATES:'        
        WRITE(43,1000)(VAL(K),K=1,NA*(NA-1)/2)
      ENDDO
      DO I=1,KTES
        IF(N1(I).GT.0)THEN
        WRITE(78,1000)(XSOL(K,I),K=1,NT)
        NNEG=0
        DO K=1,NT
         NNEG=NNEG+1
         IF(VALPT(K,I).LT.-1.0D-4)GOTO 124
        ENDDO
 124    WRITE(78,*)NNEG
        DO K=1,NT        
         WRITE(78,1000)(VECPT(K,J,I),J=1,NT)
        ENDDO
      ENDIF
      ENDDO
 1000 FORMAT(12F10.4)
 1011 FORMAT(12F10.2)
 1005 FORMAT(2X,'ENER=',F10.4,3X,'NEG=',I5,3X,'POS=',I5)
      END
 
      SUBROUTINE NR(X,IND,VALP,VV,RMM,ECUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=3,NT=3*NA,ND=NA*(NA-1)/2)
      DIMENSION X(NT),R(ND),P(NT,NT),VAL(ND),VV(NT,NT),RMAT(NT,NT+1),
     &  VALP(NT),V(NT,NT),CINT(ND),AUX1(NT,NT),VDER(NT),
     &  RAZ(NT),RM(NT,NT),RMM(NA)
      IND=0
      PREC=1.0D-4
      NITER=1000
      CALL TRANS(X,R)
      E=POTEN(R)
      IF(E.GT.ECUT)THEN
      IND=100
      RETURN
      ENDIF
      CALL DERR(X,VAL,CINT)
      CALL PROJECT(X,P,IND)
      IF(IND.EQ.100)RETURN
      CALL DDER(X,V)
      CALL DER(X,VDER)
      CALL MATP(P,V,AUX1,NT,NT,NT,NT)
      CALL MATP(AUX1,P,V,NT,NT,NT,NT)      
      CALL JACOBI(V,NT,NT,VALP,VV,NROT)
      NN=0
      DO I=1,NT
      IF(ABS(VALP(I)).GT.PREC.AND.VALP(I).LT.0.0)NN=NN+1
      ENDDO
      IF(NN.GT.1)THEN
      IND=100
      RETURN
      ENDIF
      DO KK=1,NITER
       CALL TRANS(X,R)
       DO I=1,ND
         IF(R(I).GT.12.0)THEN
         IND=100
         RETURN
         ENDIF
        ENDDO
       DO I=1,NT
        RMAT(I,NT+1)=-VDER(I)
        DO J=1,NT
         RMAT(I,J)=V(I,J)
        ENDDO
       ENDDO
       CALL SVD(RMAT,NT,RAZ,IND)
       IF(IND.EQ.100)RETURN
       NZ1=0
       NZ2=0
       DO I=1,NT
        IF(ABS(RAZ(I)).LT.PREC)NZ1=NZ1+1
        IF(ABS(VDER(I)).LT.PREC)NZ2=NZ2+1
        X(I)=X(I)+RAZ(I)
       ENDDO
       IF(NZ1.EQ.NT.AND.NZ2.EQ.NT)GOTO 20
       CALL DERR(X,VAL,CINT)
       CALL PROJECT(X,P,IND)
       IF(IND.EQ.100)RETURN
       CALL DDER(X,V)
       CALL DER(X,VDER)
       CALL MATP(P,V,AUX1,NT,NT,NT,NT)
       CALL MATP(AUX1,P,V,NT,NT,NT,NT)      
       CALL JACOBI(V,NT,NT,VALP,VV,NROT)
       NN=0
       DO I=1,NT
        IF(ABS(VALP(I)).GT.PREC.AND.VALP(I).LT.0.0)NN=NN+1
       ENDDO
       IF(NN.GT.1)THEN
        IND=100
        RETURN
       ENDIF
      ENDDO
      IND=100
      RETURN
 20   IND=0
      CALL DDER(X,V)
      CALL DER(X,VDER)
      CALL MATP(P,V,AUX1,NT,NT,NT,NT)
      CALL MATP(AUX1,P,V,NT,NT,NT,NT)
      DO I=1,NT
      DO J=1,NT
       RM(I,J)=0.0
      ENDDO
      ENDDO      
      DO I=1,NA
      DO K=1,3
      J=3*(I-1)+K
      RM(J,J)=1.0/SQRT(RMM(I))
      ENDDO
      ENDDO
      CALL MATP(RM,V,AUX1,NT,NT,NT,NT)
      CALL MATP(AUX1,RM,V,NT,NT,NT,NT)
      CALL JACOBI(V,NT,NT,VALP,VV,NROT)
      RETURN
      END

      SUBROUTINE DER(X,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT)
      DO I=1,NT
      VAL(I)=DVL(X,PREC,I)
      ENDDO
      RETURN
      END

      SUBROUTINE DERR(X,VAL,CINT)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, ND=NA*(NA-1)/2,NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(ND),CINT(ND)
      CALL CARTINT(X,CINT)
      DO I=1,ND
      VAL(I)=DVLR(CINT,PREC,I)
      ENDDO
      RETURN
      END


      FUNCTION DVLR(CINT,STEP,I)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NA=3, ND=NA*(NA-1)/2,NT=3*NA )
      DIMENSION XA(NT),CINT(ND),XD(ND),R(ND)
      DO IA=1,ND
      XD(IA)=CINT(IA)
      ENDDO
      XD(I)=CINT(I)-STEP
      CALL INTCART(XD,XA)
      CALL TRANS(XA,R)  
      FMIN1=POTEN(R)
      XD(I)=CINT(I)-2*STEP
      CALL INTCART(XD,XA)
      CALL TRANS(XA,R)  
      FMIN2=POTEN(R)
      XD(I)=CINT(I)+STEP
      CALL INTCART(XD,XA)
      CALL TRANS(XA,R)
      FMAX1=POTEN(R)
      XD(I)=CINT(I)+2*STEP
      CALL INTCART(XD,XA)
      CALL TRANS(XA,R)
      FMAX2=POTEN(R)
      DVLR=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      SUBROUTINE DDER(X,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT,NT)
      DO I=1,NT
      DO J=I,NT
      VAL(I,J)=D2VL(X,PREC,I,J)
      ENDDO
      ENDDO
      DO I=1,NT
      DO J=I,NT
      VAL(J,I)=VAL(I,J)
      ENDDO
      ENDDO
      RETURN
      END

      FUNCTION DVL(X,STEP,I)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION X(NT),XD(NT)
      DO IN=1,NT
      XD(IN)=X(IN)
      ENDDO
      XD(I)=X(I)-STEP
      FMIN1=VL(XD)
      XD(I)=X(I)+STEP
      FMAX1=VL(XD)
      XD(I)=X(I)-2*STEP
      FMIN2=VL(XD)
      XD(I)=X(I)+2*STEP
      FMAX2=VL(XD)
      DVL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      FUNCTION D2VL(X,STEP,I,J)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION X(NT),XD(NT)
      DO IN=1,NT
      XD(IN)=X(IN)
      ENDDO
      XD(J)=X(J)-STEP
      FMIN1=DVL(XD,STEP,I)
      XD(J)=X(J)+STEP
      FMAX1=DVL(XD,STEP,I)
      XD(J)=X(J)-2*STEP
      FMIN2=DVL(XD,STEP,I)
      XD(J)=X(J)+2*STEP
      FMAX2=DVL(XD,STEP,I)
      D2VL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      FUNCTION VL(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION X(NT),R(NA*(NA-1)/2)
      CALL TRANS(X,R)
      VL=POTEN(R)
      RETURN
      END


      SUBROUTINE TRANS(X,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=3,NT=3*NA)
      DIMENSION X(NT),R(NA*(NA-1)/2)
        K=0
        DO 1 I=1,NT,3
          DO 1 J=I,NT,3
            IF(I.NE.J)THEN
              K=K+1
              R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+
     &          (X(I+2)-X(J+2))**2)
            ENDIF
 1        CONTINUE
      RETURN
      END

      SUBROUTINE CARTINT(X,CINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3,ND=NA*(NA-1)/2,NT=3*NA)
      DIMENSION GEO(3,NA),COORD(3,NA)
      DIMENSION X(NT),CINT(ND)
      COMMON /CONECT/NNA(NA),NB(NA),NC(NA)
      DEGREE=1.0
      DO I=1,NA
        COORD(1,I)=X(3*(I-1)+1)
        COORD(2,I)=X(3*(I-1)+2)
        COORD(3,I)=X(3*(I-1)+3)
      ENDDO
      CALL XYZINT(COORD,NNA,NB,NC,DEGREE,GEO)
      CINT(1)=GEO(1,2)
      CINT(2)=GEO(1,3)
      CINT(3)=GEO(2,3)
      DO I=3,NA
        CINT(3*(I-3)+1)=GEO(1,I)
        CINT(3*(I-3)+2)=GEO(2,I)
        CINT(3*(I-3)+3)=GEO(3,I)
      ENDDO
      END
      
      SUBROUTINE INTCART(CINT,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3,ND=NA*(NA-1)/2,NT=3*NA)
      DIMENSION GEO(3,NA),COORD(3,NA)
      DIMENSION X(NT),CINT(ND)
      COMMON /CONECT/NNA(NA),NB(NA),NC(NA)
      GEO(1,2)=CINT(1)
      GEO(1,3)=CINT(2)
      GEO(2,3)=CINT(3)
      DO I=3,NA
        GEO(1,I)=CINT(3*(I-3)+1)
        GEO(2,I)=CINT(3*(I-3)+2)
        GEO(3,I)=CINT(3*(I-3)+3)
      ENDDO
      CALL GMETRY(GEO,COORD,NNA,NB,NC)
      DO I=1,NA
        X(3*(I-1)+1)=COORD(1,I)
        X(3*(I-1)+2)=COORD(2,I)
        X(3*(I-1)+3)=COORD(3,I)
      ENDDO
      END
      

      SUBROUTINE GMETRY(GEO,COORD,NNA,NB,NC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3)
      DIMENSION GEO(3,NA),COORD(3,NA),NNA(NA),NB(NA),NC(NA)
C***********************************************************************
C
C    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
C *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
C    (C) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATESISDONE.
C
C  ON INPUT:
C         GEO    = ARRAY OF INTERNAL COORDINATES.
C         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
C         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
C
C  ON OUTPUT:
C         COORD  = ARRAY OF CARTESIAN COORDINATES
C
C***********************************************************************
      NATOMS=NA
      COORD(1,1)=0.0D00
      COORD(2,1)=0.0D00
      COORD(3,1)=0.0D00
      COORD(1,2)=GEO(1,2)
      COORD(2,2)=0.0D00
      COORD(3,2)=0.0D00
      IF(NATOMS.EQ.2) RETURN
      CCOS=COS(GEO(2,3))
      IF(NNA(3).EQ.1)THEN
         COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
      ELSE
         COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
      ENDIF
      COORD(2,3)=GEO(1,3)*SIN(GEO(2,3))
      COORD(3,3)=0.0D00
      DO 90 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NB(I)
         MC=NNA(I)
         XB=COORD(1,MB)-COORD(1,MC)
         YB=COORD(2,MB)-COORD(2,MC)
         ZB=COORD(3,MB)-COORD(3,MC)
         RBC=XB*XB+YB*YB+ZB*ZB
         IF(RBC.LT.1.D-16)THEN
C
C     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
C
            WRITE(6,'(A,I4,A,I4,A)')' ATOMS',MB,' AND',MC,' ARE COINCIDE
     1NT'
            WRITE(6,'(A)')' THIS IS A FATAL ERROR, RUN STOPPED IN GEOMETRY
     1'
            STOP
         ELSE
            RBC=1.0D00/SQRT(RBC)
         ENDIF
         MA=NC(I)
         XA=COORD(1,MA)-COORD(1,MC)
         YA=COORD(2,MA)-COORD(2,MC)
         ZA=COORD(3,MA)-COORD(3,MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 40
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=+1
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
   40    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
         ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
         YZA=SQRT(YPA**2+ZQA**2)
         IF(YZA.LT.1.D-4)GOTO 60
         IF(YZA.LT.2.D-2)THEN
            WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
            WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',/
     110X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
            WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREEES, ARE '',
     1''IN AN ALMOST STRAIGHT'')')
            WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE
     1'',/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
            WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I4)')I
            STOP
            WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')
     1')
            WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
            DO 50 J=1,I
   50       WRITE(6,'(I6,F16.5,2F13.5)')J,(COORD(K,J),K=1,3)
            WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
     1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
     2MC,MB,MA,YZA
            STOP
         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 70
   60    CONTINUE
C
C   ANGLE TOO SMALL TO BE IMPORTANT
C
         COSKH=1.D0
         SINKH=0.D0
   70    CONTINUE
C
C     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 80
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   80    COORD(1,I)=XQD+COORD(1,MC)
         COORD(2,I)=YQD+COORD(2,MC)
         COORD(3,I)=ZQD+COORD(3,MC)
   90 CONTINUE
      RETURN
      END

      SUBROUTINE XYZINT(XYZ,NNA,NB,NC,DEGREE,GEO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3)
      DIMENSION XYZ(3,NA),NNA(NA),NB(NA),NC(NA),GEO(3,NA)
***********************************************************************
*
* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
*        IF POSSIBLE.
*
*        IF(NA(2).EQ.1 THEN THE ORIGINAL CONNECTIVITY IS USED.
*
*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
*
*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN DEGREES
*
***********************************************************************
         NUMAT=NA
         DO 30 I=1,NUMAT
            NNA(I)=2
            NB(I)=3
            NC(I)=4
            IM1=I-1
            IF(IM1.EQ.0)GOTO 30
            SUM=1.D30
            DO 20 J=1,IM1
               R=(XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2
               IF(R.LT.SUM.AND.NNA(J).NE.J.AND.NB(J).NE.J) THEN
                  SUM=R
                  K=J
               ENDIF
   20       CONTINUE
C
C   ATOM I IS NEAREST TO ATOM K
C
            NNA(I)=K
            IF(I.GT.2)NB(I)=NNA(K)
            IF(I.GT.3)NC(I)=NB(K)
C
C   FIND ANY ATOM TO RELATE TO NA(I)
C
30    CONTINUE
      NNA(1)=0
      NB(1)=0
      NC(1)=0
      NB(2)=0
      NC(2)=0
      NC(3)=0
      CALL XYZGEO(XYZ,NUMAT,NNA,NB,NC,DEGREE,GEO)
      RETURN
      END
      
      SUBROUTINE XYZGEO(XYZ,NUMAT,NNA,NB,NC,DEGREE,GEO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3)
      DIMENSION XYZ(3,NA), NNA(NA), NB(NA), NC(NA), GEO(3,NA)
***********************************************************************
*
*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
*
*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
*              NUMAT= NUMBER OF ATOMS
*              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DISTANCE
*              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY ANGLE
*              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DIHEDRAL
*
*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
*                     AND RADIANS
*
***********************************************************************
      DO 30 I=2,NUMAT
         J=NNA(I)
         K=NB(I)
         L=NC(I)
         IF(I.LT.3) GOTO 30
         II=I
         CALL BANGLE(XYZ,II,J,K,GEO(2,I))
         GEO(2,I)=GEO(2,I)*DEGREE
         IF(I.LT.4) GOTO 30
C
C   MAKE SURE DIHEDRAL IS MEANINGLFUL
C
         CALL BANGLE(XYZ,J,K,L,ANGL)
         TOL=0.2617994D0
         IF(ANGL.GT.3.1415926D0-TOL.OR.ANGL.LT.TOL)THEN
C
C  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
C  DEFINING THE DIHEDRAL.
   10       SUM=100.D0
            DO 20 I1=1,II-1
               R=(XYZ(1,I1)-XYZ(1,K))**2+
     1          (XYZ(2,I1)-XYZ(2,K))**2+
     2          (XYZ(3,I1)-XYZ(3,K))**2
               IF(R.LT.SUM.AND.I1.NE.J.AND.I1.NE.K) THEN
                  CALL BANGLE(XYZ,J,K,I1,ANGL)
                  IF(ANGL.LT.3.1415926D0-TOL.AND.ANGL.GT.TOL)THEN
                     SUM=R
                     L=I1
                     NC(II)=L
                  ENDIF
               ENDIF
   20       CONTINUE
            IF(SUM.GT.99.D0.AND.TOL.GT.0.1D0)THEN
C
C ANYTHING WITHIN 5 DEGREES?
C
               TOL=0.087266D0
               GOTO 10
            ENDIF
         ENDIF
         CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
         GEO(3,I)=GEO(3,I)*DEGREE
   30 GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1                   (XYZ(2,I)-XYZ(2,J))**2+
     2                   (XYZ(3,I)-XYZ(3,J))**2)
      GEO(1,1)=0.D0
      GEO(2,1)=0.D0
      GEO(3,1)=0.D0
      GEO(2,2)=0.D0
      GEO(3,2)=0.D0
      GEO(3,3)=0.D0
      RETURN
      END
      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3)
      DIMENSION XYZ(3,NA)
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*
*********************************************************************
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )
      RETURN
      END

      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3)
      DIMENSION XYZ(3,NA)
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-6) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
      IF (ANGLE .LT. 0.) ANGLE=4.0D0* ASIN(1.0D00)+ANGLE
      IF (ANGLE .GE. 6.2831853D0 ) ANGLE=0.D0
      RETURN
      END
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
      IMPLICIT REAL*8 (A-H,O-Z)
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      ZERO=1.0D-6
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=4.0D0* ASIN(1.0D00)-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END
      SUBROUTINE MATP(A,B,C,N,L,M,ND)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(ND,L),B(ND,M),C(ND,M)
      DO I=1,N
         DO J=1,M
            C(I,J)=0.0
            DO K=1,L
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE PROJECT(X,P,IND)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION X(NT),B(NT,NT),S(NT,NT),SINV(NT,NT),
     &  AUX1(NT,NT),BTR(NT,NT),P(NT,NT),U(NT,NT)
      DO I=1,NT
      DO J=1,6
      B(I,J)=0.0
      ENDDO
      ENDDO
      DO I=1,NT
      DO J=1,NT
      U(I,J)=0.0
      ENDDO
      U(I,I)=1.0
      ENDDO
      DO I=1,NT,3
      B(I,  1)=1.0
      B(I+1,2)=1.0
      B(I+2,3)=1.0
      B(I+1,4)= X(I+2)
      B(I+2,4)=-X(I+1)
      B(I,  5)=-X(I+2)
      B(I+2,5)= X(I)      
      B(I,  6)= X(I+1)
      B(I+1,6)=-X(I)
      ENDDO
      RNT=REAL(NA)
      DO I=1,6
      DO J=1,6
      S(I,J)=0.0
      ENDDO
      ENDDO
      S(1,1)=RNT
      S(2,2)=RNT
      S(3,3)=RNT
      DO I=1,NT,3
      S(1,5)=S(1,5)-X(I+2)
      S(1,6)=S(1,6)+X(I+1)
      S(2,4)=S(2,4)+X(I+2)
      S(2,6)=S(2,6)-X(I)
      S(3,4)=S(3,4)-X(I+1)
      S(3,5)=S(3,5)+X(I)
      S(4,5)=S(4,5)-X(I)*X(I+1)
      S(4,6)=S(4,6)-X(I)*X(I+2)
      S(5,6)=S(5,6)-X(I+1)*X(I+2)
      S(4,4)=S(4,4)+X(I+1)**2+X(I+2)**2
      S(5,5)=S(5,5)+X(I)**2  +X(I+2)**2
      S(6,6)=S(6,6)+X(I)**2  +X(I+1)**2
      ENDDO
      S(6,1)=S(1,6)
      S(6,2)=S(2,6)
      S(6,3)=S(3,6)
      S(6,4)=S(4,6)
      S(6,5)=S(5,6)
      S(5,1)=S(1,5)
      S(5,2)=S(2,5)
      S(5,3)=S(3,5)
      S(5,4)=S(4,5)
      S(4,1)=S(1,4)
      S(4,2)=S(2,4)
      S(4,3)=S(3,4)
      S(3,1)=S(1,3)
      S(3,2)=S(2,3)
      S(2,1)=S(1,2)
      CALL MINVER(S,SINV,IND)
      IF(IND.EQ.100)RETURN
      CALL MATP(B,SINV,AUX1,NT,6,6,NT)
      DO I=1,6
        DO J=1,NT
          BTR(I,J)=B(J,I)
        ENDDO
      ENDDO
      CALL MATP(AUX1,BTR,P,NT,6,NT,NT)
      DO I=1,NT
      DO J=1,NT
      P(I,J)=U(I,J)-P(I,J)
      ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE JACOBI(AA,N,NP,D,V,NROT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200)
      DIMENSION AA(NP,NP),A(NP,NP),D(NP),V(NP,NP),
     &  B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.0
          A(IP,IQ)=AA(IP,IQ)
11      CONTINUE
        V(IP,IP)=1.0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END


      SUBROUTINE MINVER(A,AINV,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION RMAT(NT,NT+1),A(NT,NT),AINV(NT,NT),X(NT)
      ND=6
      DO I=1,ND
        DO J=1,ND
          RMAT(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,ND
      DO J=1,ND
        RMAT(J,ND+1)=0.0
        ENDDO
        RMAT(I,ND+1)=1.0
        CALL SVD(RMAT,ND,X,IND)
        IF(IND.EQ.100)RETURN
        DO J=1,ND
        AINV(J,I)=X(J)
        ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE SVD(RMAT,N,A,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      PARAMETER(NMAX=1000,MMAX=100,TOL=1.0D-9)
      DIMENSION RMAT(NT,NT+1),A(NT),V(NT,NT),
     *    U(NT,NT),W(NT),B(NMAX)
      DO 12 I=1,N
        DO 11 J=1,N
          U(I,J)=RMAT(I,J)
11      CONTINUE
        B(I)=RMAT(I,N+1)
12    CONTINUE
      CALL SVDCMP(U,N,N,W,V,IND)
      IF(IND.EQ.100)RETURN
      WMAX=0.0
      DO 13 J=1,N
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,N
        IF(W(J).LT.THRESH)W(J)=0.0
14    CONTINUE
      CALL SVBKSB(U,W,V,N,N,B,A)
      RETURN
      END

      SUBROUTINE SVDCMP(A,M,N,W,V,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=3, NT=3*NA)
      PARAMETER (NMAX=1000)
      DIMENSION A(NT,NT),W(NT),V(NT,NT),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30)THEN
          IND=100
          RETURN
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END


      SUBROUTINE SVBKSB(U,W,V,M,N,B,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=1000)
      PARAMETER (NA=3, NT=3*NA)
      DIMENSION U(NT,NT),W(NT),V(NT,NT),B(NMAX),X(NT),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END

        FUNCTION RANF ( DUMMY )
        IMPLICIT REAL*8 (A-H,O-Z)

C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        SAVE        SEED
        DATA        SEED / 0 /

C    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = REAL ( SEED ) / M

        RETURN
        END


      FUNCTION POTEN(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3)  
C      R(1)=RSH
C      R(2)=RSO
C      R(3)=ROH
       POTEN=poten_cm(R)   
c       POTEN=POTSO2(R)
c      POTEN=VHO2(R(1),R(2),R(3))
c      POTEN=DMBEPOT(R)
      RETURN
      END



