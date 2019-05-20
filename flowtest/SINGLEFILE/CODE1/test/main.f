      MODULE PHYPROJFISH
      USE USERFUNC 
      IMPLICIT NONE 
      REAL(KIND=8)::X_0,X_1,Y_0,Y_1
      INTEGER(KIND=4)::M,N
*>XC,YC:ZONE GRID NODES;X,Y:COMPUTING DATA NODES
*>HKSI,HETA,VOL,AE,AK:THE DOMAIN VOLUME'S LENGTH,HEIGHT,VOLUME,
*>LEFT SIDE SURFACE AREA AND BOTTOM SURFACE AREA SEPARATELY.
      REAL(KIND=8),ALLOCATABLE::XC(:,:),YC(:,:),X(:,:),Y(:,:)
      REAL(KIND=8),ALLOCATABLE::HKSI(:,:),HETA(:,:),VOL(:,:)
      REAL(KIND=8),ALLOCATABLE::AK(:,:),AE(:,:)
*>CHORIN PROJECTION METHOD
      REAL(KIND=8),ALLOCATABLE::RHO(:,:),U(:,:),V(:,:)
      REAL(KIND=8),ALLOCATABLE::P(:,:),TEM(:,:),QST(:,:)
*>PROJSOLA CONTROL PARAMETERS****************************************
      INTEGER(KIND=4)::MODEKSI=0,SKIP=1
      REAL(KIND=8)::ALPHA=0.0,OMIC
      REAL(KIND=8)::DT=1.0E-4,DX,DY
      REAL(KIND=8)::RE=100,GR=0.0
      REAL(KIND=8)::RESIERR=1.0E-4,MAXDIV
      
      CONTAINS
      SUBROUTINE CYLINDER2DMESH
*>USED IN FLOW SIMULATING PROCESS***************************
*>ONE ADDITIONAL VIRTUAL GRID*******************************
*>X0<X1,Y0<Y1;**********************************************
*>X--->R,Y--->Z;(X,Y) RELATE TO P POINT*********************
*>FUTURE,ZONE GRID NODES CAN BE MODIFIED********************
*>AXIAL-SYSMETRY CYLINDERIAL COORDINATES********************
      IMPLICIT NONE
      REAL(KIND=8),PARAMETER::PI=3.1415926
      INTEGER(KIND=4)::I,J
      ALLOCATE(XC(0:M+2,0:N+2),YC(0:M+2,0:N+2))
      ALLOCATE(X(1:M+2,1:N+2),Y(1:M+2,1:N+2))
      ALLOCATE(HKSI(1:M+2,1:N+2),HETA(1:M+2,1:N+2))
      ALLOCATE(VOL(1:M+2,1:N+2),AK(1:M+2,1:N+2),AE(1:M+2,1:N+2))

      DO I=0,M+2
      XC(I,:)=X_0+(I-1)*(X_1-X_0)/M
      ENDDO
      DO J=0,N+2
      YC(:,J)=Y_0+(J-1)*(Y_1-Y_0)/N
      ENDDO

      DO J=1,N+2
      DO I=1,M+2
      X(I,J)=0.25*(XC(I-1,J-1)+XC(I,J-1)+XC(I-1,J)+XC(I,J))
      Y(I,J)=0.25*(YC(I-1,J-1)+YC(I,J-1)+YC(I-1,J)+YC(I,J))
      ENDDO
      ENDDO

      HKSI=0;HETA=0;VOL=0;AK=0;AE=0
      DO J=1,N+2
      DO I=1,M+2
      HKSI(I,J)=0.5*((XC(I,J)+XC(I,J-1))-(XC(I-1,J)+XC(I-1,J-1)))
      HETA(I,J)=0.5*((YC(I,J)+YC(I-1,J))-(YC(I,J-1)+YC(I-1,J-1)))
      VOL(I,J)=HKSI(I,J)*HETA(I,J)*2*PI*X(I,J)
      AK(I,J)=XC(I,J-1)-XC(I-1,J-1)
      AE(I,J)=YC(I-1,J)-YC(I-1,J-1)
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CYLINDER2DMESH 
      
      SUBROUTINE CFL
      IMPLICIT NONE
      REAL(KIND=8)::UMAX,VMAX
      UMAX=MAXVAL(U,0,M+2,1,N+2)
      VMAX=MAXVAL(V,1,M+2,0,N+2)
      DT=0.25*MIN(DX/UMAX,DY/VMAX)
      ALPHA=1.2*MAX(ABS(UMAX*DT/DX),ABS(VMAX*DT/DY))
      WRITE(*,200)DT,ALPHA
200   FORMAT(1X,'DT=',F8.4,1X,'ALPHA=',F8.4,1X)
      RETURN
      END SUBROUTINE CFL

      SUBROUTINE INITIAL
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      ALLOCATE(RHO(1:M+2,1:N+2),U(0:M+2,1:N+2),V(1:M+2,0:N+2))
      ALLOCATE(P(1:M+2,1:N+2),QST(1:M+1,1:N+1))
*>PROJSOLA CONTROL PARAMETERS****************************************
      RHO=1.0;U=0.0;V=0.0;P=1.0;QST=0
      CALL BOUNDARYU
      CALL BOUNDARYV
      CALL BOUNDARYP
      RETURN
      END SUBROUTINE INITIAL

      SUBROUTINE BOUNDARYU
      IMPLICIT NONE
      U(1,2:N+1)=0
      U(M+1,2:N+1)=0
      U(2:M,1)=0.0*2.0-U(2:M,2)
      U(2:M,N+2)=1.0*2.0-U(2:M,N+1)
      RETURN
      END SUBROUTINE BOUNDARYU

      SUBROUTINE BOUNDARYV
      IMPLICIT NONE
      V(2:M+1,N+1)=0
      V(2:M+1,1)=0
      V(M+2,2:N)=0.0*2.0-V(M+1,2:N)
      V(1,2:N)=0.0*2.0-V(2,2:N)
      RETURN
      END SUBROUTINE BOUNDARYV

      SUBROUTINE BOUNDARYP
      IMPLICIT NONE
      P(2:M+1,1)=P(2:M+1,2)
      P(2:M+1,N+2)=P(2:M+1,N+1)
      P(1,2:N+1)=P(2,2:N+1)
      P(M+2,2:N+1)=P(M+1,2:N+1)
      RETURN
      END SUBROUTINE BOUNDARYP

      SUBROUTINE FIRSTSOLVEUV
*>*************SOLVE FIRST STEP VELOCITY*******************************
*>A:COMPUTE THE FIRST VELOCITY USING DISPLAY CONVECTIVE TERM
*>B:COMPUTE THE INTERMIDDLE VELOCITY USING IMPLICIT VISCOSITY TERM(FORCE)
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD
      REAL(KIND=8)::FUX,FUY,FUC,FVX,FVY,FVC,BUOY
      REAL(KIND=8)::UT(0:M+2,1:N+2),VT(1:M+2,0:N+2)
      UT=U;VT=V
*>A:DISPLAY PART,THE FIRST VELOCITY
      DO J=2,N+1
      DO I=2,M
*>U INTERIOR POINTS WITHOUT BOUNDARY POINTS
      DX=X(I+1,J)-X(I,J)
      FUX=0.25/DX*((U(I,J)+U(I+1,J))**2+
     1             ALPHA*ABS(U(I,J)+U(I+1,J))*(U(I,J)-U(I+1,J))-
     1             (U(I-1,J)+U(I,J))**2-
     1             ALPHA*ABS(U(I-1,J)+U(I,J))*(U(I-1,J)-U(I,J)))
      DY=YC(I,J)-YC(I,J-1)
      FUY=0.25/DY*((V(I,J)+V(I+1,J))*(U(I,J)+U(I,J+1))+
     1             ALPHA*ABS(V(I,J)+V(I+1,J))*(U(I,J)-U(I,J+1))-
     1             (V(I,J-1)+V(I+1,J-1))*(U(I,J-1)+U(I,J))-
     1             ALPHA*ABS(V(I,J-1)+V(I+1,J-1))*(U(I,J-1)-U(I,J))) 
      RAD=0.5*(X(I,J)+X(I+1,J))
      FUC=MODEKSI/8.0/RAD*((U(I,J)+U(I+1,J))**2+
     1             (U(I-1,J)+U(I,J))**2+
     1             ALPHA*ABS(U(I,J)+U(I+1,J))*(U(I,J)-U(I+1,J))+
     1             ALPHA*ABS(U(I-1,J)+U(I,J))*(U(I-1,J)-U(I,J)))
      UT(I,J)=(-FUX-FUY-FUC)*DT+U(I,J)
      IF(UT(I,J).EQ.2.0)WRITE(*,*)'UT=2'
      ENDDO
      ENDDO

      DO J=2,N
      DO I=2,M+1
*>V INTERIOR POINTS
      DX=XC(I,J)-XC(I-1,J)
      FVX=0.25/DX*((U(I,J)+U(I,J+1))*(V(I,J)+V(I+1,J))+
     1             ALPHA*ABS(U(I,J)+U(I,J+1))*(V(I,J)-V(I+1,J))-
     1             (U(I-1,J)+U(I-1,J+1))*(V(I-1,J)+V(I,J))-
     1             ALPHA*ABS(U(I-1,J)+U(I-1,J+1))*(V(I-1,J)-V(I,J)))
      DY=Y(I,J+1)-Y(I,J)
      FVY=0.25/DY*((V(I,J)+V(I,J+1))**2+
     1             ALPHA*ABS(V(I,J)+V(I,J+1))*(V(I,J)-V(I,J+1))-
     1             (V(I,J-1)+V(I,J))**2-
     1             ALPHA*ABS(V(I,J-1)+V(I,J))*(V(I,J-1)-V(I,J)))
      RAD=0.5*(XC(I,J)+XC(I-1,J)) 
      FVC=MODEKSI/8.0/RAD*((U(I,J)+U(I,J+1))*(V(I,J)+V(I+1,J))+
     1             (U(I-1,J)+U(I-1,J+1))*(V(I-1,J)+V(I,J))+
     1             ALPHA*ABS(U(I,J)+U(I,J+1))*(V(I,J)-V(I+1,J))+
     1             ALPHA*ABS(U(I-1,J)+U(I-1,J+1))*(V(I-1,J)-V(I,J))) 
*c      BUOY=GR/RE/RE*TEM(I,J)
      BUOY=0
      VT(I,J)=(BUOY-FVX-FVY-FVC)*DT+V(I,J)
      ENDDO
      ENDDO
*>FIRST UPDATE VELOCITY USING BOUNDARY CONDITION
      U=UT;V=VT
      CALL BOUNDARYU
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE FIRSTSOLVEUV

      SUBROUTINE SECONDSOLVEUV
*>IMPLICIT VISCOSITY TERM FOR INTERMIDDLE VELOCITY
*>SOLVING UMAT MATRIX(WHICH CONCLUDES ALL MESH POINTS)
*>SO THE DIMENSION OF UMAT EQUALS THAT OF U
*>ONLY THE DOMAIN PIVOT ELEMENT OF UMAT ARE SETTED TO 1,
*>THE COORESPOING BU EQUALS THE BOUNDARY VALUE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
c      REAL(KIND=8)::BDA(N),BDB(N),BDC(M),BDD(M),F(M,N)
      REAL(KIND=8)::BDA(103),BDB(103),BDC(98),BDD(98),F(98,103)
      INTEGER(KIND=4)::MBDCND,NBDCND,IDIMF,IERROR
      REAL(KIND=8)::ELMBDA,PERTRB
      REAL(KIND=8)::W(5*98*103)
*>U,MBDCND:1,NBDCND=1
*>USING HSTCRT SUBROUTINE
      MBDCND=1;NBDCND=1;IDIMF=98
c      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=1.0
      DO J=1,103
      BDA(J)=0
      BDB(J)=0
      ENDDO
      DO I=1,98
      BDC(I)=0
      BDD(I)=1
      ENDDO
c      ELMBDA=-RE/DT 
      ELMBDA=-1.0E5
      DO J=1,103
      DO I=1,98
c      DO J=1,N
c      DO I=1,M
c      F(I,J)=-RE/DX*0.5*(U(I,J+1)+U(I+1,J+1))
      F(I,J)=0
      ENDDO
      ENDDO
c      CALL HSTCRT(X_0,X_1,M,MBDCND,BDA,BDB,
c     1            Y_0,Y_1,N,NBDCND,BDC,BDD,
c     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      CALL HSTCRT(1.0,3.0,98,MBDCND,BDA,BDB,
     1            -1.0,1.0,103,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,*)'SOLVE U IERROR',IERROR
      OPEN(10,FILE='F1.DAT')
      DO J=1,103
      DO I=1,98
      WRITE(10,'(ES20.5)')F(I,J)
      ENDDO
      ENDDO
      CLOSE(10)
      STOP
*>TRANSFORM F(I,J) TO U
      DO J=1,N
      DO I=1,M
      U(I+1,J+1)=0.5*(F(I,J)+F(I+1,J))
      ENDDO
      ENDDO
      CALL BOUNDARYU
*>V
      MBDCND=1;NBDCND=1;IDIMF=M
      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0.0
      ELMBDA=-RE/DT 
      DO J=1,N
      DO I=1,M
      F(I,J)=-RE/DX*0.5*(V(I+1,J+1)+V(I+1,J))
      ENDDO
      ENDDO
      CALL HSTCRT(X_0,X_1,M,MBDCND,BDA,BDB,
     1            Y_0,Y_1,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,*)'SOLVE V IERROR',IERROR
*>TRANSFORM F(I,J) TO V
      DO J=1,N
      DO I=1,M
      V(I+1,J+1)=0.5*(F(I,J)+F(I,J+1))
      ENDDO
      ENDDO
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE SECONDSOLVEUV

      SUBROUTINE SOLVEPRESSURE
*>*************USE FIRST STEP VELOCITY*******************************
*>*************SOLVE THE NEXT INTERVAL PRESSURE************************
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::BDA(N),BDB(N),BDC(M),BDD(M),F(M,N)
      INTEGER(KIND=4)::MBDCND,NBDCND,IDIMF,IERROR
      REAL(KIND=8)::ELMBDA,PERTRB
      REAL(KIND=8)::W(5*M*N)
*>P,MBDCND:3,NBDCND=3
*>USING HSTCRT SUBROUTINE
      MBDCND=3;NBDCND=3;IDIMF=M
      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0
      ELMBDA=0.0
      DO J=1,N
      DO I=1,M
      F(I,J)=-1/DT*((U(I+1,J+1)-U(I,J+1))/DX+
     1              (V(I+1,J+1)-V(I+1,J))/DY)
      ENDDO
      ENDDO
      CALL HSTCRT(X_0,X_1,M,MBDCND,BDA,BDB,
     1            Y_0,Y_1,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,*)'SOLVE P IERROR',IERROR
*>TRANSFORM F(I,J) TO P
      DO J=1,N
      DO I=1,M
      P(I+1,J+1)=F(I,J)
      ENDDO
      ENDDO
*>BOUNDARY P
      CALL BOUNDARYP
      RETURN
      END SUBROUTINE SOLVEPRESSURE

      SUBROUTINE MODIFYUV
      IMPLICIT NONE
      LOGICAL,EXTERNAL::ISCONS
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::DIV,DP,RAD
      DO J=2,N+1
      DO I=2,M+1
      DX=HKSI(I,J);DY=HETA(I,J);RAD=X(I,J)
      DIV=(U(I+1,J)-U(I,J))/DX+(V(I,J+1)-V(I,J))/DY+
     1       (U(I+1,J)+U(I,J))/2/RAD*MODEKSI 
      DP=-DIV/2/DT/(1/DX/DX+1/DY/DY)
      P(I,J)=P(I,J)+DP
      ENDDO
      ENDDO
      CALL BOUNDARYU
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE MODIFYUV

      FUNCTION ISCONS(RESIERR)
      IMPLICIT NONE
      LOGICAL::ISCONS
      REAL(KIND=8)::RESIERR
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD
      DO J=2,N+1
      DO I=2,M+1
      DX=HKSI(I,J);DY=HETA(I,J);RAD=X(I,J)
      MAXDIV=(U(I+1,J)-U(I,J))/DX+(V(I,J+1)-V(I,J))/DY+
     1       (U(I+1,J)+U(I,J))/2/RAD*MODEKSI
      IF(MAXDIV.GT.RESIERR)THEN
      ISCONS=.FALSE.
      EXIT
      ENDIF
      ENDDO
      ISCONS=.TRUE.
      EXIT
      ENDDO
      RETURN
      END FUNCTION ISCONS

      SUBROUTINE PROJSOLA
      IMPLICIT NONE
*>*************SOLVE SECOND STEP VELOCITY******************************
*>*************MODIFY THE PRESSURE AND VELOCITY************************
*>A:COMPUTE THE INTERMEDIATE VELOCITY DIVERGENCE
*>B:SOLVE POISSON EQUATION
*>C:COMPUTE THE PRESSURE GRADIENT
*>D:UPDATE VELOCITY FIELD
*>*************PROJSOLA METHOD*****************************************
      INTEGER(KIND=4)::NCOUNT
      CALL FIRSTSOLVEUV
      CALL SECONDSOLVEUV
      CALL SOLVEPRESSURE
      WRITE(*,*)'MAXPRESSURE'
      WRITE(*,*)MAXVAL(P,1,M+2,1,N+2)
      STOP
      NCOUNT=0
      WRITE(*,80)'ITERATION NCOUNT=',NCOUNT,
     1           'MAX DIVERGENCE=',MAXDIV
      DO WHILE(.NOT.(ISCONS(RESIERR).OR.NCOUNT.GT.200))
      CALL MODIFYUV
      NCOUNT=NCOUNT+1
      WRITE(*,80)'ITERATION NCOUNT=',NCOUNT,
     1           'MAX DIVERGENCE=',MAXDIV
80    FORMAT(1X,A20,1X,I6,1X,A20,ES20.5)
      ENDDO
      RETURN
      END SUBROUTINE PROJSOLA

      SUBROUTINE OUTPUT(N)
*>VELOCITY AND PRESSURE ...USEING UT,VT,PT...*******************
*>THESE VARIABLES ARE INTERPLOTED TO THE GRID POINT(XC,YC)******
*>(QST(1:M+1,1:N+1));(RHO(1:M+2,1:N+2),U(0:M+2,1:N+2),V(1:M+2,0:N+2))
*>(P(1:M+2,1:N+2));UT(0:M+2,1:N+2),VT(1:M+2,0:N+2)
      IMPLICIT NONE
      INTEGER(KIND=4),INTENT(IN)::N
      INTEGER(KIND=4)::I,J
      CHARACTER(LEN=20)::FILE1
      REAL(KIND=8)::UT(0:M+2,1:N+2),VT(1:M+2,0:N+2)
      WRITE(FILE1,'(I6)')N
      FILE1='0'//TRIM(ADJUSTL(FILE1))//'.DAT'
*>STREAMLINE FUNCTION
      CALL STREAMLINE
*>VELOCITY
      UT=0;VT=0
      DO J=2,N
      DO I=2,M
      UT(I,J)=0.5*(U(I,J)+U(I,J+1))
      VT(I,J)=0.5*(V(I,J)+V(I+1,J))
      ENDDO
      ENDDO 
      OPEN(10,FILE=FILE1)
      WRITE(10,'(A)')'VARIABLES="X","Y","STRL","U","V"'
      WRITE(10,100)M+1,N+1
      DO J=1,N+1
      DO I=1,M+1
      WRITE(10,101)XC(I,J),YC(I,J),QST(I,J),U(I,J),V(I,J)
      ENDDO
      ENDDO
      CLOSE(10)
100   FORMAT(1X,'ZONE T="FLOW",I=',I4,' , J=',I4,' , F=POINT')
101   FORMAT(1X,2F10.5,3ES20.5)
      RETURN
      END SUBROUTINE OUTPUT

      SUBROUTINE STREAMLINE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::BDA(N),BDB(N),BDC(M),BDD(M),F(M+1,N+1)
      INTEGER(KIND=4)::MBDCND,NBDCND,IDIMF,IERROR
      REAL(KIND=8)::ELMBDA,PERTRB
      REAL(KIND=8)::W(5*M*N)
*>STREAM LINE,MBDCND:3,NBDCND=3
*>USING HWSCRT SUBROUTINE
      MBDCND=1;NBDCND=1;IDIMF=M+1
      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0
      ELMBDA=0.0
      F=0
      DO J=2,N
      DO I=2,M
      F(I,J)=(U(I,J+1)-U(I,J))/DY-
     1       (V(I+1,J)-V(I,J))/DX 
      ENDDO
      ENDDO
      CALL HWSCRT(X_0,X_1,M,MBDCND,BDA,BDB,
     1            Y_0,Y_1,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
*>TRANSFORM F(I,J) TO P
      DO J=1,N
      DO I=1,M
      P(I+1,J+1)=F(I,J)
      ENDDO
      ENDDO
*>BOUNDARY P
      CALL BOUNDARYP
      RETURN
      END SUBROUTINE STREAMLINE
      END MODULE PHYPROJFISH
      
      PROGRAM MAIN
      USE PHYPROJFISH
      USE USERFUNC
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J,NOUT=1
      REAL(KIND=8)::TEND=1.0,T=0.0
      INTEGER(KIND=4)::CONTROL=0
      X_0=0.0;X_1=1.0;Y_0=0.0;Y_1=1.0
      M=100;N=100
      DX=(X_1-X_0)/M;DY=(Y_1-Y_0)/N
      CALL CYLINDER2DMESH
      WRITE(*,20)M,N,DX,DY
20    FORMAT('MESH GRID M*N',2I6,1x,'DX DY',2F8.4)
      CALL INITIAL
      WRITE(*,60)'TIME=',T,'INITIAL DT=',DT
      DO WHILE(T.LE.TEND)
      CALL CFL
      T=T+DT
      WRITE(*,60)'TIME=',T,'CFL DT=',DT
      CALL PROJSOLA
      IF(T.GE.(0.01*NOUT))THEN
      CALL OUTPUT(NOUT)
      NOUT=NOUT+1
      ENDIF
      ENDDO
60    FORMAT(1X,A10,F10.5,1X,A10,F10.5)
61    FORMAT(1X,A10,F10.5,'WRITE FILE',1X,A10,I10)

      STOP
      END PROGRAM MAIN

