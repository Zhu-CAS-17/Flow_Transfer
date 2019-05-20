*>TEST SEMI IN AN OUT VELOCITY BOUNDARY CONDITION
*>CHECK THE COUPLE OF VELOCITY AND PRESSURE
      MODULE PHYPROJFISH
      USE USERFUNC 
      IMPLICIT NONE 
      REAL(KIND=8)::A,B,C,D
      INTEGER(KIND=4)::M,N
*>XC,YC:ZONE GRID NODES;X,Y:COMPUTING DATA NODES
*>HKSI,HETA,VOL,AE,AK:THE DOMAIN VOLUME'S LENGTH,HEIGHT,VOLUME,
*>LEFT SIDE SURFACE AREA AND BOTTOM SURFACE AREA SEPARATELY.
      REAL(KIND=8),ALLOCATABLE::XC(:,:),YC(:,:),X(:,:),Y(:,:)
      REAL(KIND=8),ALLOCATABLE::HKSI(:,:),HETA(:,:),VOL(:,:)
      REAL(KIND=8),ALLOCATABLE::AK(:,:),AE(:,:)
*>CHORIN PROJECTION METHOD
      REAL(KIND=8),ALLOCATABLE::RHO(:,:),U(:,:),V(:,:),P(:,:)
      REAL(KIND=8),ALLOCATABLE::COND(:,:),TEM(:,:),QST(:,:)
      REAL(KIND=8),ALLOCATABLE::TEMWRITE(:,:)
*>PROJSOLA CONTROL PARAMETERS****************************************
      INTEGER(KIND=4)::MODEKSI=0,SKIP=1
      REAL(KIND=8)::ALPHA=1.0,OMIC
      REAL(KIND=8)::DT=1.0E-3,DX,DY
      REAL(KIND=8)::RE=1,GR=1000.0
      REAL(KIND=8)::RESIERR=1.0E-4,MAXDIV
      INTEGER(KIND=4)::MAXINT=100
      LOGICAL(KIND=4)::ISCONS
*>READIN DATAS AND PARAMETERS
      REAL(KIND=8),ALLOCATABLE::XIN(:),YIN(:),TIN(:)
      
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
      XC(I,:)=A+(I-1)*(B-A)/M
      ENDDO
      DO J=0,N+2
      YC(:,J)=C+(J-1)*(D-C)/N
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
      CALL MAXINMAT(U,0,M+2,1,N+2,UMAX)
      CALL MAXINMAT(V,1,M+2,0,N+2,VMAX)
      DT=0.25*MIN(DX/UMAX,DY/VMAX)
      ALPHA=1.2*MAX(ABS(UMAX*DT/DX),ABS(VMAX*DT/DY))
      WRITE(*,200)DT,ALPHA
200   FORMAT(1X,'DT=',F8.4,1X,'ALPHA=',F8.4,1X)
      RETURN
      END SUBROUTINE CFL

      SUBROUTINE READIN
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::NUMX,NUMY,NUMXY
      REAL(KIND=8)::YINMAX,TINMAX,TINMIN
      OPEN(10,FILE='chamber_temp.dat')
      READ(10,*)
      READ(10,*)
      READ(10,*)
*>      READ(10,30)NUMX,NUMY
      READ(10,*)NUMX,NUMY
      NUMXY=NUMX*NUMY
30    FORMAT(' zone t= "chamber",i=',I3,' , j=',I3,' , f= point') 
      ALLOCATE(XIN(NUMXY),YIN(NUMXY),TIN(NUMXY))
      DO I=1,NUMXY
      READ(10,'(3F10.4)')XIN(I),YIN(I),TIN(I)
      ENDDO
      CLOSE(10)
      CALL MAXINMAT(YIN,1,NUMXY,1,1,YINMAX)
      CALL MAXINMAT(TIN,1,NUMXY,1,1,TINMAX)
      CALL MININMAT(TIN,1,NUMXY,1,1,TINMIN)
      XIN=XIN/0.0254;YIN=(YIN-YINMAX)/0.0254+1
      CALL MATRIXTOMATRIX(TIN,XIN,YIN,NUMX,NUMY,
     1            TEM,X(2:M+1,1),Y(2,1:N+1),M,N) 
      TEM=(TEM-TINMIN)/(TINMAX-TINMIN)
      CALL MATRIXTOMATRIX(TIN,XIN,YIN,NUMX,NUMY,
     1            TEMWRITE,XC(1:M+1,1),YC(1,1:N+1),M+1,N+1) 
      TEMWRITE=(TEMWRITE-TINMIN)/(TINMAX-TINMIN) 
      RETURN
      END SUBROUTINE READIN

      SUBROUTINE INITIAL
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      ALLOCATE(RHO(1:M+2,1:N+2),U(0:M+2,1:N+2),V(1:M+2,0:N+2))
      ALLOCATE(P(1:M+2,1:N+2),QST(1:M+1,1:N+1))
      ALLOCATE(COND(1:M+2,1:N+2),TEM(2:M+1,2:N+1),TEMWRITE(1:M+1,1:N+1))
*>PROJSOLA CONTROL PARAMETERS****************************************
      RHO=1.0;U=0.0;V=0.0;P=1.0;QST=0;COND=0;TEM=0
      CALL BOUNDARYU
      CALL BOUNDARYV
      CALL BOUNDARYP
*>SPLINE INTERPLATION,TEM(1:M+1,1:N+1_
      CALL READIN
      RETURN
      END SUBROUTINE INITIAL

      SUBROUTINE BOUNDARYU
      IMPLICIT NONE
      INTEGER(KIND=4)::J
*>SYMMETRY AXIS U(1,2:N+1)=0;WALL U(M+1,2:N+1)=0
*>      U(1,2:N+1)=0
      U(M+1,2:N+1)=0
      U(2:M,1)=0.0*2.0-U(2:M,2)
      U(2:M,N+2)=0.0*2.0-U(2:M,N+1)
*>UNUIFORM BOUNDARY CONDITION
      DO J=2,N+1
      IF(0.5*(YC(1,J)+YC(1,J-1)).LE.0.3)THEN
      U(1,J)=1.0
      ELSE
      U(1,J)=0.0
      ENDIF 
      ENDDO
      RETURN
      END SUBROUTINE BOUNDARYU

      SUBROUTINE BOUNDARYV
      IMPLICIT NONE
      INTEGER(KIND=4)::I
*>WALL V(2:M+1,N+1)=0
*>      V(2:M+1,N+1)=0
      DO I=2,M+1
      IF(0.5*(XC(I,N+1)+XC(I-1,N+1)).LE.0.3)THEN
      V(I,N+1)=1.0
      ELSE
      V(I,N+1)=0.0
      ENDIF
      ENDDO
      V(2:M+1,1)=0
*>SYMMETRY AXIS V(1,2:N)=V(2,2:N);WALL V(M+2,2:N)=0.0*2.0-V(M+1,2:N)
      V(M+2,2:N)=0.0*2.0-V(M+1,2:N)
      V(1,2:N)=0.0*2.0-V(2,2:N)
*>      V(1,2:N)=V(2,2:N)
*>outlet the velocity in couple with condensation
       
      RETURN
      END SUBROUTINE BOUNDARYV

      SUBROUTINE BOUNDARYP
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      P(2:M+1,1)=P(2:M+1,2)
      DO I=2,M+1
      IF(0.5*(XC(I,N+1)+XC(I-1,N+1)).LE.0.3)THEN
*>      P(I,N+2)=P(I,N+1)+(X(I,N+2)-X(I,N+1))*V(I,N+1)
      P(I,N+2)=P(1,N+1)
      ELSE
      P(I,N+2)=P(I,N+1)
      ENDIF
      ENDDO
      DO J=2,N+1
      IF(0.5*(YC(1,J)+YC(1,J-1)).LE.0.3)THEN
*>      P(1,J)=P(2,J)-(X(J,J)-X(J-1,J))*U(1,J)
      P(1,J)=P(2,J)
      ELSE
      P(1,J)=P(2,J)
      ENDIF
      ENDDO
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
      BUOY=GR/RE/RE*0.5*(TEM(I,J)+TEM(I,J+1))
*c      BUOY=0
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
      REAL(KIND=8)::BDA(N),BDB(N),BDC(M),BDD(M),F(M,N)
c      REAL(KIND=8)::BDA(103),BDB(103),BDC(98),BDD(98),F(98,103)
      INTEGER(KIND=4)::MBDCND,NBDCND,IDIMF,IERROR
      REAL(KIND=8)::ELMBDA,PERTRB
      REAL(KIND=8)::W(5*M*N)
*>U,MBDCND:1,NBDCND=1
*>USING HSTCRT SUBROUTINE
      MBDCND=1;NBDCND=1;IDIMF=M
*>VELOCITY IN COUPLE WITH CONDENSATION
*>      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0.0
      DO J=1,N
      IF(0.5*(YC(1,J+1)+YC(1,J)).LE.0.3)THEN
      BDA(J)=1.0
      ELSE
      BDA(J)=0.0
      ENDIF
      ENDDO
      BDB(:)=0.0;BDC(:)=0.0;BDD(:)=0.0
      ELMBDA=-RE/DT 
      DO J=1,N
      DO I=1,M
      F(I,J)=-RE/DT*0.5*(U(I,J+1)+U(I+1,J+1))
      ENDDO
      ENDDO
      CALL HSTCRT(A,B,M,MBDCND,BDA,BDB,
     1            C,D,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,*)'SOLVE U IERROR',IERROR
*>TRANSFORM F(I,J) TO U
      DO J=1,N
      DO I=1,M-1
      U(I+1,J+1)=0.5*(F(I,J)+F(I+1,J))
      ENDDO
      ENDDO
      CALL BOUNDARYU
*>V
      MBDCND=1;NBDCND=1;IDIMF=M
*>      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0.0
      BDA(:)=0.0;BDB(:)=0.0;BDC(:)=0.0
      DO I=1,M
      IF(0.5*(XC(I,N+1)+XC(I+1,N+1)).LE.0.3)THEN
      BDD(I)=1.0
      ELSE
      BDD(I)=0.0
      ENDIF
      ENDDO
      ELMBDA=-RE/DT 
      DO J=1,N
      DO I=1,M
      F(I,J)=-RE/DT*0.5*(V(I+1,J+1)+V(I+1,J))
      ENDDO
      ENDDO
      CALL HSTCRT(A,B,M,MBDCND,BDA,BDB,
     1            C,D,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,*)'SOLVE V IERROR',IERROR
*>TRANSFORM F(I,J) TO V
      DO J=1,N-1
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
      CALL HSTCRT(A,B,M,MBDCND,BDA,BDB,
     1            C,D,N,NBDCND,BDC,BDD,
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

      SUBROUTINE THIRDSOLVEUV
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      DO J=2,N+1
      DO I=2,M
      U(I,J)=U(I,J)-DT/DX*(P(I+1,J)-P(I,J))
      ENDDO
      ENDDO
      CALL BOUNDARYU
      DO J=2,N
      DO I=2,M+1
      V(I,J)=V(I,J)-DT/DY*(P(I,J+1)-P(I,J))
      ENDDO
      ENDDO
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE THIRDSOLVEUV

      SUBROUTINE MODIFYUVP
*>I THINK THERE IS NO NEED TO DISTINGUISH POINT TYPES(INTERIOR OR BOUNDARY POINT)
      IMPLICIT NONE
      LOGICAL,EXTERNAL::ISCONS
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD,DIV,DP
      OPEN(10,FILE='DIVDPPUV.DAT')
      WRITE(10,'(A)')'DIV,DP,P(I,J),U(I,J),U(I-1,J),V(I,J),V(I,J-1)'
      DO J=2,N+1
      DO I=2,M+1
      RAD=X(I,J)
      DIV=(U(I,J)-U(I-1,J))/DX+(V(I,J)-V(I,J-1))/DY+
     1       (U(I,J)+U(I-1,J))/2/RAD*MODEKSI 
      DP=-DIV/2/DT/(1/DX/DX+1/DY/DY)
      P(I,J)=P(I,J)+DP
      IF(I.GT.2.AND.I.LT.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.N+1)THEN
      DP=-DIV/2/DT/(1/DX/DX+1/DY/DY)
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
      ELSEIF(I.EQ.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
      ELSEIF(I.EQ.M+1.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
      ELSEIF(I.EQ.2.AND.J.EQ.N+1)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
      ELSEIF(I.EQ.M+1.AND.J.EQ.N+1)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
      ENDIF
      WRITE(10,'(7ES15.5)')DIV,DP,P(I,J),U(I,J),U(I-1,J),V(I,J),V(I,J-1)
      ENDDO
      ENDDO
      CLOSE(10)
      CALL BOUNDARYU
      CALL BOUNDARYV
      CALL BOUNDARYP
      RETURN
      END SUBROUTINE MODIFYUVP

      SUBROUTINE JUDGE 
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD,DIV
      MAXDIV=1.0E-30
      DO J=2,N+1
      DO I=2,M+1
      RAD=X(I,J)
      DIV=ABS((U(I,J)-U(I-1,J))/DX+(V(I,J)-V(I,J-1))/DY+
     1       (U(I,J)+U(I-1,J))/2/RAD*MODEKSI)
      IF(DIV.GT.MAXDIV)THEN
      MAXDIV=DIV
      ENDIF
      ENDDO
      ENDDO
      IF(MAXDIV.GT.RESIERR)THEN
      ISCONS=.FALSE.
      ELSE
      ISCONS=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE JUDGE

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
      CALL THIRDSOLVEUV
      CALL MODIFYUVP
      NCOUNT=0
      CALL JUDGE
      WRITE(*,*)ISCONS
      WRITE(*,80)'ISCONS=',ISCONS,'ITER=',NCOUNT,
     1           'MAXDIV=',MAXDIV
      DO WHILE(.NOT.ISCONS)
      CALL MODIFYUVP
      CALL JUDGE
      NCOUNT=NCOUNT+1
      IF(MOD(NCOUNT,100).EQ.0)THEN
      WRITE(*,80)'ISCONS=',ISCONS,'ITER=',NCOUNT,
     1           'MAXDIV=',MAXDIV
      ENDIF
80    FORMAT(1X,A10,L4,A10,I6,A10,ES15.5)
      ENDDO
      RETURN
      END SUBROUTINE PROJSOLA

      SUBROUTINE OUTPUT(NUM)
*>VELOCITY AND PRESSURE ...USEING UT,VT,PT...*******************
*>THESE VARIABLES ARE INTERPLOTED TO THE GRID POINT(XC,YC)******
*>(QST(1:M+1,1:N+1));(RHO(1:M+2,1:N+2),U(0:M+2,1:N+2),V(1:M+2,0:N+2))
*>(P(1:M+2,1:N+2));UT(0:M+2,1:N+2),VT(1:M+2,0:N+2)
      IMPLICIT NONE
      INTEGER(KIND=4),INTENT(IN)::NUM
      INTEGER(KIND=4)::I,J
      CHARACTER(LEN=20)::FILE1
      REAL(KIND=8)::UT(1:M+1,1:N+1),VT(1:M+1,1:N+1),
     1              PT(1:M+1,1:N+1)
      WRITE(FILE1,'(I6)')NUM
      FILE1='0'//TRIM(ADJUSTL(FILE1))//'.DAT'
*>STREAMLINE FUNCTION
      CALL STREAMLINE
*>VELOCITY
      UT=0;VT=0;PT=0;UT(:,N+1)=1
      DO J=2,N
      DO I=2,M
      UT(I,J)=0.5*(U(I,J)+U(I,J+1))
      VT(I,J)=0.5*(V(I,J)+V(I+1,J))
      ENDDO
      ENDDO 
      DO J=1,N+1
      DO I=1,M+1
      PT(I,J)=0.25*(P(I,J)+P(I+1,J)+P(I,J+1)+P(I+1,J+1))
      ENDDO
      ENDDO 
      PT(1,1)=1.0/3*(P(2,2)+P(1,2)+P(2,1))
      PT(M+1,1)=1.0/3*(P(M+1,2)+P(M+2,2)+P(M+1,1))
      PT(M+1,N+1)=1.0/3*(P(M+1,N+1)+P(M+2,N+1)+P(M+1,N+2))
      PT(1,N+1)=1.0/3*(P(1,N+1)+P(2,N+1)+P(2,N+2)) 
      OPEN(10,FILE=FILE1)
      WRITE(10,'(A)')'VARIABLES="X","Y","T","STRL","U","V","P"'
      WRITE(10,100)M+1,N+1
      DO J=1,N+1
      DO I=1,M+1
      WRITE(10,101)XC(I,J),YC(I,J),TEMWRITE(I,J),QST(I,J),UT(I,J),
     1             VT(I,J),PT(I,J)
      ENDDO
      ENDDO
      CLOSE(10)
100   FORMAT(1X,'ZONE T="FLOW",I=',I4,' , J=',I4,' , F=POINT')
101   FORMAT(1X,3F10.5,4ES20.5)
      RETURN
      END SUBROUTINE OUTPUT

      SUBROUTINE STREAMLINE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(M+1,N+1)
      INTEGER(KIND=4)::MBDCND,NBDCND,IDIMF,IERROR
      REAL(KIND=8)::ELMBDA,PERTRB
      REAL(KIND=8)::W(5*M*N)
*>STREAM LINE,MBDCND:3,NBDCND=3
*>USING HWSCRT SUBROUTINE
      MBDCND=1;NBDCND=1;IDIMF=M+1
*>      BDA(:)=0;BDB(:)=0;BDC(:)=0;BDD(:)=0
      DO J=1,N+1
      IF(YC(1,J).LE.0.3)THEN
      BDA(J)=1.0*(YC(1,J)-0.3)
      ELSE
      BDA(J)=0.0
      ENDIF
      ENDDO
      BDB(:)=0.0;BDC(:)=0.0
      DO I=1,M+1
      IF(XC(I,N+1).LE.0.3)THEN
      BDD(I)=-1.0*(XC(I,N+1)-0.3)
      ELSE
      BDD(I)=0.0
      ENDIF
      ENDDO
      ELMBDA=0.0
      F=0
      DO J=2,N
      DO I=2,M
      F(I,J)=(U(I,J+1)-U(I,J))/DY-
     1       (V(I+1,J)-V(I,J))/DX 
      ENDDO
      ENDDO
      CALL HWSCRT(A,B,M,MBDCND,BDA,BDB,
     1            C,D,N,NBDCND,BDC,BDD,
     1            ELMBDA,F,IDIMF,PERTRB,IERROR,W) 
      WRITE(*,'(A30,I4)')'SOLVE STREAM FUNCTION ERROR=',IERROR
*>TRANSFORM F(I,J) TO QST 
      DO J=1,N+1
      DO I=1,M+1
      QST(I,J)=F(I,J)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE STREAMLINE
      END MODULE PHYPROJFISH
      
      PROGRAM MAIN
      USE PHYPROJFISH
      USE USERFUNC
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J,NOUT
      REAL(KIND=8)::TEND=1.0,T=0.0
      INTEGER(KIND=4)::CONTROL=0
      A=0.0;B=1.0;C=0.0;D=1.0
      M=80;N=80;NOUT=1
      DX=(B-A)/M;DY=(D-C)/N
      CALL CYLINDER2DMESH
      WRITE(*,20)M,N,DX,DY
20    FORMAT('MESH GRID M*N',2I6,1x,'DX DY',2F8.4)
      CALL INITIAL
      WRITE(*,60)'TIME=',T,'INITIAL DT=',DT
      DO WHILE(T.LE.TEND)
*>      CALL CFL
      T=T+DT
      WRITE(*,60)'TIME=',T,'CFL DT=',DT
      CALL PROJSOLA
      IF(T.GE.(0.001*NOUT))THEN
      CALL OUTPUT(NOUT)
      NOUT=NOUT+1
      WRITE(*,50)'OUTPUT NOUT=',NOUT,'TIME=',T
50    FORMAT(1X,A15,I6,1X,A5,1X,F8.4)
      ENDIF
      ENDDO
60    FORMAT(1X,A10,F10.5,1X,A10,F10.5)
61    FORMAT(1X,A10,F10.5,'WRITE FILE',1X,A10,I10)
      STOP
      END PROGRAM MAIN

