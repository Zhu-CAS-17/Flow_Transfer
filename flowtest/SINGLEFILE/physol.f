      MODULE PHYSOLVE
      USE MESH
      USE USERFUNC
      IMPLICIT NONE
      REAL(KIND=8),ALLOCATABLE::RHO(:,:),U(:,:),V(:,:)
      REAL(KIND=8),ALLOCATABLE::P(:,:),TEM(:,:),QST(:,:)
*>PROJSOLA CONTROL PARAMETERS****************************************
      INTEGER(KIND=4)::MODEKSI=0,SKIP=1
      REAL(KIND=8)::ALPHA=0.0,OMIC
      REAL(KIND=8)::DT=1.0E-5,DX,DY
      REAL(KIND=8)::RE=100,GR=0.0
      REAL(KIND=8)::RESIERR=1.0E-4,MAXDIV
      CONTAINS
*>***************BOUNDARY SETTINGS*************************************
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
      IMPLICIT NONE
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
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD
      INTEGER(KIND=8)::NODE(5)
      REAL(KIND=8)::UMAT((M-1)*N,(M-1)*N),VMAT(M*(N-1),M*(N-1))
      REAL(KIND=8)::BU(1:(M-1)*N),BV(1:M*(N-1))
      INTEGER(KIND=4)::IPIV(M*N),INFO
*>IMPLICIT VISCOSITY TERM FOR INTERMIDDLE VELOCITY
      UMAT=0
*>U
      DO J=2,N+1
      DO I=2,M
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=(I-1)+(M-1)*(J-2)
      NODE(2)=(I-1-1)+(M-1)*(J-2)
      NODE(3)=(I+1-1)+(M-1)*(J-2)
      NODE(4)=(I-1)+(M-1)*(J-1-2)
      NODE(5)=(I-1)+(M-1)*(J+1-2)
      DX=X(I+1,J)-X(I,J)
      DY=YC(I,J)-YC(I,J-1)
      RAD=0.5*(X(I,J)+X(I+1,J))
*      VISX=1/RE*(1/DX/DX*(U(I+1,J)-2*U(I,J)+U(I-1,J))+
*     1           1/DY/DY*(U(I,J+1)-2*U(I,J)+U(I,J-1))+
*     1           MODEKSI/2.0/DX/RAD*(U(I+1,J)-U(I-1,J))-
*     1           MODEKSI*U(I,J)/RAD/RAD)
      IF(I.GT.2.AND.I.LT.M.AND.J.GT.2.AND.J.LT.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BU(NODE(1))=U(I,J)
*>SOUTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M.AND.J.EQ.2)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
*>NORTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M.AND.J.EQ.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BU(NODE(1))=U(I,J)
*>EAST WALL
      ELSEIF(I.EQ.M.AND.J.GT.2.AND.J.LT.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BU(NODE(1))=U(I,J)
*>SW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
*>SE POINT
      ELSEIF(I.EQ.M.AND.J.EQ.2)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
*>NE POINT
      ELSEIF(I.EQ.M.AND.J.EQ.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
*>NW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.N+1)THEN
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)-
     1                     (-DT/RE*(1/DY/DY))
      BU(NODE(1))=U(I,J)
      ENDIF
      ENDDO
      ENDDO
      WRITE(*,*)'SECOND U CALL PLU'
      CALL PLU(UMAT,BU,(M-1)*N)
      WRITE(*,*)MAXVAL(BU,1,(M-1)*N,1,1)
*>TRANSFORM BU TO U
      DO J=2,N+1
      DO I=2,M
*>NODE1:(I,J)
      NODE(1)=(I-1)+(M-1)*(J-2)
      U(I,J)=BU(NODE(1))
      ENDDO
      ENDDO
      CALL BOUNDARYU
*>V
      DO J=2,N
      DO I=2,M+1
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=(I-1)+M*(J-2)
      NODE(2)=(I-1-1)+M*(J-2)
      NODE(3)=(I+1-1)+M*(J-2)
      NODE(4)=(I-1)+M*(J-1-2)
      NODE(5)=(I-1)+M*(J+1-2)
      DX=XC(I,J)-XC(I-1,J)
      DY=Y(I,J+1)-Y(I,J)
      RAD=0.5*(XC(I,J)+XC(I-1,J)) 
*      VISY=1/RE*(1/DX/DX*(V(I+1,J)-2*V(I,J)+V(I-1,J))+
*     1           1/DY/DY*(V(I,J+1)-2*V(I,J)+V(I,J-1))+
*     1           MODEKSI/DX/RAD/2.0*(V(I+1,J)-V(I-1,J)))
      IF(I.GT.2.AND.I.LT.M+1.AND.J.GT.2.AND.J.LT.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>SOUTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.2)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>NORTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX))
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>EAST WALL
      ELSEIF(I.EQ.M+1.AND.J.GT.2.AND.J.LT.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX))
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>SW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX))
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>SE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.2)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX))
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>NE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX))
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>NW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.N)THEN
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)-
     1                     (-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX))
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
      ENDIF
      ENDDO
      ENDDO
      CALL PLU(VMAT,BV,M*(N-1))
*>TRANSFORM BV TO V 
      DO J=2,N
      DO I=2,M+1
      NODE(1)=(I-1)+M*(J-2)
      V(I,J)=BV(NODE(1))
      ENDDO
      ENDDO
*>BOUNDARY V
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE SECONDSOLVEUV

      SUBROUTINE SOLVEPRESSURE
*>*************SOLVE FIRST STEP VELOCITY*******************************
*>*************SOLVE THE NEXT INTERVAL PRESSURE************************
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::NODE(5)
      REAL(KIND=8)::RAD
      REAL(KIND=8)::PMAT(M*N,M*N),BP(1:M*N)
      INTEGER(KIND=4)::IPIV(M*N),INFO
*>*************SOLVE THE NEXT INTERVAL PRESSURE************************
      DO J=2,N+1
      DO I=2,M+1
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=(I-1)+M*(J-2)
      NODE(2)=(I-1-1)+M*(J-2)
      NODE(3)=(I+1-1)+M*(J-2)
      NODE(4)=(I-1)+M*(J-1-2)
      NODE(5)=(I-1)+M*(J+1-2)
      DX=HKSI(I,J)
      DY=HETA(I,J)
      RAD=X(i,j)
      IF(I.GT.2.AND.I.LT.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(4))=1/DY/DY
      PMAT(NODE(1),NODE(5))=1/DY/DY
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>SOUTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.2)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(5))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>NORTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(4))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX-MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(4))=1/DY/DY
      PMAT(NODE(1),NODE(5))=1/DY/DY
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>EAST WALL
      ELSEIF(I.EQ.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX+MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(4))=1/DY/DY
      PMAT(NODE(1),NODE(5))=1/DY/DY
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>SW POINT
      ELSEIF(I.GT.2.AND.J.EQ.2)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX-MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(5))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>SE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.2)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX+MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(5))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>NE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX+MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(4))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>NW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.N+1)THEN
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX+
     1                     (1/DX/DX-MODEKSI*1/RAD/2/DX)
      PMAT(NODE(1),NODE(4))=1/DY/DY+
     1                     (1/DY/DY)
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
      ENDIF
      ENDDO
      ENDDO
      CALL PLU(PMAT,BP,(M-1)*(N-1))
      DO J=2,N+1
      DO I=2,M+1
      NODE(1)=(I-1)+M*(J-2)
      P(I,J)=BP(NODE(1))
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
      IF(I.GT.2.AND.I.LT.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
*>SOUTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
*>NORTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M+1.AND.J.EQ.N+1)THEN
      U(I,J)=U(I,J)+DT*DP/DX
      U(I-1,J)=U(I-1,J)-DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
*>EAST WALL
      ELSEIF(I.EQ.M+1.AND.J.GT.2.AND.J.LT.N+1)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+DT*DP/DY
      V(I,J-1)=V(I,J-1)-DT*DP/DY
*>SW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
*>SE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+2.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-0.0*DT*DP/DY
*>NE POINT
      ELSEIF(I.EQ.M+1.AND.J.EQ.N+1)THEN
      U(I,J)=U(I,J)+0.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-2.0*DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
*>NW POIONT
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      U(I,J)=U(I,J)+2.0*DT*DP/DX
      U(I-1,J)=U(I-1,J)-0.0*DT*DP/DX
      V(I,J)=V(I,J)+0.0*DT*DP/DY
      V(I,J-1)=V(I,J-1)-2.0*DT*DP/DY
      ENDIF
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
      WRITE(*,*)'SECONDSOLVE'
      CALL SOLVEPRESSURE
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

      SUBROUTINE STREAMLINE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::NODE(5)
      REAL(KIND=8)::RAD
      REAL(KIND=8)::QSTMAT((M-1)*(N-1),(M-1)*(N-1))
      REAL(KIND=8)::BQST(1:(M-1)*(N-1))
      INTEGER(KIND=4)::IPIV(M*N),INFO
      QSTMAT=0;BQST=0
      DO J=2,N
      DO I=2,M 
      DX=0.5*((XC(I,J)+XC(I+1,J))-(XC(I-1,J)+XC(I,J)))
      DY=0.5*((YC(I,J)+YC(I,J+1))-(YC(I,J-1)+YC(I,J)))
      RAD=XC(I,J)
      NODE(1)=I-1+(M-1)*(J-2)
      NODE(2)=I-1-1+(M-1)*(J-2)
      NODE(3)=I+1-1+(M-1)*(J-2)
      NODE(4)=I-1+(M-1)*(J-1-2)
      NODE(5)=I-1+(M-1)*(J+1-2)
      IF(I.GT.2.AND.I.LT.M.AND.J.GT.2.AND.J.LT.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>SOUTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M.AND.J.EQ.2)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>NORTH WALL
      ELSEIF(I.GT.2.AND.I.LT.M.AND.J.EQ.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GT.2.AND.J.LT.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>EAST WALL
      ELSEIF(I.EQ.M.AND.J.GT.2.AND.J.LT.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>SW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.2)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>SE POINT
      ELSEIF(I.EQ.M.AND.J.EQ.2)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(5))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>NE POINT
      ELSEIF(I.EQ.M.AND.J.EQ.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
*>NW POINT
      ELSEIF(I.EQ.2.AND.J.EQ.N)THEN
      QSTMAT(NODE(1),NODE(1))=-2.0*(1/DX/DX+1/DY/DY)
      QSTMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      QSTMAT(NODE(1),NODE(4))=1/DY/DY
      IF(MODEKSI.EQ.1)THEN
      BQST(NODE(1))=-2.0*0.5*(V(I,J)+V(I+1,J))+RAD*
     1                    ((U(I,J+1)-U(I,J))/DY-
     1                     (V(I+1,J)-V(I,J))/DX)
      ELSE
      BQST(NODE(1))=(U(I,J+1)-U(I,J))/DY-(V(I+1,J)-V(I,J))/DX
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      CALL PLU(QSTMAT,BQST,(M-1)*(N-1))
      DO J=2,N
      DO I=2,M
      NODE(1)=I-1+(M-1)*(J-2)
      QST(I,J)=BQST(NODE(1))
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE STREAMLINE

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

      SUBROUTINE CFL
      IMPLICIT NONE
      REAL(KIND=8)::UMAX,VMAX
      DX=X(3,2)-X(2,2);DY=Y(2,3)-Y(2,2)
      UMAX=MAXVAL(U,0,M+2,1,N+2)
      VMAX=MAXVAL(V,1,M+2,0,N+2)
      DT=0.25*MIN(DX/UMAX,DY/VMAX)
      ALPHA=1.2*MAX(ABS(UMAX*DT/DX),ABS(VMAX*DT/DY))
      WRITE(*,200)DT,ALPHA
200   FORMAT(1X,'DT=',F8.4,1X,'ALPHA=',F8.4,1X)
      RETURN
      END SUBROUTINE CFL
      
      END MODULE PHYSOLVE
