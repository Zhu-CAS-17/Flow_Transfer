      MODULE MESH
      IMPLICIT NONE
      INTEGER(KIND=4)::M,N
*>XC,YC:ZONE GRID NODES;X,Y:COMPUTING DATA NODES
*>HKSI,HETA,VOL,AE,AK:THE DOMAIN VOLUME'S LENGTH,HEIGHT,VOLUME,
*>LEFT SIDE SURFACE AREA AND BOTTOM SURFACE AREA SEPARATELY.
      REAL(KIND=8),ALLOCATABLE::XC(:,:),YC(:,:),X(:,:),Y(:,:)
      REAL(KIND=8),ALLOCATABLE::HKSI(:,:),HETA(:,:),VOL(:,:)
      REAL(KIND=8),ALLOCATABLE::AK(:,:),AE(:,:)

      CONTAINS 

      SUBROUTINE CYLINDER2DMESH(X0,X1,Y0,Y1)
*>USED IN FLOW SIMULATING PROCESS***************************
*>ONE ADDITIONAL VIRTUAL GRID*******************************
*>X0<X1,Y0<Y1;**********************************************
*>X--->R,Y--->Z;(X,Y) RELATE TO P POINT*********************
*>FUTURE,ZONE GRID NODES CAN BE MODIFIED********************
*>AXIAL-SYSMETRY CYLINDERIAL COORDINATES********************
      IMPLICIT NONE
      REAL(KIND=8),INTENT(IN)::X0,X1,Y0,Y1
      REAL(KIND=8),PARAMETER::PI=3.1415926
      INTEGER(KIND=4)::I,J
      ALLOCATE(XC(0:M+2,0:N+2),YC(0:M+2,0:N+2))
      ALLOCATE(X(1:M+2,1:N+2),Y(1:M+2,1:N+2))
      ALLOCATE(HKSI(1:M+2,1:N+2),HETA(1:M+2,1:N+2))
      ALLOCATE(VOL(1:M+2,1:N+2),AK(1:M+2,1:N+2),AE(1:M+2,1:N+2))

      DO I=0,M+2
      XC(I,:)=X0+(I-1)*(X1-X0)/M
      ENDDO
      DO J=0,N+2
      YC(:,J)=Y0+(J-1)*(Y1-Y0)/N
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
      
      END MODULE MESH
