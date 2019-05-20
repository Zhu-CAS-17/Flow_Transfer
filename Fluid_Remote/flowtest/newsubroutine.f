      SUBROUTINE SECONDSOLVEUV
*>IMPLICIT VISCOSITY TERM FOR INTERMIDDLE VELOCITY
*>SOLVING UMAT MATRIX(WHICH CONCLUDES ALL MESH POINTS)
*>SO THE DIMENSION OF UMAT EQUALS THAT OF U
*>ONLY THE DOMAIN PIVOT ELEMENT OF UMAT ARE SETTED TO 1,
*>THE COORESPOING BU EQUALS THE BOUNDARY VALUE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::RAD
      INTEGER(KIND=8)::NODE(5)
      REAL(KIND=8)::UMAT((M+3)*(N+2),(M+3)*(N+2))
      REAL(KIND=8)::VMAT((M+2)*(N+3),(M+2)*(N+3))
      REAL(KIND=8)::BU(1:(M+3)*(N+2)),BV(1:(M+2)*(N+3))
      INTEGER(KIND=4)::IPIV((M+3)*(N+3)),INFO
      UMAT=0;VMAT=0;BU=0;BV=0
      DO I=1,(M+3)*(N+2)
      UMAT(I,I)=1
      ENDDO
      DO I=1,(M+2)*(N+3)
      VMAT(I,I)=1
      ENDDO
*>U
      DO J=2,N+1
      DO I=2,M
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=I+1+(M+3)*(J-1)
      NODE(2)=I-1+1+(M+3)*(J-1)
      NODE(3)=I+1+1+(M+3)*(J-1)
      NODE(4)=I+1+(M+3)*(J-1-1)
      NODE(5)=I+1+(M+3)*(J+1-1)
      DX=X(I+1,J)-X(I,J)
      DY=YC(I,J)-YC(I,J-1)
      RAD=0.5*(X(I,J)+X(I+1,J))
*      VISX=1/RE*(1/DX/DX*(U(I+1,J)-2*U(I,J)+U(I-1,J))+
*     1           1/DY/DY*(U(I,J+1)-2*U(I,J)+U(I,J-1))+
*     1           MODEKSI/2.0/DX/RAD*(U(I+1,J)-U(I-1,J))-
*     1           MODEKSI*U(I,J)/RAD/RAD)
      UMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY-MODEKSI/RAD/RAD)
      UMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      UMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      UMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BU(NODE(1))=U(I,J)
*>SOUTH WALL
      IF(I.GE.2.AND.I.LE.M.AND.J.EQ.2)THEN
      NODE(1)=I+1+(M+3)*(J-2)
      NODE(5)=I+1+(M+3)*(J-2+1)
      U(NODE(1),NODE(5))=1
      BU(NODE(1))=0
*>NORTH WALL
      ELSEIF(I.GE.2.AND.I.LE.M.AND.J.EQ.N+1)THEN
      NODE(1)=I+1+(M+3)*J
      NODE(4)=I+1+(M+3)*(J-1)
      U(NODE(1),NODE(4))=1
      BU(NODE(1))=2.0*1.0
      ENDIF
*>WEST WALL AND EAST WALL SATISFY THE BOUNDARY REALATION
      ENDDO
      ENDDO
      WRITE(*,*)'SECOND U CALL PLU'
      CALL PLU(UMAT,BU,(M+3)*(N+2))
*>TRANSFORM BU TO U
      DO J=2,N+1
      DO I=2,M
*>NODE1:(I,J)
      NODE(1)=(I+1)+(M+3)*(J-1)
      U(I,J)=BU(NODE(1))
      ENDDO
      ENDDO
      CALL BOUNDARYU
*>V
      DO J=2,N
      DO I=2,M+1
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=I+(M+2)*J
      NODE(2)=I-1+(M+2)*J
      NODE(3)=I+1+(M+2)*J
      NODE(4)=I+(M+2)*(J-1)
      NODE(5)=I+(M+2)*(J+1)
      DX=XC(I,J)-XC(I-1,J)
      DY=Y(I,J+1)-Y(I,J)
      RAD=0.5*(XC(I,J)+XC(I-1,J)) 
*      VISY=1/RE*(1/DX/DX*(V(I+1,J)-2*V(I,J)+V(I-1,J))+
*     1           1/DY/DY*(V(I,J+1)-2*V(I,J)+V(I,J-1))+
*     1           MODEKSI/DX/RAD/2.0*(V(I+1,J)-V(I-1,J)))
      VMAT(NODE(1),NODE(1))=1-DT/RE*(-2/DX/DX-2/DY/DY)
      VMAT(NODE(1),NODE(2))=-DT/RE*(1/DX/DX-MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(3))=-DT/RE*(1/DX/DX+MODEKSI/RAD/2/DX)
      VMAT(NODE(1),NODE(4))=-DT/RE*(1/DY/DY)
      VMAT(NODE(1),NODE(5))=-DT/RE*(1/DY/DY)
      BV(NODE(1))=V(I,J)
*>SOUTH WALL AND NORTH WALL SATISIFY THE BOUNDARY CONDITION NATURALLY
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GE.2.AND.J.LE.N)THEN
      NODE(1)=(I-1)+(M+2)*J
      NODE(3)=(I+1-1)+(M+2)*J
      VMAT(NODE(1),NODE(3))=1
      BV(NODE(1))=0
*>EAST WALL
      ELSEIF(I.EQ.M+1.AND.J.GE.2.AND.J.LE.N)THEN
      NODE(1)=I+1+(M+2)*J
      NODE(2)=I+(M+2)*J
      VMAT(NODE(1),NODE(2))=1
      BV(NODE(1))=0
      ENDIF
      ENDDO
      ENDDO
      CALL PLU(VMAT,BV,(M+2)*(N+3))
*>TRANSFORM BV TO V 
      DO J=2,N
      DO I=2,M+1
      NODE(1)=I+(M+2)*J
      V(I,J)=BV(NODE(1))
      ENDDO
      ENDDO
*>BOUNDARY V
      CALL BOUNDARYV
      RETURN
      END SUBROUTINE SECONDSOLVEUV

      SUBROUTINE SOLVEPRESSURE
*>*************USE FIRST STEP VELOCITY*******************************
*>*************SOLVE THE NEXT INTERVAL PRESSURE************************
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::NODE(5)
      REAL(KIND=8)::RAD
      REAL(KIND=8)::PMAT((M+1)*(N+1),(M+1)*(N+1))
      REAL(KIND=8)::BP(1:(M+1)*(N+1))
      INTEGER(KIND=4)::IPIV((M+1)*(N+1)),INFO
      PMAT=0;BP=0
      DO I=1,(M+1)*(N+1)
      PMAT(I,I)=1
      ENDDO
      DO J=2,N+1
      DO I=2,M+1
*>NODE1:(I,J)2:(I-1,J)3:(I+1,J)4:(I,J-1)5:(I,J+1)
      NODE(1)=I+(M+2)*(J-1)
      NODE(2)=I-1+(M+2)*(J-1)
      NODE(3)=I+1+(M+2)*(J-1)
      NODE(4)=I+(M+2)*(J-1-1)
      NODE(5)=I+(M+2)*(J+1-1)
      DX=HKSI(I,J)
      DY=HETA(I,J)
      RAD=X(i,j)
      PMAT(NODE(1),NODE(1))=-2/DX/DX-2/DX/DX
      PMAT(NODE(1),NODE(2))=1/DX/DX-MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(3))=1/DX/DX+MODEKSI*1/RAD/2/DX
      PMAT(NODE(1),NODE(4))=1/DY/DY
      PMAT(NODE(1),NODE(5))=1/DY/DY
      BP(NODE(1))=-1/DT*((U(I,J)-U(I-1,J))/DX+
     1           MODEKSI*(U(I,J)+U(I-1,J))/2/RAD+
     1                   (V(I,J)-V(I,J-1))/DY)
*>SOUTH WALL
      IF(I.GE.2.AND.I.LE.M+1.AND.J.EQ.2)THEN
      NODE(1)=I+(M+2)*(J-2)
      NODE(5)=I+(M+2)*(J-2+1)
      PMAT(NODE(1),NODE(5))=-1
      BP(NODE(1))=0
*>NORTH WALL
      ELSEIF(I.GE.2.AND.I.LE.M+1.AND.J.EQ.N+1)THEN
      NODE(1)=I+(M+2)*J
      NODE(4)=I+(M+2)*(J-1)
      PMAT(NODE(1),NODE(4))=-1
      BP(NODE(1))=0
*>WEST WALL
      ELSEIF(I.EQ.2.AND.J.GE.2.AND.J.LE.N+1)THEN
      NODE(1)=I-1+(M+2)*(J-1)
      NODE(3)=I+(M+2)*(J-1)
      PMAT(NODE(1),NODE(3))=-1
      BP(NODE(1))=0
*>EAST WALL
      ELSEIF(I.EQ.M+1.AND.J.GE.2.AND.J.LE.N+1)THEN
      NODE(1)=I+1+(M+2)*(J-1)
      NODE(2)=I+(M+2)*(J-1)
      PMAT(NODE(1),NODE(2))=-1
      BP(NODE(1))=0
      ENDIF
      ENDDO
      ENDDO
      CALL PLU(PMAT,BP,(M+1)*(N+1))
      DO J=2,N+1
      DO I=2,M+1
      NODE(1)=I+(M+2)*(J-1)
      P(I,J)=BP(NODE(1))
      ENDDO
      ENDDO
*>BOUNDARY P
      CALL BOUNDARYP
      RETURN
      END SUBROUTINE SOLVEPRESSURE

      SUBROUTINE STREAMLINE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::NODE(5)
      REAL(KIND=8)::RAD
      REAL(KIND=8)::QSTMAT((M+1)*(N+1),(M+1)*(N+1))
      REAL(KIND=8)::BQST(1:(M+1)*(N+1))
      INTEGER(KIND=4)::IPIV(M*N),INFO
      QSTMAT=0;BQST=0
      DO I=1,(M+1)*(N+1)
      QSTMAT(I,I)=1
      ENDDO
      DO J=2,N
      DO I=2,M 
      DX=0.5*((XC(I,J)+XC(I+1,J))-(XC(I-1,J)+XC(I,J)))
      DY=0.5*((YC(I,J)+YC(I,J+1))-(YC(I,J-1)+YC(I,J)))
      RAD=XC(I,J)
      NODE(1)=I+(M+1)*(J-1)
      NODE(2)=I-1+(M+2)*(J-1)
      NODE(3)=I+1+(M+2)*(J-1)
      NODE(4)=I+(M+2)*(J-1-1)
      NODE(5)=I+(M+2)*(J+1-1)
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
*>SOUTH NORTH WEST EAST WALL:QST=0
      ENDDO
      ENDDO
      CALL PLU(QSTMAT,BQST,(M+1)*(N+1))
      DO J=2,N
      DO I=2,M
      NODE(1)=I+(M+1)*(J-1)
      QST(I,J)=BQST(NODE(1))
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE STREAMLINE

