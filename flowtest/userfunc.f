      MODULE USERFUNC

      CONTAINS 
*>****************USER SUBROUTINE****************************
      SUBROUTINE  PLU(A,B,N)
*>1.FIND THE COLUMN OF THE MAXIMUM IN EACH ROW;
*>2.TRANSFORM THIS COLUMN TO THAT OF THE CURRENT ROW VALUE.
*>3.LU FRACTORIZATION
*>4.WHEN SOLVING EQUS:AX=B-->(PA`)X=B-->P(LU)X=B--->LUX=P-1B
*>5.I=II=(IC1C2...CNRN...R2R1I)=P(RN..R2R1I)=PP-1
*>6.LUX=P-1B=(RN...R2R1I)B
*>7.B:DO ROW OPERATION CORRESPONDING TO P'S COLUMN OPERATIONS.
      IMPLICIT NONE
      INTEGER(KIND=4),INTENT(IN)::N
      REAL(KIND=8)::A(1:N,1:N),B(1:N,1)
      INTEGER(KIND=4)::I,J,K
      REAL(KIND=8)::MAXVALUE,TMP(1:N,1),BTMP(1:N,1)
      REAL(KIND=8)::P(1:N,1:N),L(1:N,1:N),U(1:N,1:N)
      INTEGER(KIND=4)::MAXINDEX
      BTMP=B
      P=0.0;L=0.0;U=0.0
      DO I=1,N
      P(I,I)=1.0;L(I,I)=1.0
      ENDDO
*>P MATRIX:PERMUTATION MATRIX
*>I ROW <----> MAXINDEX ROW----->ROW OPERATION
      DO I=1,N
      MAXVALUE=1E-6
      MAXINDEX=1
      DO J=I,N
      IF(ABS(A(J,I)).GT.MAXVALUE)THEN
      MAXVALUE=A(J,I)
      MAXINDEX=J
      ENDIF
      ENDDO
      DO K=1,N
      TMP(K,1)=A(I,K)
      A(I,K)=A(MAXINDEX,K)
      A(MAXINDEX,K)=TMP(K,1)
      TMP(K,1)=P(K,I)
*>P DO COLUMN OPERATION RELATED TO A ROW OPERATION
      P(K,I)=P(K,MAXINDEX)
      P(K,MAXINDEX)=TMP(K,1)
      ENDDO
*>B ONLY NEEDED TO BE DONE ROW OPERATION
      TMP(1,1)=B(I,1)
      B(I,1)=B(MAXINDEX,1)
      B(MAXINDEX,1)=TMP(1,1)
      ENDDO
*>LU FRACTORIZATION
      U(1,:)=A(1,:)
      L(:,1)=A(:,1)/U(1,1)
      DO I=2,N
*>U MATRIX UPPER TRIANGULAR
      DO J=I,N
      TMP(1,1)=0.0
      DO K=1,I-1
      TMP(1,1)=TMP(1,1)+L(I,K)*U(K,J)
      ENDDO
      U(I,J)=A(I,J)-TMP(1,1)
      ENDDO
*>L MATRIX LOWER TRIANGULAR
      DO J=I+1,N
      TMP(1,1)=0
      DO K=1,I-1
      TMP(1,1)=TMP(1,1)+L(J,K)*U(K,I)
      ENDDO
      L(J,I)=(A(J,I)-TMP(1,1))/U(I,I)
      ENDDO
      ENDDO
*>LUX=B,LY=B,UX=Y
*>SOLVE Y
      B(1,1)=B(1,1)/L(1,1)
      DO I=2,N
      TMP(1,1)=0
      DO J=1,I-1
      TMP(1,1)=TMP(1,1)+L(I,J)*B(J,1)
      ENDDO 
      B(I,1)=B(I,1)-TMP(1,1) 
      ENDDO
*>SOLVE X
      B(N,1)=B(N,1)/U(N,N)
      DO I=1,N-1
      TMP(1,1)=0
      DO J=1,I
      TMP(1,1)=TMP(1,1)+U(N-I,N-I+J)*B(N-I+J,1)
      ENDDO
      B(N-I,1)=(B(N-I,1)-TMP(1,1))/U(N-I,N-I)
      ENDDO
      RETURN
      END SUBROUTINE PLU


      SUBROUTINE SPLINE(X,Y,N,X0,Y0,N0)
*>RETURN A ARRAY NOT A VALUE*******IMPORTANT POINT**********
*>X,Y KNOWN POINT DATAS,X0 INTERPOLATION POINT,Y0 UNKONWN
*>SOLVE NATURAL BOUNDARY CONDITION**************************
      IMPLICIT NONE
*>INPUT AND OUTPUT PARAMETERS
*>NOT USED EXTERNAL IN THE SAME MODULE
*c      EXTERNAL:: PLU
      INTEGER(KIND=4),INTENT(IN)::N,N0
      REAL(KIND=8),INTENT(IN)::X(N),Y(N),X0(N0)
      REAL(KIND=8),INTENT(OUT)::Y0(N0)
*>INTERNAL VARIABLES
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::H(N-1),M(N),A(N,N),B(N)
      REAL(KIND=8)::CO_1,CO_2,CO_3,CO_4
      INTEGER(KIND=4)::INDEX(N0)
*>SOLVE SECOND-ORDER DERIVATIVE MATRIX
      DO I=1,N-1
      H(I)=X(I+1)-X(I)
      ENDDO
      DO I=1,N
      DO J=1,N
      A(I,J)=0
      ENDDO
      B(I)=0
      M(I)=0
      ENDDO
      DO I=2,N-1
      A(I,I)=2*(H(I-1)+H(I))
      A(I,I-1)=H(I-1)
      A(I,I+1)=H(I)
      B(I)=6*((Y(I+1)-Y(I))/H(I)-(Y(I)-Y(I-1))/H(I-1))
      ENDDO

c      WRITE(*,'(7F6.2)') (H(I),I=1,7)
c      WRITE(*,'(7F6.2)') (B(I),I=1,7)
*>**************************************************************      
      A(1,1)=1
      A(N,N)=1
      B(1)=0
      B(N)=0
*>***********LAPACK*********************************************       
      CALL PLU(A,B,N)
      M=B
*>INDEX(:) HELP SELECTING WHICH SEGMENT FUNCTION
      DO I=1,N0
      INDEX(I)=1
      DO J=2,N
      IF(X0(I).GT.X(J))THEN
      INDEX(I)=INDEX(I)+1
      ELSE
      CYCLE
      ENDIF
      ENDDO
*>************THIRD-ORDER INTERPOLATION*************************
*>******Y0=CO_1+CO_2*(X0-X)+CO_3*(X0-X)**2+CO_4*(X0-X)**3*******
      CO_4=(M(INDEX(I)+1)-M(INDEX(I)))/6/H(INDEX(I))
      CO_2=(Y(INDEX(I)+1)-Y(INDEX(I)))/H(INDEX(I))
     1   -H(INDEX(I))/2*M(INDEX(I))
     1   -H(INDEX(I))/6*(M(INDEX(I)+1)-M(INDEX(I)))
      CO_3=M(INDEX(I))/2
      CO_1=Y(INDEX(I))
      Y0(I)=CO_1+CO_2*(X0(I)-X(INDEX(I)))
     1     +CO_3*(X0(I)-X(INDEX(I)))**2
     1     +CO_4*(X0(I)-X(INDEX(I)))**3
      ENDDO

      RETURN
      END SUBROUTINE SPLINE 
*>********************USER FUNCTION******************************** 
      FUNCTION SPLINT(X,Y,N,X0)
*>RETURN A VALUE CORRESPONDING TO X0************************
*>X,Y KNOWN POINT DATAS,X0 INTERPOLATION POINT
*>SOLVE NATURAL BOUNDARY CONDITION**************************
      IMPLICIT NONE
*>INPUT AND OUTPUT PARAMETERS
*>NOT USED EXTERNAL IN THE SAME MODULE
*c      EXTERNAL:: PLU
      REAL(KIND=8)::SPLINT
      INTEGER(KIND=4),INTENT(IN)::N
      REAL(KIND=8),INTENT(IN)::X(N),Y(N)
      REAL(KIND=8),INTENT(IN)::X0
*>INTERNAL VARIABLES
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::H(N-1),M(N),A(N,N),B(N)
      REAL(KIND=8)::CO_1,CO_2,CO_3,CO_4
      INTEGER(KIND=4)::INDEX
*>SOLVE SECOND-ORDER DERIVATIVE MATRIX
      DO I=1,N-1
      H(I)=X(I+1)-X(I)
      ENDDO
      DO I=1,N
      DO J=1,N
      A(I,J)=0
      ENDDO
      B(I)=0
      M(I)=0
      ENDDO
      DO I=2,N-1
      A(I,I)=2*(H(I-1)+H(I))
      A(I,I-1)=H(I-1)
      A(I,I+1)=H(I)
      B(I)=6*((Y(I+1)-Y(I))/H(I)-(Y(I)-Y(I-1))/H(I-1))
      ENDDO

c      WRITE(*,'(7F6.2)') (H(I),I=1,7)
c      WRITE(*,'(7F6.2)') (B(I),I=1,7)
*>*NATURAL BOUNDARY:THE 2-ORDER DERIVATION AT ENDPOINTS EQUALS 0*      
      A(1,1)=1
      A(N,N)=1
      B(1)=0
      B(N)=0
*>**************************************************************       
      CALL PLU(A,B,N)
      M=B
*>INDEX HELP SELECTING WHICH SEGMENT FUNCTION
      INDEX=1
      DO J=2,N
      IF(X0.GT.X(J))THEN
      INDEX=INDEX+1
      ELSE
      CYCLE
      ENDIF
      ENDDO
*>************THIRD-ORDER INTERPOLATION*************************
*>******Y0=CO_1+CO_2*(X0-X)+CO_3*(X0-X)**2+CO_4*(X0-X)**3*******
      CO_4=(M(INDEX+1)-M(INDEX))/6/H(INDEX)
      CO_2=(Y(INDEX+1)-Y(INDEX))/H(INDEX)
     1   -H(INDEX)/2*M(INDEX)
     1   -H(INDEX)/6*(M(INDEX+1)-M(INDEX))
      CO_3=M(INDEX)/2
      CO_1=Y(INDEX)
      SPLINT=CO_1+CO_2*(X0-X(INDEX))
     1     +CO_3*(X0-X(INDEX))**2
     1     +CO_4*(X0-X(INDEX))**3

      RETURN
      END FUNCTION SPLINT 

      FUNCTION MAXVAL(A,X1,X2,Y1,Y2)
      IMPLICIT NONE
      REAL(KIND=8)::MAXVAL
      INTEGER(KIND=4),INTENT(IN)::X1,X2,Y1,Y2
      REAL(KIND=8)::A(*)
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::REF=1.0E-30
      DO I=1,(X2-X1+1)*(Y2-Y1+1)
      IF(A(I).GT.REF)REF=A(I)
      ENDDO
      MAXVAL=REF
      RETURN
      END FUNCTION MAXVAL

      FUNCTION MINVAL(A,X1,X2,Y1,Y2)
      IMPLICIT NONE
      REAL(KIND=8)::MINVAL
      INTEGER(KIND=4),INTENT(IN)::X1,X2,Y1,Y2
      REAL(KIND=8)::A(*)
      INTEGER(KIND=4)::I,J
      REAL(KIND=8)::REF=1.0E+30
      DO I=1,(X2-X1+1)*(Y2-Y1+1)
      IF(A(I).LT.REF)REF=A(I)
      ENDDO
      MINVAL=REF
      RETURN
      END FUNCTION MINVAL

      END MODULE USERFUNC
