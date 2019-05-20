      PROGRAM TESTMATRIX
      USE USERFUNC
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J
      INTEGER(KIND=4)::M=4,N=3
      REAL(KIND=8),ALLOCATABLE::A(:),B(:,:)
      REAL(KIND=8),ALLOCATABLE::X0(:),Y0(:),X1(:),Y1(:)
      ALLOCATE(A(M*M),X0(M*M),Y0(M*M),B(M,N),X1(M),Y1(N))
      DO I=1,M*M
      X0(I)=I/2.0+0.5;Y0(I)=I/2.0+0.5
      A(I)=SIN(I/2.0+0.5)
      ENDDO
      WRITE(*,*)'*******WRITE X0******'
      WRITE(*,'(4F10.5)')X0
      WRITE(*,*)'*******WRITE Y0******'
      WRITE(*,'(4F10.5)')Y0
      WRITE(*,*)'*******WRITE A*******'
      WRITE(*,'(4F10.5)')A
      X1(1:M)=X0(1:M)/2.0+0.2;Y1=Y0(1:N)/2.0+0.1 
      WRITE(*,*)'*******WRITE X1*******'
      WRITE(*,'(4F10.5)')X1
      WRITE(*,*)'*******WRITE Y1*******'
      WRITE(*,'(4F10.5)')Y1
      CALL MATRIXTOMATRIX(A,X0,Y0,M,M,B,X1,Y1,M,N)
      WRITE(*,'(4F10.5)')B
      STOP
      END PROGRAM TESTMATRIX
