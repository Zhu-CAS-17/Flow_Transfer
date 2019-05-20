      PROGRAM MAIN
      USE PHYSOLVE
      IMPLICIT NONE
      INTEGER(KIND=4)::I,J,NOUT=1
      REAL(KIND=8)::TEND=1.0,T=0.0
      INTEGER(KIND=4)::CONTROL=0
      REAL(KIND=8)::X0=0.0,X1=1.0,Y0=0.0,Y1=1.0
*******************************************
******************************************
      M=100;N=100
      CALL CYLINDER2DMESH(X0,X1,Y0,Y1)
      WRITE(*,60)'TIME=',T,'INTERVAL=',DT
      CALL INITIAL
      DO WHILE(T.LE.TEND)
      T=T+DT
      WRITE(*,60)'TIME=',T,'INTERVAL=',DT
      CALL PROJSOLA
      STOP
      CALL CFL
      IF(T.GE.(0.005*NOUT))THEN
      CALL OUTPUT(NOUT)
      NOUT=NOUT+1
      WRITE(*,61)'TIME=',T,'NOUT=',NOUT
      ENDIF
      ENDDO
60    FORMAT(1X,A10,F10.5,1X,A10,F10.5)
61    FORMAT(1X,A10,F10.5,'WRITE FILE',1X,A10,I10)

      STOP
      END PROGRAM MAIN
