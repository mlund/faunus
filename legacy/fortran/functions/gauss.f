************************************************************************
*                                                                      *
************************************************************************
      SUBROUTINE GAUSS(N,LDIM,A,X,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(LDIM,LDIM),X(LDIM),C(LDIM)
      DO 10 I=1,N
         X(I)=C(I)
10    CONTINUE
      DO 100 I=1,N-1
         K=I
         DO 110 J=I+1,N
            IF(ABS(A(K,I)) .LT. ABS(A(J,I))) K=J
110      CONTINUE
         IF(K .NE. I) THEN
C            WRITE(*,'(A,2I3)') ' SWAPPING:',I,K
            DO 120 J=I,N
               SWAP=A(I,J)
               A(I,J)=A(K,J)
               A(K,J)=SWAP
120         CONTINUE
            SWAP=X(I)
            X(I)=X(K)
            X(K)=SWAP
         END IF
         DO 130 K=I+1,N
            FACT=A(K,I)/A(I,I)
            DO 131 J=I+1,N
               A(K,J)=A(K,J)-FACT*A(I,J)
131         CONTINUE
            X(K)=X(K)-FACT*X(I)
130      CONTINUE
100   CONTINUE
      X(N)=X(N)/A(N,N)
      DO 200 I=N-1,1,-1
         DO 210 K=I+1,N
            X(I)=X(I)-A(I,K)*X(K)
210      CONTINUE
         X(I)=X(I)/A(I,I)
200   CONTINUE
      RETURN
      END
