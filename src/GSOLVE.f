C
C =====================================================================
      SUBROUTINE GSOLVE(A, LDA, M, N, X, B)
C
C     Solve the system of linear equations AX=B using the Gauss
C     Elimination algorithm.
C
C     .. Scalar Arguments ..
C     INTEGER*4         LDA: leading dimension of the LHS
C                         M: number of rows 
C                         N: number of columns
C     ..
C     .. Array Arguments ..
C     REAL*8      A(LDA, *): LHS matrix
C                      X(*): Vector of unknowns
C                      B(*): RHS vector
C     ..
C     .. Local Scalars ..
C     REAL*8          RATIO: Ratio of eliminated equation's element 
C                            corresponding to the eliminated unknown
C                            to the pivot element
C                     PIVOT: Pivot element 
C     ,,
      INTEGER*4        LDA, M, N
      REAL*8           A(LDA, *), X(*), B(*)
C     ..
C =====================================================================
C     .. Local Scalars ..
      REAL*8           RATIO, PIVOT
C     ..
C     .. Executable statements ..
C
C     Forward elimination 
C
      DO 10 I = 1, M-1
         PIVOT = A(I,I)
         B(I) = B(I) / PIVOT
         DO 20 K = I, N
            A(I,K) = A(I,K) / PIVOT
   20    CONTINUE
         DO 30 J = I+1, M
            IF( DABS( A(J,I) ).LT.( 1.D-14 ) ) GO TO 30
            RATIO = A(J,I) / A(I,I)
            B(J) = B(J) - RATIO*B(I)
            DO 40 K = I, N
               A(J,K) = A(J,K) - RATIO*A(I,K)
   40       CONTINUE
   30    CONTINUE
   10 CONTINUE
      B(M) = B(M)/A(M,M)
      A(M,M)=1.D0
C
C     Back Substitution
C
      X(M) = B(M) / A(M,M)
      DO 50 I = M-1, 1, -1
         X(I) = B(I)
         DO 60 J = I+1, N
            X(I) = X(I) - A(I,J)*X(J)
   60    CONTINUE
         X(I) = X(I) / A(I,I)
   50 CONTINUE
      RETURN
C      
C     .. End of GSOLVE ..
C      
      END 


       
