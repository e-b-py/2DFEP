C
C =====================================================================
      SUBROUTINE GPTS(NGP, GP, W)
C
C     Return the Gauss points and their corresponding weights 
C
C     .. Scalar Arguments ..
C     INTEGER*4        NGP    : Number of Gauss points
C     ..
C     .. Array Arguments ..
C     REAL*8           GP(4)  : Vector of Gauss points
C                      W(4)   : Vector of weights of GPs
C     ..
C     .. Scalar Arguments ..
      INTEGER*4        NGP
C     ..
C     .. Array Arguments ..
      REAL*8           GP(4), W(4)
C     ..
C     .. Executable statements ..
C
      IF ( NGP.EQ.1 ) GO TO 1
      IF ( NGP.EQ.2 ) GO TO 2
      IF ( NGP.EQ.3 ) GO TO 3
      IF ( NGP.EQ.4 ) GO TO 4
C
    1 GP(1) =  0.D0
      W(1)  =  2.D0
      GO TO 100
C
    2 GP(1) = -0.5773502691896257D0
      GP(2) =  0.5773502691896257D0
      W(1)  =  1.D0
      W(2)  =  1.D0
      GO TO 100
C
    3 GP(1) =  0.D0
      GP(2) = -0.7745966692414834D0
      GP(3) =  0.7745966692414834D0
      W(1)  =  0.8888888888888888D0
      W(2)  =  0.5555555555555556D0
      W(3)  =  0.5555555555555556D0
      GO TO 100
C
    4 GP(1) = -0.3399810435848563D0
      GP(2) =  0.3399810435848563D0
      GP(3) = -0.8611363115940526D0
      GP(4) =  0.8611363115940526D0
      W(1)  =  0.6521451548625461D0
      W(2)  =  0.6521451548625461D0
      W(3)  =  0.3478548451374538D0
      W(4)  =  0.3478548451374538D0
  100 RETURN
C
C     .. End of GPTS ..
C
      END
