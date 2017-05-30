C
C =====================================================================
      SUBROUTINE SFQ4(KSI, ETA, N, DNK, DNE)
C
C     Evaluate the shape functions and their derivatives wrt the natural
C     coordinates of ISOQ4 at the Gauss Point (KSI, ETA)
C
C     .. Scalar Arguments ..
C     INTEGER*4        KSI    : Value of the horizontal N.C
C                      ETA    : Value of the vertical N.C
C     ..
C     .. Array Arguments ..
C     REAL*8           N(4)   : Shape functions
C                      DNK(4) : Derivaties of S.F. wrt KSI
C                      DNE(4) : Derivaties of S.F. wrt ETA
C     ..
C     .. Scalar Arguments ..
      REAL*8           KSI, ETA
C     ..
C     .. Array Arguments ..
      REAL*8           N(4), DNK(4), DNE(4)
C     ..
C     .. Executable statements ..
C
      N(1) = 0.25D0*(1.D0 - KSI)*(1.D0 - ETA)
      N(2) = 0.25D0*(1.D0 + KSI)*(1.D0 - ETA)
      N(3) = 0.25D0*(1.D0 + KSI)*(1.D0 + ETA)
      N(4) = 0.25D0*(1.D0 - KSI)*(1.D0 + ETA)
C
      DNK(1) = -0.25D0*(1.D0 - ETA)
      DNK(2) =  0.25D0*(1.D0 - ETA)
      DNK(3) =  0.25D0*(1.D0 + ETA)
      DNK(4) = -0.25D0*(1.D0 + ETA)
C
      DNE(1) = -0.25D0*(1.D0 - KSI)
      DNE(2) = -0.25D0*(1.D0 + KSI)
      DNE(3) =  0.25D0*(1.D0 + KSI)
      DNE(4) =  0.25D0*(1.D0 - KSI)
C
      RETURN
C
C     .. End of SFQ4 ..
C
      END
