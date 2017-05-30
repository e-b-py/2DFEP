C
C =====================================================================
      SUBROUTINE T3(ELCOR, CMAT, D, KL, VLDS, LDVLDS, FL)
C
C     Element routine for the Constant Strain Triangle (T3)
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDVLDS
C     ..
C     .. Array Arguments ..
C     INTEGER*4    
C
C     REAL*8           ELCOR(6, 2)    : Matrix of element nodal coordinates
C                      CMAT(4)        : Vector of material properties
C                      D(3, 3)        : Material stiffness matrix
C                      KL(12, 12)     : Local stiffness matrix
C                      VLDS(LDVLDS, *): Matrix of volume loads
C                      FL(12)         : Local force matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        
C                      
C     REAL*8           GAMM: Unit weight
C                      A   : Area of the element
C                      B1, B2, B3, C1, C2, C3, D11, D12, D21, D22, D33
C     ..
C     .. Local Arrays ..
C     INTEGER*4        PRM(6)         : Permutation vector
C     
C     REAL*8           COEFF(3, 3)    : Matrix that stores a_i, b_i, c_i    
C     ..
C     .. Common Scalars ..
C     INTEGER*4        NNODE: number of nodes
C                      NELE : number of elements
C                      NRBEL: number of restrained boundary elements
C                      NLBEL: number of loaded boundary elements
C                      NRNOD: number of restrained points
C                      NLNOD: number of loaded points
C                      NERST: number of element restraints
C                      NPRST: number of point restraints
C                      NDOF : number of degrees of freedom
C                      NNDL : number of nodal loads
C                      NSL  : number of surface loads
C                      NVL  : number of volume loads
C                      ETYPE: element type specifier
C                      PTYPE: problem type specifier
C                      NBCO : number of coord. for boundary elements
C                      NECO : number of coord. for domain elements
C     ,,
C     .. Scalar Arguments ..
      INTEGER*4        LDVLDS
C     ..
C     .. Array Arguments ..
      REAL*8           ELCOR(6, 2), CMAT(4), VLDS(LDVLDS, *),
     ;                 KL(12, 12), FL(12), D(3,3)
C     ..
C =====================================================================
C     .. Local Scalars ..
      REAL*8           GAMM, A, B1, B2, B3, C1, C2, C3,
     ;                 D11, D12, D21, D22, D33
C     ..
C     .. Local Arrays ..
      INTEGER*4        PRM(6)
      REAL*8           COEFF(3, 3)
C     .. Common Scalars ..
      INTEGER*4        NNODE, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                 NPRST, NDOF, NNDL, NSL, NVL, ETYPE, PTYPE, NBCO,
     ;                 NECO
      COMMON /CONFIG/  NNODE, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                 NPRST, NDOF, NNDL, NSL, NVL, ETYPE, PTYPE, NBCO,
     ;                 NECO
C     ..
C     .. Executable statements ..
C
C     Compute the coefficients a_i, b_i, c_i, {i = 1, 2, 3}
C
      WRITE(30, '(A)') 'ELCOR vector (Element coordinates)'
      DO 901 K = 1, NECO
         WRITE(30, *)
         DO 902 J = 1, 2
            WRITE(30, '(F8.2)', ADVANCE='NO') ELCOR(K, J)
  902    CONTINUE
  901 CONTINUE
      PRM = (/1, 2, 3, 1, 2, 3/)
c      PRM = (/3, 2, 1, 3, 2, 1/) ? 
      DO 10 J = 1, 3
         COEFF(J, 1) = ELCOR(PRM(J+1), 1)*ELCOR(PRM(J+2), 2)
     ;                 - ELCOR(PRM(J+1), 2)*ELCOR(PRM(J+2), 1)
         COEFF(J, 2) = ELCOR(PRM(J+1), 2) - ELCOR(PRM(J+2), 2)
         COEFF(J, 3) = ELCOR(PRM(J+2), 1) - ELCOR(PRM(J+1), 1)
   10 CONTINUE
      WRITE(92) COEFF
      WRITE(30, '(//A)') 'Element coefficients'
      WRITE(30, '(/A, F10.4)') 'A1 =', COEFF(1, 1)
      WRITE(30, '(A, F10.4)') 'A2 =', COEFF(2, 1)
      WRITE(30, '(A, F10.4)') 'A3 =', COEFF(3, 1)
      WRITE(30, '(A, F10.4)') 'B1 =', COEFF(1, 2)
      WRITE(30, '(A, F10.4)') 'B2 =', COEFF(2, 2)
      WRITE(30, '(A, F10.4)') 'B3 =', COEFF(3, 2)
      WRITE(30, '(A, F10.4)') 'C1 =', COEFF(1, 3)
      WRITE(30, '(A, F10.4)') 'C2 =', COEFF(2, 3)
      WRITE(30, '(A, F10.4)') 'C3 =', COEFF(3, 3)
C
C     Compute the area of the element
C
      A =  ( ELCOR(2, 1)*ELCOR(3, 2) - ELCOR(3, 1)*ELCOR(2, 2) )
     ;   + ( ELCOR(3, 1)*ELCOR(1, 2) - ELCOR(1, 1)*ELCOR(3, 2) )
     ;   + ( ELCOR(1, 1)*ELCOR(2, 2) - ELCOR(2, 1)*ELCOR(1, 2) )
      A = DABS(A) / 2.D0
      WRITE(92) A
      WRITE(30, '(/A, F10.4)') 'Element Area, A =', A
      
C
C     Stiffness matrix
C
      B1 = COEFF(1, 2)
      B2 = COEFF(2, 2)
      B3 = COEFF(3, 2)
      C1 = COEFF(1, 3)
      C2 = COEFF(2, 3)
      C3 = COEFF(3, 3)
      D11 = D(1, 1)
      D12 = D(1, 2)
      D21 = D(2, 1)
      D22 = D(2, 2)
      D33 = D(3, 3)
      KL(1, 1) = D11*B1**2 + D33*C1**2
      KL(1, 2) = ( D12 + D33 )*B1*C1
      KL(2, 2) = D22*C1**2 + D33*B1**2
      KL(1, 3) = D11*B1*B2 + D33*C1*C2
      KL(2, 3) = D12*B2*C1 + D33*B1*C2
      KL(3, 3) = D11*B2**2 + D33*C2**2
      KL(1, 4) = D12*B1*C2 + D33*B2*C1
      KL(2, 4) = D22*C1*C2 + D33*B1*B2
      KL(3, 4) = ( D12 + D33 ) *B2*C2
      KL(4, 4) = D22*C2**2 + D33*B2**2
      KL(1, 5) = D11*B1*B3 + D33*C1*C3 
      KL(2, 5) = D21*B3*C1 + D33*B1*C3 
      KL(3, 5) = D11*B2*B3 + D33*C2*C3 
      KL(4, 5) = D21*B3*C2 + D33*B2*C3 
      KL(5, 5) = D11*B3**2 + D33*C3**2 
      KL(1, 6) = D12*B1*C3 + D33*B3*C1 
      KL(2, 6) = D22*C1*C3 + D33*B1*B3 
      KL(3, 6) = D12*B2*C3 + D33*B3*C2
      KL(4, 6) = D22*C2*C3 + D33*B2*B3
      KL(5, 6) = ( D12 + D33 )*B3*C3
      KL(6, 6) = D22*C3**2 + D33*B3**2
C     Symmetry
      KL = KL * CMAT(4) / (4*A)
      DO 20 J = 1, 6
         DO 30 K = J+1, 6
            KL(K, J) = KL(J, K)
   30    CONTINUE
   20 CONTINUE
C
C     Compute the volume loads and self weight
C
      GAMM = CMAT(3)
      DO 40 J = 1, NVL
         DO 50 K = 1, 5, 2
            FL(K) = FL(K) + VLDS(J, 1)
            FL(K+1) = FL(K+1) + VLDS(J, 2) - GAMM
   50    CONTINUE
   40 CONTINUE
      FL = FL*A*CMAT(4)/3.D0
      WRITE(30, '(/A)') 'Volume Loads'
      WRITE(30, '(/A, F10.4)') 'FX = ', FL(1)
      WRITE(30, '(A, F10.4)') 'FY = ', FL(2)
      WRITE(30, '(/)')
      RETURN
C
C     .. End of T3 ..
C
      END
