C
C =====================================================================
      SUBROUTINE ISOQ4(ELCOR, CMAT, D, KL, VLDS, LDVLDS, FL)
C
C     Element routine for the 4-Node Isoparametric Quadrilateral
C     (ISOQ4)
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
C     INTEGER*4        NGP   : Number of Guass points for integration
C                      
C     REAL*8           KSI   : Horizontal natural coordinate
C                      ETA   : Vertical natural coordinate
C                      DETJ  : Determinant of the Jacobian
C                      GAMM  : Self weight
C                      THCK  : Element thickness
C                      FX    : Volume force in x-direction
C                      FY    : Volume force in y-direction
C     ..
C     .. Local Arrays ..
C     INTEGER*4        NODE(6)        : Nodes of the element
C     
C     REAL*8           GP(4)          : Vector of Gauss points
C                      W(4)           : Vector of weights of GPs
C                      N(4)           : Shape functions
C                      DNK(4)         : Derivatives of S.F. wrt KSI
C                      DNE(4)         : Derivatives of S.F. wrt ETA
C                      JAC(2, 2)      : Jacobian matrix
C                      INVJAC(2, 2)   : Inverse of the Jacobian
C                      B(3, 8)        : Strain interpolation matrix
C                      KTMP(8, 8)     : Temporary stiffness matrix
C                      FTMPX(4)       : Temporary force matrix (x-dir)
C                      FTMPY(4)       : Temporary force matrix (y-dir)
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
      INTEGER*4        NGP
      REAL*8           KSI, ETA, GAMM, THCK, FX, FY
     ;          
C     ..
C     .. Local Arrays ..
C      INTEGER*4
      REAL*8           GP(4), W(4), N(4), DNK(4), DNE(4), JAC(2, 2),
     ;                 INVJAC(2, 2), B(3, 8), KTMP(8, 8), FTMPX(4),
     ;                 FTMPY(4)
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
      WRITE(28, '(/A//)') 'ISOQ4 DEBUG FILE'
C
C     Get the position of Gauss Points and the weights
C
      WRITE(28, '(A)') 'Element Coordinates'
      DO 901 I = 1, NECO
         WRITE(28, *)
         DO 902 J = 1, 2
            WRITE(28, '(F6.2)', ADVANCE='NO') ELCOR(I, J)
  902    CONTINUE
  901 CONTINUE
      NGP = 2
      CALL GPTS(NGP, GP, W)
      WRITE(28, '(A/)') 'Gauss Points and Weights'
      WRITE(28, '(A, F18.10, F18.10)') 'GP', GP(1), GP(2)
      WRITE(28, '(A, F18.10, F18.10)') 'W', W(1), W(2)
C
C     Loop through the Gauss Points
C
      GAMM = CMAT(3)
      THCK = CMAT(4)
      KL = 0.D0
      FL = 0.D0
      DO 10 I = 1, NGP
         KSI = GP(I)
         DO 20 J = 1, NGP
            ETA = GP(J)
            WRITE(28, '(/A, F14.10, 5X, A, F14.10/)') 
     ;            'KSI = ', KSI, 'ETA = ', ETA
            CALL SFQ4(KSI, ETA, N, DNK, DNE)
            WRITE(28, '(/A, F14.10, F14.10, F14.10, F14.10)')
     ;      'N = ', N(1),N(2),N(3),N(4)
            WRITE(28, '(/A, F14.10, F14.10, F14.10, F14.10)')
     ;      'DNK = ',DNK(1),DNK(2),DNK(3),DNK(4)
            WRITE(28, '(/A, F14.10, F14.10, F14.10, F14.10)')
     ;      'DNE = ',DNE(1),DNE(2),DNE(3),DNE(4)
C
C           Compute Jacobian and its inverse
C
            JAC = 0.D0
            DO 30 K = 1, NECO
               JAC(1, 1) = JAC(1, 1) + DNK(K) * ELCOR(K, 1)
               JAC(1, 2) = JAC(1, 2) + DNK(K) * ELCOR(K, 2)
               JAC(2, 1) = JAC(2, 1) + DNE(K) * ELCOR(K, 1)
               JAC(2, 2) = JAC(2, 2) + DNE(K) * ELCOR(K, 2)
   30       CONTINUE
            WRITE(28, '(/A/)') 'Jacobian matrix'
            WRITE(28, '(F14.10, F14.10)') JAC(1, 1), JAC(1, 2)
            WRITE(28, '(F14.10, F14.10)') JAC(2, 1), JAC(2, 2)
            DETJ = JAC(1, 1)*JAC(2, 2) - JAC(1, 2)*JAC(2, 1)
            WRITE(28, '(/A, F14.10/)') 'Determinant of Jacobian', DETJ
            INVJAC(1, 1) =  JAC(2, 2) / DETJ
            INVJAC(1, 2) = -JAC(1, 2) / DETJ
            INVJAC(2, 1) = -JAC(2, 1) / DETJ
            INVJAC(2, 2) =  JAC(1, 1) / DETJ
            WRITE(28, '(/A/)') 'Inverse of the Jacobian matrix'
            WRITE(28, '(F14.10, F14.10)') INVJAC(1, 1), INVJAC(1, 2)
            WRITE(28, '(F14.10, F14.10)') INVJAC(2, 1), INVJAC(2, 2)
C
C           Fill the strain interpolation matrix
C
            B = 0.D0
            B(1, 1) = INVJAC(1, 1)*DNK(1) + INVJAC(1, 2)*DNE(1)
            B(2, 2) = INVJAC(2, 1)*DNK(1) + INVJAC(2, 2)*DNE(1)
            B(3, 1) = B(2, 2)
            B(3, 2) = B(1, 1)
            B(1, 3) = INVJAC(1, 1)*DNK(2) + INVJAC(1, 2)*DNE(2)
            B(2, 4) = INVJAC(2, 1)*DNK(2) + INVJAC(2, 2)*DNE(2)
            B(3, 3) = B(2, 4)
            B(3, 4) = B(1, 3)
            B(1, 5) = INVJAC(1, 1)*DNK(3) + INVJAC(1, 2)*DNE(3)
            B(2, 6) = INVJAC(2, 1)*DNK(3) + INVJAC(2, 2)*DNE(3)
            B(3, 5) = B(2, 6)
            B(3, 6) = B(1, 5)
            B(1, 7) = INVJAC(1, 1)*DNK(4) + INVJAC(1, 2)*DNE(4)
            B(2, 8) = INVJAC(2, 1)*DNK(4) + INVJAC(2, 2)*DNE(4)
            B(3, 7) = B(2, 8)
            B(3, 8) = B(1, 7)
            WRITE(28, '(/A)') 'B matrix'
            DO 903 K = 1, 3
               WRITE(28, *)
               DO 904 L = 1, 8
                  WRITE(28, '(F6.2)', ADVANCE='NO') B(K, L)
  904          CONTINUE
  903       CONTINUE
                  
C
C           Compute  (B_T*D*B*detJ) and add to the stiffness matrix
C
            KTMP = MATMUL(TRANSPOSE(B), MATMUL(D, B))*DETJ*W(I)*W(J)
            WRITE(28, '(/A)') 'KTMP'
            DO 905 K = 1, 8
               WRITE(28, *)
               DO 906 L = 1, 8
                  WRITE(28, '(E14.4)', ADVANCE='NO') KTMP(K, L)
  906          CONTINUE
  905       CONTINUE

            DO 40 K = 1, 8
               DO 50 L = 1, 8
                  KL(K, L) = KL(K, L) + KTMP(K, L)
   50          CONTINUE
   40       CONTINUE
C
C           Integrate volume forces
C
            DO 60 K = 1, NVL
               FX = VLDS(K, 1)
               FY = VLDS(K, 2) - GAMM
               FTMPX = N*FX*DETJ*W(I)*W(J)
               FTMPY = N*FY*DETJ*W(I)*W(J)
               DO 70 L = 1, 4
                  IODD  = 2*L - 1
                  IEVEN = 2*L
                  FL(IODD)  = FL(IODD)  + FTMPX(L)
                  FL(IEVEN) = FL(IEVEN) + FTMPY(L)
   70          CONTINUE
   60       CONTINUE
C
   20    CONTINUE
   10 CONTINUE
      KL = KL*THCK
      FL = FL*THCK
      RETURN
C
C     .. ISOQ4 ..
C
      END
