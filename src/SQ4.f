C =====================================================================
      SUBROUTINE SQ4(CON, LDCON, CRD, LDCRD, UCRD, IDOF, LDIDOF, 
     ;               DSPG, D)
C
C     Compute sig_xx, sig_yy, sig_xy, principal min and max stresses, 
C     their orientations and Von Mises stresses for ISOQ4 elements
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCON, LDCRD, LDIDOF
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Connectivity matrix
C                      IDOF(LDIDOF, *): Matrix of the degrees of freedom
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of coordinates
C                      UCRD(NNODE, 2) : Matrix of updated coordinates
C                      DSPG(*)        : Global nodal displacement vector
C                      D(3, 3)        : Material stiffness matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        NGP      : Number of Gauss points
C                      NLDOF    : Number of local DOFs
C                      CNODE    : Current node during GP loop
C
C     REAL*8           A        : Element area
C                      SGX      : Sigma_xx
C                      SGY      : Sigma_yy
C                      SGXY     : Sigma_xy
C                      SGVM     : Von Mises stress
C                      SGMAX    : Maximum principle stress
C                      SGMIN    : Minimum principle stress
C                      ORMAX    : Orientation of SGMAX
C                      ORMIN    : Orientation of SGMIN
C                      KSI      : Value of the horizontal N.C
C                      ETA      : Value of the vertical N.C
C                      DETJ     : Determinant of the Jacobian
C     ..
C     .. Local Arrays ..
C     INTEGER*4        ELDOF(8)  : Global DOFs corresponding to locals 
C                      NODE(4)   : Nodes of the element
C
C     REAL*8           DSPL(8)        : Local displacements of the el.
C                      ELCOR(4, 2)    : Nodal coordinates of the element
C                      GP(4)          : Vector of Gauss points
C                      W(4)           : Vector of weights of GPs
C                      N(4)           : Shape functions
C                      DNK(4)         : Derivatives of S.F wrt KSI
C                      DNE(4)         : Derivatives of S.F. wrt ETA
C                      B(3, 8)        : Strain interpolation matrix
C                      JAC(2, 2)      : Jacobian matrix
C                      INVJAC(2, 2)   : Inverse of the Jacobian
C                      SIG(NNODE, 4)  : Stress matrix                        
C                      SIGE(4, 3)     : Element stresses at each node
C            
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
      INTEGER*4        LDCON, LDIDOF, LDCRD
C     ..
C     .. Array Arguments ..
      INTEGER*4        CON(LDCON, *), IDOF(LDIDOF, *)
      REAL*8           CRD(LDCRD, *), DSPG(*), D(3, 3)
C     ..
C =====================================================================
      PARAMETER        ( SQRT3 = 1.7320508075688772D0 )
C     .. Local Scalars ..
      INTEGER*4        NGP, NLDOF, CNODE
      REAL*8           A, SGX, SGY, SGXY, SGVM, SGMAX, SGMIN, ORMAX,
     ;                 ORMIN, KSI, ETA, DETJ
     ;                 
C     ..
C     .. Local Arrays ..
      INTEGER*4        ELDOF(8), NODE(4)
      REAL*8           DSPL(8), GP(4), W(4), N(4), DNK(4), DNE(4),
     ;                 B(3, 8), JAC(2, 2), INVJAC(2, 2), SIG(NNODE, 4),
     ;                 SIGGP(4, 3), SIGNP(4, 3), ELCOR(4, 2),
     ;                 UCRD(NNODE, 2)
C     ..
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
      OPEN(UNIT=98, FILE='../debug/stress_iq4.dat')
      WRITE(98, '(/A//)') 'ISOQ4 STRESS RECOVERY - DEBUG' 
C
C     Prepare the vtk files
C
      DO 10 I = 48, 62, 2
         WRITE(I, '(A)') '# vtk DataFile Version 2.0'
         WRITE(I, '(A)') 'undeformed mesh'
         WRITE(I, '(A)') 'ASCII'
         WRITE(I, '(A)') 'DATASET UNSTRUCTURED_GRID'
         WRITE(I, '(A, I6, A)') 'POINTS', NNODE, ' float'
         DO 20 J = 1, NNODE
            WRITE(I, '(F22.16, F22.16, F22.16)') 
     ;                  UCRD(J, 1), UCRD(J, 2), 0.D0
   20    CONTINUE
         WRITE(I, '(/A, I6, I6)') 'CELLS', NELE, NELE*(NECO+1)
         DO 30 J = 1, NELE
            WRITE(I, '(I6)', ADVANCE='NO') NECO
            DO 40 K = 1, NECO
               WRITE(I, '(I6)', ADVANCE='NO') CON(J, K)-1
   40       CONTINUE
            WRITE(I, *)
   30    CONTINUE
         WRITE(I, '(/A, I6)') 'CELL_TYPES', NELE
         DO 50 J = 1, NELE
            WRITE(I, '(I6)') 5
   50    CONTINUE
   10 CONTINUE
      WRITE(48, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(48, '(A)') 'SCALARS sig_xx float 1'
      WRITE(48, '(A)') 'LOOKUP_TABLE default'
      WRITE(50, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(50, '(A)') 'SCALARS sig_yy float 1'
      WRITE(50, '(A)') 'LOOKUP_TABLE default'
      WRITE(52, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(52, '(A)') 'SCALARS sig_xy float 1'
      WRITE(52, '(A)') 'LOOKUP_TABLE default'
      WRITE(54, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(54, '(A)') 'SCALARS sig_vm float 1'
      WRITE(54, '(A)') 'LOOKUP_TABLE default'
      WRITE(56, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(56, '(A)') 'SCALARS sig_max float 1'
      WRITE(56, '(A)') 'LOOKUP_TABLE default'
      WRITE(58, '(/A, I6)') 'POINT_DATA', NNODE
      WRITE(58, '(A)') 'SCALARS sig_min float 1'
      WRITE(58, '(A)') 'LOOKUP_TABLE default'
C
C     Retrieve Gauss points and initiate other parameters
C
      NGP = 2
      CALL GPTS(NGP, GP, W)
      NLDOF = NECO*2
      SIG = 0.D0
C
C     <<<    E L E M E N T   L O O P   >>>
C
      DO 60 I = 1, NELE
C
C        Retrieve local displacements and element coordinates
C
         WRITE(98, '(A, I3/)') 'Element', I
         DO 70 M = 1, NECO
            NODE(M) = CON(I, M)
   70    CONTINUE
         WRITE(98, '(A/)') 'Vector: NODE'
         DO 901 M = 1, NECO
            WRITE(98, '(I5)') NODE(M)
  901    CONTINUE
         DO 71 M = 1, NECO
            ELCOR(M, 1) = CRD(NODE(M), 1)
            ELCOR(M, 2) = CRD(NODE(M), 2)
   71    CONTINUE
         WRITE(98, '(/A/)') 'Matrix: ELCOR'
         DO 902 K = 1, NECO
            WRITE(98, '(F6.2, F6.2)') ELCOR(K, 1), ELCOR(K, 2)
  902    CONTINUE
         L = 1
         DO 80 M = 1, NECO
            ELDOF(L)   = IDOF(NODE(M), 1)
            ELDOF(L+1) = IDOF(NODE(M), 2)
            L = L + 2
   80    CONTINUE
         WRITE(98, '(/A/)') 'Vector: ELDOF'
         DO 903 K = 1, NECO
            WRITE(98, '(I5)') ELDOF(K)
  903    CONTINUE
         DSPL  = 0.D0
         DO 90 M = 1, NLDOF
            IF ( ELDOF(M).NE.-1 ) THEN
               DSPL(M) = DSPG(ELDOF(M))
            ENDIF
   90    CONTINUE
         WRITE(98, '(/A/)') 'Vector: DSPL'
         DO 904 K = 1, NECO
            WRITE(98, '(F12.10)') DSPL(K)
  904    CONTINUE
C
C        Reset related parameters before Gauss loop
C
         A     = 0.D0
         SIGGP = 0.D0
         SIGNP = 0.D0
         CNODE = 0
C
C        << G A U S S   L O O P >>
C
         WRITE(98, '(//A/)') 'GAUSS LOOP'
         DO 100 J = 1, NGP
            KSI = GP(J)
            DO 110 K = 1, NGP
               ETA = GP(K)
               WRITE(98, '(A, F12.6, A, F12.6)') 'KSI=',KSI,'ETA=',ETA
               CNODE = CNODE + 1
               CALL SFQ4(KSI, ETA, N, DNK, DNE)
C
C              Compute Jacobian, its determinant (for A)  and inverse
C
               JAC = 0.D0
               DO 120 M = 1, NECO
                  JAC(1, 1) = JAC(1, 1) + DNK(M) * ELCOR(M, 1)
                  JAC(1, 2) = JAC(1, 2) + DNK(M) * ELCOR(M, 2)
                  JAC(2, 1) = JAC(2, 1) + DNE(M) * ELCOR(M, 1)
                  JAC(2, 2) = JAC(2, 2) + DNE(M) * ELCOR(M, 2)
  120          CONTINUE
               DETJ = JAC(1, 1)*JAC(2, 2) - JAC(1, 2)*JAC(2, 1)
               A = A + DETJ
               INVJAC(1, 1) =  JAC(2, 2) / DETJ
               INVJAC(1, 2) = -JAC(1, 2) / DETJ
               INVJAC(2, 1) = -JAC(2, 1) / DETJ
               INVJAC(2, 2) =  JAC(1, 1) / DETJ
               WRITE(98, '(/A/)') 'Jacobian'
               DO 921 M = 1, 2
                  WRITE(98, '(F10.6, F10.6)') JAC(M, 1), JAC(M, 2)
  921          CONTINUE
               WRITE(98, '(/A, F10.6/)') 'Determinant =', DETJ
               WRITE(98, '(/A/)') 'Inverse Jacobian'
               DO 922 M = 1, 2
                  WRITE(98, '(F10.6, F10.6)') INVJAC(M, 1), INVJAC(M, 2)
  922          CONTINUE
C
C              Fill the B matrix
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
               WRITE(98, '(/A/)') 'B matrix'
               DO 923 M = 1, 3
                  WRITE(98, *)
                  DO 924 L = 1, 8
                     WRITE(98, '(F6.2)', ADVANCE='NO') B(M, L)
  924             CONTINUE
  923          CONTINUE

C
C              Compute Stresses
C        
               SIGGP(CNODE, :) = MATMUL( D, MATMUL(B, DSPL) )
               WRITE(98, '(/A, I3/)') 'SIGGP, row:', CNODE
               WRITE(98, '(F16.8, F16.8, F16.8)')
     ;               SIGGP(CNODE, 1), SIGGP(CNODE, 2), SIGGP(CNODE, 3)
  925          CONTINUE
  110       CONTINUE
  100    CONTINUE
         WRITE(98, '(/A/)') 'Matrix: SIGGP'
         DO 911 K = 1, 4
            WRITE(98, *)
            DO 912 L = 1, 3
               WRITE(98, '(F16.10)', ADVANCE='NO') SIGGP(K, L)
  912       CONTINUE
  911    CONTINUE
C
C        Extrapolate stresses to the nodes (numeration: 
C        in GP loop => 1-4-2-3, consider when projecting)
C
C         CALL SFQ4(-SQRT3, -SQRT3, N, DNK, DNE)
         CALL SFQ4(DSQRT(3.D0), -DSQRT(3.D0), N, DNK, DNE)
         SIGNP(1, 1) =  SIGGP(1, 1)*N(1) + SIGGP(2, 1)*N(4)
     ;                + SIGGP(3, 1)*N(2) + SIGGP(4, 1)*N(3)
         SIGNP(1, 2) =  SIGGP(1, 2)*N(1) + SIGGP(2, 2)*N(4)
     ;                + SIGGP(3, 2)*N(2) + SIGGP(4, 2)*N(3)
         SIGNP(1, 3) =  SIGGP(1, 3)*N(1) + SIGGP(2, 3)*N(4)
     ;                + SIGGP(3, 3)*N(2) + SIGGP(4, 3)*N(3)
C         CALL SFQ4(SQRT3, -SQRT3, N, DNK, DNE)
         CALL SFQ4(-DSQRT(3.D0), DSQRT(3.D0), N, DNK, DNE)
         SIGNP(2, 1) =  SIGGP(1, 1)*N(1) + SIGGP(2, 1)*N(4)
     ;                + SIGGP(3, 1)*N(2) + SIGGP(4, 1)*N(3)
         SIGNP(2, 2) =  SIGGP(1, 2)*N(1) + SIGGP(2, 2)*N(4)
     ;                + SIGGP(3, 2)*N(2) + SIGGP(4, 2)*N(3)
         SIGNP(2, 3) =  SIGGP(1, 3)*N(1) + SIGGP(2, 3)*N(4)
     ;                + SIGGP(3, 3)*N(2) + SIGGP(4, 3)*N(3)
C         CALL SFQ4(SQRT3, SQRT3, N, DNK, DNE)
         CALL SFQ4(DSQRT(3.D0), DSQRT(3.D0), N, DNK, DNE)
         SIGNP(3, 1) =  SIGGP(1, 1)*N(1) + SIGGP(2, 1)*N(4)
     ;                + SIGGP(3, 1)*N(2) + SIGGP(4, 1)*N(3)
         SIGNP(3, 2) =  SIGGP(1, 2)*N(1) + SIGGP(2, 2)*N(4)
     ;                + SIGGP(3, 2)*N(2) + SIGGP(4, 2)*N(3)
         SIGNP(3, 3) =  SIGGP(1, 3)*N(1) + SIGGP(2, 3)*N(4)
     ;                + SIGGP(3, 3)*N(2) + SIGGP(4, 3)*N(3)
C         CALL SFQ4(-SQRT3, SQRT3, N, DNK, DNE)
         CALL SFQ4(-DSQRT(3.D0), DSQRT(3.D0), N, DNK, DNE)
         SIGNP(4, 1) =  SIGGP(1, 1)*N(1) + SIGGP(2, 1)*N(4)
     ;                + SIGGP(3, 1)*N(2) + SIGGP(4, 1)*N(3)
         SIGNP(4, 2) =  SIGGP(1, 2)*N(1) + SIGGP(2, 2)*N(4)
     ;                + SIGGP(3, 2)*N(2) + SIGGP(4, 2)*N(3)
         SIGNP(4, 3) =  SIGGP(1, 3)*N(1) + SIGGP(2, 3)*N(4)
     ;                + SIGGP(3, 3)*N(2) + SIGGP(4, 3)*N(3)
          WRITE(98, '(//A/)') 'Matrix: SIGNP'
         DO 913 K = 1, 4
            WRITE(98, *)
            DO 914 L = 1, 3
               WRITE(98, '(F16.4)', ADVANCE='NO') SIGNP(K, L)
  914       CONTINUE
  913    CONTINUE
         WRITE(98, '(/A, F16.6)'), 'lement Area=', A
         WRITE(98, '(//)')

C
C        Feed these stresses to the nodal stress matrix
C        
         DO 130 K = 1, NECO 
            CNODE = NODE(K)
            DO 140 L = 1, 3
               SIG(CNODE, L) =  SIG(CNODE, L) + SIGNP(K, L)*A
  140       CONTINUE
            SIG(CNODE, 4) = SIG(CNODE, 4) + A
  130    CONTINUE
   60 CONTINUE
      WRITE(98, '(/A/)') 'AT THE END OF THE ELEMENT LOOP'
      WRITE(98, '(A/)') 'Matrix: SIG'
      DO 915 K = 1, NNODE
         WRITE(98, *)
         DO 916 L = 1, 4
            WRITE(98, '(F16.10)', ADVANCE='NO') SIG(K, L)
  916    CONTINUE
  915 CONTINUE
C
C     Loop through nodes, compute derived stresses and write to file
C        
      DO 150 I = 1, NNODE
         SGX   = SIG(I, 1)/ SIG(I, 4)
         SGY   = SIG(I, 2)/ SIG(I, 4)
         SGXY  = SIG(I, 3)/ SIG(I, 4)
         SGMAX = (SGX + SGY)/2.D0
     ;          + DSQRT( ((SGX - SGY)/2.D0)**2 + SGXY**2 )
         SGMIN = (SGX + SGY)/2.D0
     ;          - DSQRT( ((SGX - SGY)/2.D0)**2 + SGXY**2 )
         SGVM  = DSQRT( SGMAX**2 - SGMAX*SGMIN + SGMIN**2 )
               
         WRITE(48, '(E26.16)') SGX
         WRITE(50, '(E26.16)') SGY
         WRITE(52, '(E26.16)') SGXY
         WRITE(54, '(E26.16)') SGVM
         WRITE(56, '(E26.16)') SGMAX
         WRITE(58, '(E26.16)') SGMIN
  150 CONTINUE
      RETURN
C
C     .. End of SQ4 ..
C        
      END
         
