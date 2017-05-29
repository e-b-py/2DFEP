C =====================================================================
      SUBROUTINE ST3(CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, 
     ;                DSPG, D)
C
C     Compute sig_xx, sig_yy, sig_xy, principal min and max stresses, 
C     their orientations and Von Mises stresses for T3 elements
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCON, LDCRD, LDIDOF
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Connectivity matrix
C                      IDOF(LDIDOF, *): Matrix of the degrees of freedom
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of coordinates
C                      DSPG(*)        : Global nodal displacement vector
C                      D(3, 3)        : Material stiffness matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        VTKCT    : Cell type for vtk file
C                      DSCALE   : Deformation scale
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
C     ..
C     .. Local Arrays ..
C     INTEGER*4        ELDOF(12)  : Global DOFs corresponding to locals 
C
C     REAL*8           COEFF(3, 3): Coefficients of the element
C                      DSPL(6)    : Local displacements of the element
C                      B(3, 6)    : Strain interpolation matrix
C                      SIG(3)     : Stress vector                        
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
C     .. Local Scalars ..
      INTEGER*4        VTKCT
      REAL*8           A, B1, B2, B3, C1, C2, C3, D11, D12, D21, D22,
     ;                 D33, SGX, SGY, SGXY, SGVM, SGMAX, SGMIN, ORMAX,
     ;                 ORMIN
C     ..
C     .. Local Arrays ..
      INTEGER*4        ELDOF(12)
      REAL*8           COEFF(3, 3), DSPL(6), B(3, 6), SIG(3)
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
      REWIND(92)
      ELDOF = 0
      COEFF = 0.D0
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
     ;                  CRD(J, 1), CRD(J, 2), 0.D0
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
      WRITE(48, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(48, '(A)') 'SCALARS sig_xx float 1'
      WRITE(48, '(A)') 'LOOKUP_TABLE default'
      WRITE(50, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(50, '(A)') 'SCALARS sig_yy float 1'
      WRITE(50, '(A)') 'LOOKUP_TABLE default'
      WRITE(52, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(52, '(A)') 'SCALARS sig_xy float 1'
      WRITE(52, '(A)') 'LOOKUP_TABLE default'
      WRITE(54, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(54, '(A)') 'SCALARS sig_vm float 1'
      WRITE(54, '(A)') 'LOOKUP_TABLE default'
      WRITE(56, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(56, '(A)') 'SCALARS sig_max float 1'
      WRITE(56, '(A)') 'LOOKUP_TABLE default'
      WRITE(58, '(/A, I6)') 'CELL_DATA', NELE
      WRITE(58, '(A)') 'SCALARS sig_min float 1'
      WRITE(58, '(A)') 'LOOKUP_TABLE default'
C
C     <<<  ELEMENT LOOP >>>
C
      DO 60 I = 1, NELE
         READ(92) ELDOF
         READ(92) COEFF
         READ(92) A
         B(1, 1) = COEFF(1, 2)
         B(1, 3) = COEFF(2, 2)
         B(1, 5) = COEFF(3, 2)
         B(2, 2) = COEFF(1, 3)
         B(2, 4) = COEFF(2, 3)
         B(2, 6) = COEFF(3, 3)
         B(3, 1) = B(2, 2)
         B(3, 2) = B(1, 1)
         B(3, 3) = B(2, 4)
         B(3, 4) = B(1, 3)
         B(3, 5) = B(2, 6)
         B(3, 6) = B(1, 5)
C
C        Retrieve local displacements
C
         DSPL = 0.D0
         DO 70 J = 1, 6
            IF ( ELDOF(J).NE.-1 ) THEN
               DSPL(J) = DSPG(ELDOF(J))
            ENDIF
   70    CONTINUE
C
C        Compute Stresses
C        
         SIG = MATMUL( D, MATMUL(B, DSPL) )
         SGX = SIG(1)
         SGY = SIG(2)
         SGXY= SIG(3)
         SGMAX = (SGX + SGY)/2.D0 
     ;          + DSQRT( ((SGX - SGY)/2.D0)**2 + SGXY**2 )
         SGMIN = (SGX + SGY)/2.D0 
     ;          - DSQRT( ((SGX - SGY)/2.D0)**2 + SGXY**2 )
         SGVM = DSQRT( SGMAX**2 - SGMAX*SGMIN + SGMIN**2 )
         WRITE(48, '(E26.16)') SGX
         WRITE(50, '(E26.16)') SGY
         WRITE(52, '(E26.16)') SGXY
         WRITE(54, '(E26.16)') SGVM
         WRITE(56, '(E26.16)') SGMAX
         WRITE(58, '(E26.16)') SGMIN
   60 CONTINUE
      RETURN
      END
         
