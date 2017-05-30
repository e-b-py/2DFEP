
C =====================================================================
      SUBROUTINE POST(CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, 
     ;                ELDS, LDELDS, DSPG, D, DOUT)
C
C     Output the requested nodal displacements, compute the gradients of
C     the solution
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCON, LDCRD, LDIDOF, LDELDS
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Connectivity matrix
C                      IDOF(LDIDOF, *): Matrix of the degrees of freedom
C                      DOUT(*)        : Displacement output nodes
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of coordinates
C                      ELDS(LDELDS, *): Matrix of element loads
C                      DSPG(*)        : Global nodal displacement vector
C                      D(3, 3)        : Material stiffness matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        VTKCT    : Cell type for vtk file
C                      DSCALE   : Deformation scale
C                      NODN     : Number of the node 
C                      CDOF     : Current DOF in generic loop
C
C     REAL*8           UX       : Displacement along x
C                      UY       : Displacement along y
C
C     ..
C     .. Local Arrays ..
C     INTEGER*4        ELDOF(6) : Global DOFs corresponding to locals 
C
C     REAL*8    
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
      INTEGER*4        LDCON, LDIDOF, LDCRD, LDELDS
C     ..
C     .. Array Arguments ..
      INTEGER*4        CON(LDCON, *), IDOF(LDIDOF, *), DOUT(*)
      REAL*8           CRD(LDCRD, *), ELDS(LDELDS, *), DSPG(*), D(3, 3)
C     ..
C =====================================================================
C     .. Local Scalars ..
      INTEGER*4        VTKCT, NODN, CDOF
      REAL*8           DSCALE, UX, UY
C     ..
C     .. Local Arrays ..

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
C     Get cell type for Vtk
C
      IF ( ETYPE.EQ.0 ) VTKCT = 5
      IF ( ETYPE.EQ.1 ) VTKCT = 22
      IF ( ETYPE.EQ.2 ) VTKCT = 9
C
C     Write the undeformed mesh to vtk file
C
      WRITE(44, '(A)') '# vtk DataFile Version 2.0'
      WRITE(44, '(A)') 'undeformed mesh'
      WRITE(44, '(A)') 'ASCII'
      WRITE(44, '(A)') 'DATASET UNSTRUCTURED_GRID'
      WRITE(44, '(A, I6, A)') 'POINTS', NNODE, ' float'
      DO 10 I = 1, NNODE
         WRITE(44, '(F22.16, F22.16, F22.16)') 
     ;               CRD(I, 1), CRD(I, 2), 0.D0
   10 CONTINUE
      WRITE(44, '(/A, I6, I6)') 'CELLS', NELE, NELE*(NECO+1)
      DO 20 I = 1, NELE
         WRITE(44, '(I6)', ADVANCE='NO') NECO
         DO 30 J = 1, NECO
            WRITE(44, '(I6)', ADVANCE='NO') CON(I, J)-1
   30    CONTINUE
         WRITE(44, *)
   20 CONTINUE
      WRITE(44, '(/A, I6)') 'CELL_TYPES', NELE
      DO 40 I = 1, NELE
         WRITE(44, '(I6)') VTKCT
   40 CONTINUE
C
C     Update coordinates
C
      DSCALE = 500.D0
      DO 50 I = 1, NNODE
         DO 60 J = 1, 2
            IF ( IDOF(I, J).EQ.-1 ) GO TO 50
            CRD(I, J) = CRD(I, J) + DSPG(IDOF(I, J))*DSCALE
   60    CONTINUE
   50 CONTINUE
C
C     Write deformed mesh to vtk file
C
      WRITE(46, '(A)') '# vtk DataFile Version 2.0'
      WRITE(46, '(A)') 'undeformed mesh'
      WRITE(46, '(A)') 'ASCII'
      WRITE(46, '(A)') 'DATASET UNSTRUCTURED_GRID'
      WRITE(46, '(A, I6, A)') 'POINTS', NNODE, ' float'
      DO 70 I = 1, NNODE
         WRITE(46, '(F22.16, F22.16, F22.16)') 
     ;               CRD(I, 1), CRD(I, 2), 0.D0
   70 CONTINUE
      WRITE(46, '(/A, I6, I6)') 'CELLS', NELE, NELE*(NECO+1)
      DO 80 I = 1, NELE
         WRITE(46, '(I6)', ADVANCE='NO') NECO
         DO 90 J = 1, NECO
            WRITE(46, '(I6)', ADVANCE='NO') CON(I, J)-1
   90    CONTINUE
         WRITE(46, *)
   80 CONTINUE
      WRITE(46, '(/A, I6)') 'CELL_TYPES', NELE
      DO 100 I = 1, NELE
         WRITE(46, '(I6)') VTKCT
  100 CONTINUE
C
C     Write the displacements of the specified output nodes to file
C
      WRITE(42, '(A, T15, A, T40, A)') 'Node', 'u(mm)', 'v(mm)'
      DO 110 I = 1, 50
         IF ( DOUT(I).EQ.0 ) GO TO 799
         WRITE(42, '(I4, T10, F15.12, T35, F15.12)')
     ;   DOUT(I), DSPG(IDOF(DOUT(I), 1))*1000.D0,
     ;            DSPG(IDOF(DOUT(I), 2))*1000.D0
  110 CONTINUE
C     Compute and write the stresses
C
  799 IF ( ETYPE.EQ.0 ) GO TO 800
      IF ( ETYPE.EQ.2 ) GO TO 802
  800 CALL ST3(CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, DSPG, D)
  802 RETURN
      END
         
