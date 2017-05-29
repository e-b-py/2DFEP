C
C =====================================================================
      SUBROUTINE LOADS(IDOF, LDIDOF, NLDS, SLDS, LDSLDS, VLDS, LDVLDS,
     ;                 LNOD, LDLNOD )
C
C     Reads the input data, processes the loads and assigns them to the
C     correct positions in the corresponding arrays.
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDIDOF, LDSLDS, LDVLDS, LDLNOD
C     ..
C     .. Array Arguments ..
C     REAL*8           NLDS(*)        : Vector of nodal loads
C                      SLDS(LDSLDS, *): Matrix of element loads
C                      VLDS(LDVLDS, *): Matrix of volume loads
C              
C     INTEGER*4        IDOF(LDIDOF, *): Matrix of the degrees of freedom
C                      LNOD(LDLNOD, *): Loaded points
C     ..
C     .. Local Scalars ..
C     REAL*8           CMAGN: Magnitude of the nodal load
C                      FX   : Force in x direction
C                      FY   : Force in y direction
C   
C     INTEGER*4        CNODE: Node to which load is applied
C                      CDIR : Direction of the applied nodal load
C                      CDOF : DOF corresponding to the node and
C                             direction of a given nodal load
C                      LABEL: Generic group label
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
C                      PTYPE: constitutive relation specifier
C                      NBCO : number of coord. for boundary elements
C                      NECO : number of coord. for domain elements
C     ,,
C     .. Scalar Arguments ..
      INTEGER*4        LDIDOF, LDSLDS, LDVLDS, LDLNOD
C     ..
C     .. Array Arguments ..
      INTEGER*4        IDOF(LDIDOF, *), LNOD(LDLNOD, *)
      REAL*8           NLDS(*), SLDS(LDSLDS, *), VLDS(LDVLDS, *)
C     ..
C =====================================================================
C     .. Local Scalars ..
      REAL*8           FX, FY
      INTEGER*4        CNODE, CDIR, CDOF, LABEL
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
C     Read the nodal loads, assign them to the nodal load vector.
C
      WRITE(12, '(/A)') '# NODAL LOADS'
      WRITE(12, '(A, T14, A, T32, A)') 'Label', 'Direction', 'Magnitude'
      READ(11, *)
      READ(11, *) NNDL
      DO 10 I = 1, NNDL
         READ(11, *) LABEL, CDIR, CMAGN
         WRITE(12, '(I3, T16, I3, T32, F6.2)') LABEL, CDIR, CMAGN
         DO 20 J = 1, NLNOD
            IF ( LNOD(J, 1).EQ.LABEL ) THEN
               CNODE = LNOD(J, 2)
               CDOF  = IDOF(CNODE, CDIR)
               NLDS(CDOF) = CMAGN
            ENDIF
   20    CONTINUE
   10 CONTINUE
C
C     Read the volume loads, assign them to the volume load matrix.
C
      WRITE(12, '(/A)') '# VOLUME LOADS'
      WRITE(12, '(A, T15, A)')
     ;            'Fx', 'Fy'
      READ(11, *)
      READ(11, *) NVL
      DO 30 I = 1, NVL
         READ(11, *) VLDS(I, 1), VLDS(I, 2)
         WRITE(12, '(F6.2, T15, F6.2)')
     ;               VLDS(I, 1), VLDS(I, 2)
   30 CONTINUE
C
C     Read the surface loads, assign them to the surface load matrix.
C
      WRITE(12, '(/A)') '# SURFACE LOADS'
      WRITE(12, '(A, T14, A, T24, A, T34, A, T44, A, T54, A, T64, A)')
     ;            'Label', 'Fx1', 'Fx2', 'Fx3', 'Fy1', 'Fy2', 'Fy3'
      READ(11, *)
      READ(11, *) NSL
      DO 40 I = 1, NSL
         READ(11, *) (SLDS(I, J), J=1,7)
C         READ(11, *) SLDS(I, 1), SLDS(I, 2), SLDS(I, 3)
         WRITE(12, '(I3, T8, 6F10.2)')
     ;              IDINT(SLDS(I, 1)), (SLDS(I, J), J=2,7)
   40 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DEBUG LOG CCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(24, '(A//)') 'XLOADS DEBUG FILE'
      WRITE(24, '(A/)') 'All common parameters'
      WRITE(24, '(A, T7, A, I3)') 'NNODE', ':', NNODE
      WRITE(24, '(A, T7, A, I3)') 'NELE', ':', NELE
      WRITE(24, '(A, T7, A, I3)') 'NRBEL', ':', NRBEL
      WRITE(24, '(A, T7, A, I3)') 'NLBEL', ':', NLBEL
      WRITE(24, '(A, T7, A, I3)') 'NRNOD', ':', NRNOD
      WRITE(24, '(A, T7, A, I3)') 'NLNOD', ':', NLNOD
      WRITE(24, '(A, T7, A, I3)') 'NERST', ':', NERST
      WRITE(24, '(A, T7, A, I3)') 'NPRST', ':', NPRST
      WRITE(24, '(A, T7, A, I3)') 'NDOF', ':', NDOF
      WRITE(24, '(A, T7, A, I3)') 'NNDL', ':', NNDL
      WRITE(24, '(A, T7, A, I3)') 'NSL', ':', NSL
      WRITE(24, '(A, T7, A, I3)') 'NVL', ':', NVL
      WRITE(24, '(A, T7, A, I3)') 'ETYPE', ':', ETYPE
      WRITE(24, '(A, T7, A, I3)') 'PTYPE', ':', PTYPE
      WRITE(24, '(A, T7, A, I3)') 'NBCO', ':', NBCO
      WRITE(24, '(A, T7, A, I3)') 'NECO', ':', NECO
      WRITE(24, '(//A/)') 'IDOF MATRIX'
      DO 901 K = 1, NNODE
         WRITE(24, '(I5, I5)') IDOF(K, 1), IDOF(K, 2)
  901 CONTINUE
      RETURN
C      
C     .. End of LOADS ..      
C      
      END
