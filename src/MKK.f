C
C =====================================================================
      SUBROUTINE MKK(I, CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, LBEL,
     ;               LDLBEL, CMAT, SLDS, LDSLDS, VLDS, LDVLDS, D, KL,
     ;               FL, ELDOF)
C
C     Computes the local stiffness and force matrix of the given element
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCON, LDCRD, LDIDOF, LDLBEL, LDSLDS, LDVLDS, I
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Connectivity matrix
C                      IDOF(LDIDOF, *): Matrix of the degrees of freedom
C                      LBEL(LDLBEL, *): Loaded boundary elements
C                      ELDOF(12)      : Global DOFs of the element
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of coordinates
C                      CMAT(4)        : Vector of material properties
C                      SLDS(LDSLDS, *): Matrix of surface loads
C                      VLDS(LDVLDS, *): Matrix of volume loads
C                      D(3, 3)        : Material stiffness matrix
C                      KL(12, 12)     : Local stiffness matrix
C                      FL(12)         : Local force matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        
C                      
C     REAL*8   
C     ..
C     .. Local Arrays ..
C     INTEGER*4        NODE(6)        : Nodes of the element
C     REAL*8           ELCOR(6, 2)    : Nodal coordinates of the element
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
      INTEGER*4        I, LDCON, LDCRD, LDIDOF, LDLBEL, LDSLDS, LDVLDS
C     ..
C     .. Array Arguments ..
      INTEGER*4        CON(LDCON, *), IDOF(LDIDOF, *), LBEL(LDLBEL, *),
     ;                 ELDOF(12)
      REAL*8           CRD(LDCRD, *), CMAT(4), SLDS(LDSLDS, *),
     ;                 VLDS(LDVLDS, *), D(3,3),  KL(12, 12), FL(12)
C     ..
C =====================================================================
C     .. Local Scalars ..
      INTEGER*4        GTL
C      REAL*8  
C     ..
C     .. Local Arrays ..
      INTEGER*4        NODE(6)
      REAL*8           ELCOR(6, 2)
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
C
C     Associate global degrees of freedom with the locals
C     
      NODE = 0
      DO 10 J = 1, NECO
         NODE(J) = CON(I, J)
   10 CONTINUE
      K = 1 
      DO 20 J = 1, NECO
         ELDOF(K)   = IDOF(NODE(J), 1)
         ELDOF(K+1) = IDOF(NODE(J), 2)
         K = K + 2
   20 CONTINUE
      WRITE(92) ELDOF
C
C     Retrieve nodal coordinates 
C
      ELCOR = 0.D0
      DO 30 J = 1, NECO
         ELCOR(J, 1) = CRD(NODE(J), 1)
         ELCOR(J, 2) = CRD(NODE(J), 2)
   30 CONTINUE
C
C     Go to the related element's routine
C
      WRITE(30, '(A, I5//)') 'Elemenent:', I
      IF ( ETYPE.EQ.0 ) GO TO 800
C      IF ( ETYPE.EQ.1 ) GO TO 801
C      IF ( ETYPE.EQ.2 ) GO TO 802
  800 CALL T3(ELCOR, CMAT, D, KL, VLDS, LDVLDS, FL)
C
C    1 
C
C    2 
  107 RETURN
C
C     .. End of MKK ..
C
      END
