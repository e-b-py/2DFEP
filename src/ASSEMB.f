C
C =====================================================================
      SUBROUTINE ASSEMB(CON, LDCON, IDOF, LDIDOF, CRD, LDCRD, LBEL,
     ;                  LDLBEL, LNOD, LDLNOD, CMAT, SLDS, LDSLDS, NLDS,
     ;                  VLDS, LDVLDS, KG, LDKG, FG, D)
C
C     Assemble the global stiffness matrix and force vector.
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCON, LDCRD, LDIDOF, LDLBEL, LDLNOD, LDSLDS,
C                      LDVLDS, LDKG
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Connectivity matrix
C                      IDOF(LDIDOF, *): Matrix of the degrees of freedom
C                      LBEL(LDLBEL, *): Loaded boundary elements
C                      LNOD(LDLNOD, *): Loaded points
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of coordinates
C                      CMAT(4)        : Vector of material properties
C                      SLDS(LDSLDS, *): Matrix of surface loads
C                      NLDS(*)        : Vector of nodal loads
C                      VLDS(LDVLDS, *): Matrix of volume loads
C                      KG(LDKG, *)    : Global stiffness matrix
C                      FG(*)          : Global force matrix
C                      D(3, 3)        : Material stiffness matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        SNOD : Start node of the element
C                      ENOD : End node of the element
C                      RI   : Raw index  
C                      CI   : Column index  
C                      
C     REAL*8           E    : Young's Modulus
C                      NU   : Poisson's Ratio
C                      NELDOF  : Number of element DOFS
C     ..
C     .. Local Arrays ..
C     INTEGER*4        ELDOF(12): Global DOFs of the element 
C                     
C     REAL*8           KL(12,12)      : Local stiffness matrix
C                      FL(12)         : Local force matrix
C       
C     .. Common Scalars ..
C     INTEGER*4        NNODE: number of nodes
C                      NELE : number of elements
C                      NRBEL: number of restrained boundary elements
C                      NLBEL: number of loaded boundary elements
C                      NRNOD: number of restrained points
C                      NLNOD: number of loaded points
C                      NERST: number of element restraints
C                      NPRST: number of element restraints
C                      NDOF : number of degrees of freedom
C                      NNDL : number of nodal loads
C                      NSL  : number of element loads
C                      NVL  : number of volume loads
C                      ETYPE: element type specifier
C                      PTYPE: constitutive relation specifier
C                      NBCO : number of coord. for boundary elements
C                      NECO : number of coord. for domain elements
C     ,,
C     .. Scalar Arguments ..
      INTEGER*4        LDCON, LDIDOF, LDCRD, LDLBEL, LDLNOD, LDSLDS,
     ;                 LDVLDS, LDKG
C     ..
C     .. Array Arguments ..
      INTEGER*4        CON(LDCON, *), IDOF(LDIDOF, *), LBEL(LDLBEL, *),
     ;                 LNOD(LDLNOD, *)
      REAL*8           CRD(LDCRD, *), CMAT(4), SLDS(LDSLDS, *), NLDS(*),
     ;                 VLDS(LDVLDS, *), KG(LDKG, *), FG(*), D(3, 3)
C     ..
C =====================================================================
C     .. Local Scalars ..
      INTEGER*4        SNOD, ENOD, RI, CI, NELDOF
      REAL*8           E, NU
C     ..
C     .. Local Arrays ..
      INTEGER*4        ELDOF(12)
      REAL*8           KL(12, 12), FL(12)
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
CCCCCCCCCCCCCCCCCCCCCC   D E B U G   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(26, '(/A//)') 'XASSEMB DEBUG FILE'
      WRITE(26, '(A/)') 'In  ASSEMB.f'
      WRITE(26, *) 'Matrix: CON'
      WRITE(26, '(3I5)')
      DO 901 K = 1, NELE
         WRITE(26, *) CON(K, 1), CON(K, 2), CON(K, 3)
  901 CONTINUE
      WRITE(26, '(/A/)') 'Matrix: CRD'
      WRITE(26, *)
      DO 902 K = 1, NNODE
         WRITE(26, '(F6.2, F6.2)') CRD(K, 1), CRD(K, 2)
  902 CONTINUE
      WRITE(26, '(/A/)') 'Matrix: IDOF'
      WRITE(26, *)
      DO 903 K = 1, NNODE
         WRITE(26, '(I6, I6)') IDOF(K, 1), IDOF(K, 2)
  903 CONTINUE
CCCCCCCCCCCCCCCCCCCCCC   D E B U G   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Initialize D according to specified PTYPE
C
      E  = CMAT(1)
      NU = CMAT(2)
      IF ( PTYPE.EQ.0 ) THEN
         D(1, 1) = E*(1 - NU) / ( (1 + NU)*(1 - 2*NU) )
         D(2, 2) = D(1, 1)
         D(3, 3) = E / ( 2*(1 + NU) )
         D(1, 2) = E*NU / ( (1 + NU)*(1 - 2*NU) )
         D(2, 1) = D(1, 2)
      ELSEIF ( PTYPE.EQ.1 ) THEN
         D(1, 1) = E / (1 - NU**2)
         D(2, 2) = D(1, 1)
         D(3, 3) = E*(1 - NU) / (2*(1 - NU**2))
         D(1, 2) = D(1, 1)*NU
         D(2, 1) = D(1, 2)
      ENDIF
      WRITE(26, '(/A/)') 'Matrix: D'
      WRITE(26, *)
      DO 904 K = 1, 3
         WRITE(26, '(F20.4, F20.4, F20.4)') D(K, 1), D(K, 2), D(K, 3)
  904 CONTINUE
C
C     Loop over elements
C
      NELDOF = NECO*2
      WRITE(26, '(//A/)') 'LOOP OVER ELEMENTS:'
      DO 10 I = 1, NELE
         ELDOF = 0
         WRITE(26, '(//A)') '------------------------------------------'
         WRITE(26, '(/A, I3/)') 'ELEMENT', I
C
C     Compute element stiffness matrices and integrate volume forces
C
         CALL MKK(I, CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, LBEL, LDLBEL,
     ;            CMAT, SLDS, LDSLDS, VLDS, LDVLDS, D, KL, FL, ELDOF)
         
C
C        Assemble 
C
         WRITE(26, '(/A/)') 'ELDOF vector'
         WRITE(26, *)
         WRITE(26, '(I3)') (ELDOF(K), K=1,NELDOF)
         WRITE(26, '(/A/)') 'Matrix: KL (Local Stiffness Matrix)'
         DO 905 K = 1, NELDOF
            WRITE(26, *)
            DO 906 L = 1, NELDOF
               WRITE(26, '(E10.2)', ADVANCE='NO') KL(K, L)
  906       CONTINUE
  905    CONTINUE
         WRITE(26, '(/A/)') 'Matrix: FL (Local Force Matrix)'
         DO 907 K = 1, NELDOF
            WRITE(26, '(F8.2)') FL(K)
  907    CONTINUE
C
         DO 40 M = 1, NELDOF
            RI = ELDOF(M)
            IF( RI.EQ.-1 ) GO TO 40
            FG(RI) = FG(RI) + FL(M)
            DO 50 N = 1, NELDOF
               CI = ELDOF(N)
               IF ( CI.EQ.-1) GO TO 50
               KG(RI, CI) = KG(RI, CI) + KL(M, N)
   50       CONTINUE
   40    CONTINUE
   10 CONTINUE
C
C     Integrate surface forces and add them to FG
C
      CALL SURFL(IDOF, LDIDOF, CRD, LDCRD, SLDS, LDSLDS, LBEL, LDLBEL,
     ;           FG)
C
C     Add nodal loads to FG <TODO><TODO> DO THIS IN VLOADS <TODO><TODO>
C
      DO 60 I = 1, NDOF
         FG(I) = FG(I) + NLDS(I)
   60 CONTINUE
      RETURN
C
C     .. End of ASSEMB ..
C
      END

         
