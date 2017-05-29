C
C =====================================================================
      SUBROUTINE SCODE(IDOF, LDIDOF, RBEL, LDRBEL, RNOD, LDRNOD, ERST,
     ;                 LDERST, PRST, LDPRST)
C
C     Numerates the degrees of freedom considering the restraints
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDIDOF, LDRBEL, LDRNOD, LDERST, LDPRST
C     ..
C     .. Array Arguments ..
C     INTEGER*4        IDOF(LDIDOF, *): Matrix of the degrees of freedom
C                      RBEL(LDRBEL, *): Matrix of restrained b. elements
C                      RNOD(LDRNOD, *): Matrix of restrained points
C                      ERST(LDERST, *): Matrix of element restraints
C                      PRST(LDPRST, *): Matrix of point restraints
C     ..
C     .. Local Scalars ..
C     INTEGER*4        LABEL: Generic group label
C                      
C     ..
C     .. Local Arrays ..
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
C                      NNDL : number of nodel loads
C                      NSL  : number of surface loads
C                      NVL  : number of volume loads
C                      ETYPE: element type specifier
C                      PTYPE: constitutive relation specifier
C                      NBCO : number of coord. for boundary elements
C                      NECO : number of coord. for domain elements
C     ,,
C     .. Scalar Arguments ..
      INTEGER*4        LDIDOF, LDRBEL, LDRNOD, LDERST, LDPRST
C     ..
C     .. Array Arguments ..
      INTEGER*4        IDOF(LDIDOF, *),RBEL(LDRBEL, *),RNOD(LDRNOD, *),
     ;                 ERST(LDERST, *), PRST(LDPRST, *)
C     ..
C =====================================================================
C     .. Local Scalars ..
      INTEGER*4        LABEL
C     ,,
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
C     Loop through the element restraints, assign -1 to corresponding DOFs
C
      DO 10 I = 1, NERST
         LABEL = ERST(I, 1)
         DO 20 J = 1, NRBEL
            IF ( RBEL(J, 1).EQ.LABEL ) THEN
               DO 30 K = 1, NBCO
                  IDOF( RBEL(J, K+1), ERST(I, 2) ) = -1
   30          CONTINUE
            ENDIF
   20    CONTINUE
   10 CONTINUE
C
C     Loop through the point restraints, assign -1 to corresponding DOFs
C
      DO 40 I = 1, NPRST
         LABEL = PRST(I, 1)
         DO 50 J = 1, NRNOD
            IF ( RNOD(J, 1).EQ.LABEL ) THEN
               IDOF( RNOD(J, 2), PRST(I, 2) ) = -1
            ENDIF
   50    CONTINUE              
   40 CONTINUE
C
C     Numerate the DOFS
C
      NDOF = 0
      DO 60 I = 1, NNODE
         DO 70 J = 1, 2
            IF( IDOF(I, J).EQ.-1 ) GO TO 70
            NDOF = NDOF + 1
            IDOF(I, J) = NDOF
   70    CONTINUE
   60 CONTINUE
      RETURN
C
C     .. End of SCODE ..
C
      END
