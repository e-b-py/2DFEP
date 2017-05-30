C
C =====================================================================
      SUBROUTINE GEOMET(CRD, LDCRD, CON, LDCON, LBEL, LDLBEL, RBEL,
     ;                  LDRBEL, LNOD, LDLNOD, RNOD, LDRNOD, CMAT, ERST, 
     ;                  LDERST, PRST, LDPRST)
C
C     Read the node, element, material, restraint, and load data
C     data provided in the input and mesh files into the related 
C     arrays.
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCRD, LDCON, LDLBEL, LDRBEL, LDLNOD, LDRNOD,
C                      LDERST, LDPRST
C     ..
C     .. Array Arguments ..
C     INTEGER*4        CON(LDCON, *)  : Element connectivity matrix
C                      ERST(LDERST, *): Matrix of element restraints
C                      PRST(LDPRST, *): Matrix of point restraints
C                      LBEL(LDLBEL, *): Loaded boundary elements
C                      RBEL(LDRBEL, *): Restrained boundary elements
C                      LNOD(LDLNOD, *): Loaded points
C                      RNOD(LDRNOD, *): Restrained points
C
C     REAL*8           CRD(LDCRD, *)  : Matrix of nodal coordinates
C                      CMAT(4)        : Vector of material properties
C     ..
C     .. Local Scalars ..
C     INTEGER*4        NRL  : Number of restraint group labels (b. el.)
C                      NLL  : Number of load group labels (b. el.)
C                      NRPL : Number of restraint group labels (points)
C                      NLPL : Number of load group labels (points)
C                      NT   : Generic, number of something
C                      LABEL: Generic label
C                      LLB  : Line length of boundary elements
C                      LLD  : Line length of domain elements
C                      LLP  : Line length of point elements
C                      TBE  : Tag of boundary elements
C                      TDE  : Tag of domain elements
C                      TPE  : Tag of point elements
C     .. 
C     .. Local Arrays  ..
C     INTEGER*4        GRB(20)        : Vector to read unnecessary data 
C                      LGRP(20)       : Loaded boundary group labels
C                      RGRP(20)       : Restrained boundary group labels
C                      LPGRP(20)      : Loaded point group labels
C                      RPGRP(20)      : Restrained point group labels
C
C     CHARACTER        MESH*80        : Mesh file's name
C                      ESTR*80        : Element type string
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
C                      NBCO : Number of coord. for boundary elements
C                      NECO : Number of coord. for domain elements
C     ,,
C     .. Scalar Arguments ..
      INTEGER*4        LDCRD, LDCON, LDLBEL, LDRBEL, LDLNOD, LDRNOD, 
     ;                 LDERST, LDPRST
C     ..
C     .. Array Arguments ..
      INTEGER*4        CON(LDCON, *), ERST(LDERST, *), PRST(LDPRST, *),
     ;                 RBEL(LDRBEL, *), LBEL(LDLBEL, *),RNOD(LDRNOD, *),
     ;                 LNOD(LDLNOD, *)
      REAL*8           CRD(LDCRD, *), CMAT(4)
C     ..
C =====================================================================
C     .. Local Scalars
      INTEGER*4        NRL, NLL, NRPL, NLPL, NT, LABEL, LLB, LLD, LLP,
     ;                 TBE, TDE, TPE
C     .. 
C     .. Local Arrays
      INTEGER*4        GRB(20), LGRP(20), RGRP(20), LPGRP(20), RPGRP(20)
      CHARACTER        MESH*80, ESTR*80
C     ..
C     .. Common Scalars ..
      INTEGER*4        NNODE, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                 NPRST, NDOF, NNDL, NSL, NVL, ETYPE, PTYPE, NBCO,
     ;                 NECO
      COMMON /CONFIG/  NNODE, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                 NPRST, NDOF, NNDL, NSL, NVL, ETYPE, PTYPE, NBCO,
     ;                 NECO
C     ..
C
C     .. Executable statements ..
C
C     Open the files for I/O
C
C      OPEN(UNIT=11, FILE='input.dat')
      OPEN(UNIT=11, FILE='/dev/stdin')
      OPEN(UNIT=12, FILE='../out/output.dat')
C
C     Read the name of the meshfile, element type and problem type
C
      READ(11, '(//)')
      READ(11, '(A)') MESH
      READ(11, *)
      READ(11, *) ETYPE
      READ(11, *)
      READ(11, *) PTYPE
C
C     Read the material properties and thickness of elements
C
      READ(11, *)
      READ(11, *) CMAT(1), CMAT(2), CMAT(3), CMAT(4)
C
C     Read the restrained and loaded boundaries (groups)
C
      RGRP = 0
      LGRP = 0
      RPGRP = 0
      LPGRP = 0
      READ(11, *)
      READ(11, *) NRL
      DO 10 I = 1, NRL
         READ(11, *) RGRP(I)
   10 CONTINUE
      READ(11, *)
      READ(11, *) NRPL
      DO 20 I = 1, NRPL
         READ(11, *) RPGRP(I)
   20 CONTINUE
      READ(11, *)
      READ(11, *) NLL
      DO 30 I = 1, NLL
         READ(11, *) LGRP(I)
   30 CONTINUE
      READ(11, *)
      READ(11, *) NLPL
      DO 40 I = 1, NLPL
         READ(11, *) LPGRP(I)
   40 CONTINUE
C
C     Open the mesh file, read the nodal coordinates
C
      OPEN(UNIT=13, FILE=MESH)
      READ(13, '(///)')
      READ(13, *), NNODE
      DO 50 I = 1, NNODE
         READ(13, *) GRB(1), CRD(I, 1), CRD(I, 2)
   50 CONTINUE
C
C     READ THE ELEMENTS
C
C     Take the element type and adjust parameters
C
      LLP = 6
      TPE = 15
      IF (ETYPE.EQ.0) THEN
         ESTR = 'CST (T3)'
         LLB  = 7
         LLD  = 8
         TBE  = 1
         TDE  = 2
         NBCO = 2
         NECO = 3
      ELSEIF (ETYPE.EQ.1) THEN
         ESTR = 'LST (T6)'
         LLB  = 8
         LLD  = 11
         TBE  = 8
         TDE  = 9
         NBCO = 3
         NECO = 6
      ELSEIF (ETYPE.EQ.2) THEN
         ESTR = 'ISO-Q4'
         LLB  = 7
         LLD  = 9
         TBE  = 1
         TDE  = 3
         NBCO = 2
         NECO = 4
      ENDIF
C
      GRB = 0
      NELE  = 0
      NRBEL = 0
      NLBEL = 0
      NRNOD = 0
      NLNOD = 0
      K = 1 
      READ(13, '(/)')
      READ(13, *), NT
C     Read points 
      DO 60 I = 1, NT
         READ(13, *) (GRB(J), J=1,LLP)
         IF ( GRB(2).NE.TPE ) GO TO 900
         K = K + 1
         LABEL = GRB(4)
         DO 70 J = 1, NRPL
            IF ( RPGRP(J).EQ.LABEL ) THEN
               NRNOD = NRNOD + 1
               RNOD(NRNOD, 1) = LABEL
               RNOD(NRNOD, 2) = GRB(6)
            ENDIF
   70    CONTINUE
         DO 80 J = 1, NLPL
            IF ( LPGRP(J).EQ.LABEL) THEN
               NLNOD = NLNOD + 1
               LNOD(NLNOD, 1) = LABEL
               LNOD(NLNOD, 2) = GRB(6)
            ENDIF
   80    CONTINUE
   60 CONTINUE
  900 BACKSPACE(13)
C     Read elements
      DO 90 I = K, NT
         READ(13, *) (GRB(J), J=1,LLB)
         IF ( GRB(2).EQ.TBE ) THEN
            LABEL = GRB(4)
            DO 100 J = 1, NRL
               IF ( RGRP(J).EQ.LABEL ) THEN
                  NRBEL = NRBEL + 1
                  RBEL(NRBEL, 1) = LABEL
                  DO 110 K = 1, NBCO
                     RBEL(NRBEL, K+1) = GRB(5+K)
  110             CONTINUE
               ENDIF
  100       CONTINUE
            DO 120 J = 1, NLL
               IF ( LGRP(J).EQ.LABEL ) THEN
                  NLBEL = NLBEL + 1
                  LBEL(NLBEL, 1) = LABEL
                  DO 130 K = 1, NBCO
                     LBEL(NLBEL, K+1) = GRB(5+K)
  130             CONTINUE
               ENDIF
  120       CONTINUE
         ELSEIF (GRB(2).EQ.TDE ) THEN
            BACKSPACE(13)
            READ(13, *) (GRB(J), J=1,LLD)
            NELE = NELE + 1
            DO 140 K = 1, NECO
               CON(NELE, K) = GRB(5+K)
  140       CONTINUE
         ENDIF
   90 CONTINUE
C
C     Read the restraints
C
      READ(11, *)
      READ(11, *) NERST
      DO 150 I = 1, NERST
         READ(11, *) ERST(I, 1), ERST(I, 2)
  150 CONTINUE
      READ(11, *)
      READ(11, *) NPRST
      DO 160 I = 1, NPRST
         READ(11, *) PRST(I, 1), PRST(I, 2)
  160 CONTINUE

C =====================================================================
C 
C Write all to output file
C 
      WRITE(12, '(A)') 'INPUT DATA: [kN/m/C]'
      WRITE(12, '(A)') '-----------------------------------------------'
      WRITE(12, '(A)') '# ELEMENT TYPE'
      WRITE(12, '(A)') ESTR
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# NODE COORDINATES'
      WRITE(12, '(A, T14, A, T26, A)') 'Node', 'x-coord', 'y-coord'
      DO 170 I = 1, NNODE
         WRITE(12, '(I3, T8, 2F12.2)') I, CRD(I, 1), CRD(I, 2)
  170 CONTINUE
      WRITE(12, '(/A)') '# MATERIAL PROPERTIES'
      WRITE(12, '(T5, A, T22, A, T36, A, T52, A)') 
     :           'E', 'v', 'gamma', 't'
      WRITE(12, '(D9.4, T18, F6.2, T34, F6.2, T49, F6.2)')
     ;           (CMAT(I), I=1,4)
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# RESTRAINED POINTS'
         WRITE(12, '(A, T16, A)')
     ;              'Label', 'Node'
         DO 180 I = 1, NRNOD
            WRITE(12, '(I3, T14, I3)')
     ;                  RNOD(I, 1), RNOD(I, 2)
  180    CONTINUE
C ---------------------------------------------------------------------

      WRITE(12, '(/A)') '# RESTRAINED BOUNDARY ELEMENTS'
      IF ( ETYPE.EQ.1 ) THEN
         WRITE(12, '(A, T16, A, T30, A, T44, A)')
     ;              'Label', 'N1', 'N2', 'N3'
         DO 190 I = 1, NRBEL
            WRITE(12, '(I3, T12, I3, T28, I3, T44, I3)')
     ;                  RBEL(I, 1), RBEL(I, 2), RBEL(I, 3), RBEL(I, 4)
  190    CONTINUE
      ELSE
         WRITE(12, '(A, T16, A, T30, A)')
     ;              'Label', 'N1', 'N2'
         DO 200 I = 1, NRBEL
            WRITE(12, '(I3, T14, I3, T28, I3)')
     ;                  RBEL(I, 1), RBEL(I, 2), RBEL(I, 3)
  200    CONTINUE
      ENDIF
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# LOADED POINTS'
      WRITE(12, '(A, T16, A)')
     ;              'Label', 'Node'
      DO 210 I = 1, NLNOD
         WRITE(12, '(I3, T14, I3)')
     ;               LNOD(I, 1), LNOD(I, 2)
  210    CONTINUE
C ---------------------------------------------------------------------

      WRITE(12, '(/A)') '# LOADED BOUNDARY ELEMENTS'
      IF ( ETYPE.EQ.1 ) THEN
         WRITE(12, '(A, T16, A, T30, A, T44, A)')
     ;              'Label', 'N1', 'N2', 'N3'
         DO 220 I = 1, NLBEL
            WRITE(12, '(I3, T12, I3, T28, I3, T44, I3)')
     ;                  LBEL(I, 1), LBEL(I, 2), LBEL(I, 3), LBEL(I, 4)
  220    CONTINUE
      ELSE
         WRITE(12, '(A, T16, A, T30, A)')
     ;              'Label', 'N1', 'N2'
         DO 230 I = 1, NLBEL
            WRITE(12, '(I3, T14, I3, T28, I3)')
     ;                  LBEL(I, 1), LBEL(I, 2), LBEL(I, 3)
  230    CONTINUE
      ENDIF
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# ELEMENTS OF THE BODY'
      IF ( ETYPE.EQ.0 ) THEN
         WRITE(12, '(A, T16, A, T24, A, T32, A)')
     ;              'EL.NO', 'N1', 'N2', 'N3'
         DO 240 I = 1, NELE
            WRITE(12, '(I3, T14, I3, T22, I3, T30, I3)')
     ;                 I, CON(I, 1), CON(I, 2), CON(I, 3)
  240    CONTINUE
      ELSEIF ( ETYPE.EQ.1 ) THEN
         WRITE(12, '(A,T16,A,T24,A,T32,A,T40,A,T48,A,T56,A)')
     ;             'EL.NO', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6'
         DO 250 I = 1, NELE
            WRITE(12, '(I3,T14,I3,T22,I3,T30,I3,T38,I3,T46,I3,T54,I3)')
     ;                 I, (CON(I, K), K=1,6)
  250    CONTINUE
      ELSEIF ( ETYPE.EQ.2 ) THEN
            WRITE(12, '(A, T16, A, T24, A, T32, A, T40, A)')
     ;                'EL.NO', 'N1', 'N2', 'N3', 'N4'
         DO 260 I = 1, NELE
            WRITE(12, '(I3, T14, I3, T22, I3, T30, I3, T38, I3)')
     ;                 I, (CON(I, K), K=1,4)
  260    CONTINUE
      ENDIF
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# RESTRAINTS (POINTS)'
      WRITE(12, '(A, T20, A)')
     ;           'Group', 'Direction'
      DO 270 I = 1, NPRST
         WRITE(12, '(I3, T20, I3)')
     ;              PRST(I, 1), PRST(I, 2)
  270 CONTINUE
C ---------------------------------------------------------------------
      WRITE(12, '(/A)') '# RESTRAINTS (BOUNDARY ELEMENTS)'
      WRITE(12, '(A, T20, A)')
     ;           'Group', 'Direction'
      DO 280 I = 1, NERST
         WRITE(12, '(I3, T20, I3)')
     ;              ERST(I, 1), ERST(I, 2)
  280 CONTINUE
      RETURN
C
C     .. End of GEOMET ..
C
      END

