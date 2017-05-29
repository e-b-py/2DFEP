C
C =====================================================================
      SUBROUTINE SURFL(IDOF, LDIDOF, CRD, LDCRD, SLDS, LDSLDS, LBEL,
     ;                 LDLBEL, FG)
C
C     Integrate surface forces and add them to the global force matrix
C
C     .. Scalar Arguments ..
C     INTEGER*4        LDCRD, LDSLDS, LDLBEL
C     ..
C     .. Array Arguments ..
C     INTEGER*4        LBEL(LDLBEL, *): Loaded boundary elements
C
C     REAL*8           IDOF(LDIDOF, *): Matrix of degrees of freedom
C                      CRD(LDCRD, *)  : Matrix of coordinates
C                      SLDS(LDSLDS, *): Matrix of surface loads
C                      FG(*)          : Global force matrix
C     ..
C     .. Local Scalars ..
C     INTEGER*4        STATE          : Used to designate label change
C                     
C     REAL*8           BELCR(NLBEL*3) : Coordinates of the boundary els.
C                  
C     ..
C     .. Local Arrays ..
C     INTEGER*4        
C     
C     REAL*8          
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
      INTEGER*4        LDCRD, LDSLDS, LDLBEL
C     ..
C     .. Array Arguments ..
      INTEGER*4        IDOF(LDIDOF, *), LBEL(LDLBEL, *) 
      REAL*8           CRD(LDCRD, *), SLDS(LDSLDS, *), FG(*) 
     ;                
C     ..
C =====================================================================
C     .. Local Scalars ..
      INTEGER*4        IND, AIND, ECOUNT, RI
      REAL*8           FX1,FX2, FX3, FY1, FY2, FY3, XDIFF, YDIFF, CDIFF,
     ;                 FXDIFF, FYDIFF, XSLOP, YSLOP, PX1, PX2, PY1, PY2
     ;                 ADIFF, LENGTH
C     ..
C     .. Local Arrays ..
      INTEGER*4        NODE(NLBEL*3)
      REAL*8           VCR(NLBEL*6)
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
C      IF ( ETYPE.EQ. 2 ) GO TO 98
C
C     CST and ISOQ4: 2 node line elements at the boundary
C
      WRITE(32, '(A//)') 'XSURFL DEBUG FILE'
      WRITE(32, '(A, I3)') 'NSL  :', NSL
      WRITE(32, '(A, I3/)') 'NLBEL:', NLBEL
      DO 10 J = 1, NSL
         WRITE(32, '(A, I3, A, I3/)') 'SURFL MAIN LOOP, J =', J,' -',NSL
         WRITE(32, '(A, T15, A, T30, A, T45, A, T60, A)') 
     ;         'LABEL', 'FX1', 'FX2', 'FY1', 'FY2'
         LABEL = SLDS(J, 1)
         FX1 = SLDS(J, 2)
         FX2 = SLDS(J, 3)
         FY1 = SLDS(J, 5)
         FY2 = SLDS(J, 6)
         WRITE(32, '(I3, T13, F6.2, T28, F6.2, T43, F6.2, T58, F6.2/)') 
     ;         LABEL, FX1, FX2, FY1, FY2
         VCR = 0.D0
         M = 1
         N = 1
         ECOUNT = 0
         WRITE(32, '(A, I3/)') 
     ;       'Inner Loop to collect b.el. of label SURFL(J, 1) =', LABEL
         DO 20 K = 1, NLBEL
            IF ( LBEL(K, 1).EQ.LABEL ) THEN
               NODE(N) = LBEL(K, 2)
               NODE(N+1) = LBEL(K, 3)
               VCR(M)   = CRD(LBEL(K, 2), 1)
               VCR(M+1) = CRD(LBEL(K, 2), 2)
               VCR(M+2) = CRD(LBEL(K, 3), 1)
               VCR(M+3) = CRD(LBEL(K, 3), 2)
               M = M + 4
               N = N + 2
               ECOUNT = ECOUNT + 1
            ELSEIF ( (LBEL(K, 1).NE.LABEL).AND.(M.NE.1) ) THEN
               GO TO 21
            ENDIF
C     Break the loop if it once entered the above if but now fails
C     (i.e. now k is at the next boundary element group)
   20    CONTINUE
         WRITE(32, '(A/)') 'After the inner loop:'
         WRITE(32, '(A, I3/)') 'ECOUNT=', ECOUNT
         WRITE(32, '(A/)') 'NODE vector'
         DO 991 K = 1, ECOUNT*2
            WRITE(32, '(I3)') NODE(K)
  991    CONTINUE
         WRITE(32, '(/A/)') 'VCR vector (N1x, N1y, N2x, N2y ...)'
         DO 992 K = 1, ECOUNT*4
            WRITE(32, '(F6.2)') VCR(K)
  992    CONTINUE
   21    YDIFF = VCR(ECOUNT*4) - VCR(2)
         XDIFF = VCR(ECOUNT*4-1) - VCR(1)
         LXDIFF= FX2 - FX1
         LYDIFF= FY2 - FY1
         IF ( DABS(XDIFF).GE.DABS(YDIFF) ) THEN
            XSLOP = LXDIFF / XDIFF
            YSLOP = LYDIFF / XDIFF
            IND = 0
            AIND = 1
         ELSE
            XSLOP = LXDIFF / YDIFF
            YSLOP = LYDIFF / YDIFF
            IND = 1
            AIND = 0
         ENDIF
         CFX1 = FX1
         WRITE(32, '(/A/)') 'Load eval.loop for collected the elements'
         WRITE(32, '(A, F6.2)') 'CFX1=',  CFX1
         CFY1 = FY1
         WRITE(32, '(A, F6.2/)') 'CFY1=',  CFY1
         M = 1
         N = 1
         DO 30 K = 1, ECOUNT
            WRITE(32, '(A, I3, A, I3/)') 'K =',  K, ' -', ECOUNT
            CDIFF = DABS(VCR(M+IND+2) - VCR(M+IND))
            WRITE(32, '(A, T10, A,  F6.2)') 'CDIFF', '=',  CDIFF
            ADIFF = VCR(M+AIND+2) - VCR(M+AIND)
            WRITE(32, '(A, T10, A, F6.2)') 'ADIFF', '=',  ADIFF
            LENGTH= DSQRT(CDIFF**2 + ADIFF**2)
            WRITE(32, '(A, T10, A, F6.2)') 'LENGTH', '=', LENGTH
            CFX2  = CFX1 + XSLOP*CDIFF
            WRITE(32, '(A, T10, A, F6.2)') 'CFX2', '=',  CFX2
            CFY2  = CFY1 + YSLOP*CDIFF
            WRITE(32, '(A, T10, A, F6.2)') 'CFY2', '=',  CFY2
            PX1   = LENGTH / 6 *(2*CFX1  + CFX2)
            WRITE(32, '(A, T10, A, F6.2)') 'PX1', '=', PX1
            PX2   = LENGTH / 6 *(CFX1  + 2*CFX2)
            WRITE(32, '(A, T10, A, F6.2)') 'PX2', '=', PX2
            PY1   = LENGTH / 6 *(2*CFY1  + CFY2)
            WRITE(32, '(A, T10, A, F6.2)') 'PY1', '=', PY1
            PY2   = LENGTH / 6 *(CFY1  + 2*CFY2)
            WRITE(32, '(A, T10, A, F6.2)') 'PY2', '=', PY2
            RI = NODE(N)
            WRITE(32, *) 'RI=', RI
            IF ( IDOF(RI, 1).EQ.-1 ) GO TO 901
            FG( IDOF(RI, 1) ) = FG( IDOF(RI, 1) ) + PX1
  901       IF ( IDOF(RI, 2).EQ.-1 ) GO TO 902
            FG( IDOF(RI, 2) ) = FG( IDOF(RI, 2) ) + PY1
  902       RI = NODE(N+1)
            WRITE(32, *) 'RI=', RI
            IF ( IDOF(RI, 1).EQ.-1 ) GO TO 903
            FG( IDOF(RI, 1) ) = FG( IDOF(RI, 1) ) + PX2
  903       IF ( IDOF(RI, 2).EQ.-1 ) GO TO 904
            FG( IDOF(RI, 2) ) = FG( IDOF(RI, 2) ) + PY2
            CFX1 = CFX2
  904       M = M + 4
            N = N + 2
   30    CONTINUE     
   10 CONTINUE  
      GO TO 99
     
   99 RETURN
C
C     .. End of SURFL ..
C
      END
