C     ..
C     .. Parameters ..
C     ..
      PARAMETER        ( MAXND = 3000,
     ;                   MAXEL = 3000,
     ;                   MAXRT = 20,
     ;                   MAXDF = 2*MAXND,
     ;                   MAXSL = 20,
     ;                   MAXVL = 20,
     ;                   MAXNLD= 20,
     ;                   MAXBE = 1000,
     ;                   MAXBN = 1000 )
C     ..
C     .. Scalar Variables ..
C     ..
      INTEGER*4        LDCRD, LDCON, LDIDOF, LDVLDS, LDSLDS, LEDRST,
     ;                 LDPRST, LDKG, LDLBEL, LDRBEL, LDLNOD, LDRNOD,
     ;                 NDOF, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                 NPRST, NNDL, NSL, NVL, NNODE, ETYPE, PTYPE, NBCO,
     ;                 NECO
C     ..
C     .. Array Variables ..
C     ..
      INTEGER*4        CON(MAXEL, 8), IDOF(MAXND, 2), ERST(MAXRT, 2),
     ;                 PRST(MAXRT, 2), RBEL(MAXBE, 4), LBEL(MAXBE, 4), 
     ;                 RNOD(MAXBN, 2), LNOD(MAXBN, 2), DOUT(50)
     
      REAL*8           CRD(MAXND, 2), CMAT(4), D(3, 3),
     ;                 SLDS(MAXSL, 7), NLDS(MAXDF), VLDS(MAXVL, 2),
     ;                 DSPG(MAXDF), DSPL(12), KG(MAXDF, MAXDF),FG(MAXDF)
      COMMON /CONFIG/ NNODE, NELE, NRBEL, NLBEL, NRNOD, NLNOD, NERST,
     ;                NPRST, NDOF, NNDL, NSL, NVL, ETYPE, PTYPE, NBCO,
     ;                NECO
C 
C     .. Executable Statements .. 
C
C     Open all debugging and output files
C
      OPEN(UNIT=24, FILE='../debug/debug_loads.log')
      OPEN(UNIT=26, FILE='../debug/debug_assemb.log')
      OPEN(UNIT=28, FILE='../debug/debug_ISOQ4.log')
      OPEN(UNIT=30, FILE='../debug/debug_T3.log')
      OPEN(UNIT=32, FILE='../debug/debug_surfl.log')
      OPEN(UNIT=34, FILE='../debug/stiff.dat')
      OPEN(UNIT=36, FILE='../debug/force.dat')
      OPEN(UNIT=42, FILE='../out/disp.out')
      OPEN(UNIT=44, FILE='../out/mesh.vtk')
      OPEN(UNIT=46, FILE='../out/deformed.vtk')
      OPEN(UNIT=48, FILE='../out/sigxx.vtk')
      OPEN(UNIT=50, FILE='../out/sigyy.vtk')
      OPEN(UNIT=52, FILE='../out/sigxy.vtk')
      OPEN(UNIT=54, FILE='../out/sigvm.vtk')
      OPEN(UNIT=56, FILE='../out/sigmax.vtk')
      OPEN(UNIT=58, FILE='../out/sigmin.vtk')
      OPEN(UNIT=60, FILE='../out/ormax.vtk')
      OPEN(UNIT=62, FILE='../out/ormin.vtk')
      OPEN(UNIT=92, FILE='../debug/T3_coeff.dat', FORM='UNFORMATTED')
C
C     Initialize the leading dimensions
C
      LDCRD = MAXND   
      LDCON = MAXEL  
      LDLBEL= MAXBE
      LDRBEL= MAXBE
      LDLNOD= MAXBN
      LDRNOD= MAXBN
      LDIDOF= MAXND   
      LDSLDS= MAXSL
      LDVLDS= MAXVL
      LDERST= MAXRT   
      LDPRST= MAXRT   
      LDKG  = MAXDF
C
C     Initialize the arrays
C
      CON  = 0
      LBEL = 0
      RBEL = 0
      LNOD = 0
      RNOD = 0
      IDOF = 0
      ERST = 0
      PRST = 0
      DOUT = 0
C
      CRD  = 0.D0
      CMAT = 0.D0
      SLDS = 0.D0
      NLDS = 0.D0
      VLDS = 0.D0
      DSPG = 0.D0
      DSPL = 0.D0
      D    = 0.D0
      KG   = 0.D0
      FG   = 0.D0
C
C     Call the subroutines
C     
      PRINT *,
      PRINT *, "Starting GEOMET ..."
      CALL GEOMET(CRD, LDCRD, CON, LDCON, LBEL, LDLBEL, RBEL, LDRBEL,
     ;            LNOD, LDLNOD, RNOD, LDRNOD, CMAT, ERST, LDERST, 
     ;            PRST, LDPRST, DOUT)
      PRINT *, "GEOMET finished successfully"
      PRINT *, 
      PRINT *, "Starting SCODE ..."
      CALL SCODE(IDOF, LDIDOF, RBEL, LDRBEL, RNOD, LDRNOD, ERST, LDERST,
     ;           PRST, LDPRST)
      PRINT *, "SCODE finished successfully"
      PRINT *, 
      PRINT *, "Starting LOADS ..."
      CALL LOADS(IDOF, LDIDOF, NLDS, SLDS, LDSLDS, VLDS, LDVLDS, LNOD,
     ;           LDLNOD)
      PRINT *, "LOADS finished successfully"
      PRINT *, 
      PRINT *, "Starting ASSEMB ..."
      CALL ASSEMB(CON, LDCON, IDOF, LDIDOF, CRD, LDCRD, LBEL, LDLBEL,
     ;            LNOD, LDLNOG, CMAT, SLDS, LDSLDS, NLDS, VLDS, LDVLDS,
     ;            KG, LDKG, FG, D)
      PRINT *, "ASSEMB finished successfully"
      PRINT *, 
C      DO 80 I = 1, NDOF
C         WRITE(34,*)
C         DO 90 J = 1, NDOF
C            WRITE(34, '(E26.16)', ADVANCE='NO') KG(I, J)
C   90    CONTINUE
C   80    CONTINUE
C      DO 100 I = 1, NDOF
C         WRITE(36, '(F26.16)') FG(I)
C  100 CONTINUE
      PRINT *, "Starting GSOLVE ..."
      CALL GSOLVE(KG, LDKG, NDOF, NDOF, DSPG, FG)
      PRINT *, "GSOLVE finished succesfully"
      PRINT *, 
C      WRITE(42, '(/A//)') 'GLOBAL DISPLACEMENTS'
C      DO 110 I = 1, NDOF
C         WRITE(42, '(F18.10)') DSPG(I)
C  110 CONTINUE
      PRINT *, "Starting POST ..."
      CALL POST(CON, LDCON, CRD, LDCRD, IDOF, LDIDOF, ELDS, LDELDS, 
     ;          DSPG, D, DOUT)
      PRINT *, "POST finished succesfully"
      PRINT *, 
      PRINT *, "Program stopped without errors"
      PRINT *,
      STOP
      END

