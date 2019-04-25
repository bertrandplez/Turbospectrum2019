C
      SUBROUTINE TRANEQ
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 3. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 4. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQ'S.
C 79.06.21 *NORD*
C
C
      INCLUDE 'spectrum.inc'
C
      PARAMETER (ITMAX=12000)
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU   
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) 
      COMMON /RHOC/RHO(NDP)
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  NIMPAC,KIMPAC(NRAYS),PIMPAC(NRAYS),MMU(NDP),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON /TRDBUG/IDEBUG
      LOGICAL debug
      DIMENSION A(ITMAX)
      logical hydrovelo,computeIplus
      real velocity
      common/velo/velocity(ndp),hydrovelo,computeIplus

      DATA debug/.false./
C
C INITIATE, XJ IS SET TO THE DIFFUSION LIMIT VALUE
109   IDEBUG=0
      DO 100 K=1,JTAU
      IF (K.GT.1) GO TO 101
      DTAUB=(TAU(2)-TAU(1))*0.5*(X(2)+S(2)+X(1)+S(1))
      DBPLB=(BPLAN(2)-BPLAN(1))/DTAUB
      D2BPL=0.
      GO TO 102
101   IF (K.EQ.JTAU) GO TO 102
      DTAUA=DTAUB
      DTAUB=(TAU(K+1)-TAU(K))*0.5*(X(K)+S(K)+X(K+1)+S(K+1))
      DBPLA=DBPLB
      DBPLB=(BPLAN(K+1)-BPLAN(K))/DTAUB
      DTAUC=0.5*(DTAUA+DTAUB)
      D2BPL=(DBPLB-DBPLA)/DTAUC
102   XH(K)=D2BPL
      XJ(K)=BPLAN(K)+0.333333*(X(K)+S(K))/X(K)*D2BPL
      XK(K)=0.333333*BPLAN(K)+(0.2+0.111111*S(K)/X(K))*D2BPL
      FJ(K)=1.
100   SOURCE(K)=BPLAN(K)
      IF (DEBUG) PRINT 103,X,S,BPLAN,XJ,XH,XK
103   FORMAT('0X,S,B,XJ,D2B,XK='/(/10(1X,1P,10E12.4/)))
C
C CALCULATE THE MATRIX ELEMENTS
      CALL TRRAYS
      CALL TRANFR
      CALL FORMAL
      NIMP1=NIMPAC+1
      IF (DEBUG) PRINT 132,XJ,SOURCE,ERROR,FJ
     & ,((PFEAU(I,K),K=1,NDP),I=1,NIMP1)
      IF (IDEBUG.GT.1) GO TO 150
      if (computeIplus) then
* skip iteration on S, BPz 06/06-2018
        itm=1
        goto 141
      endif
C
C ITERATION LOOP
      DO 110 IT=1,ITMAX
110   A(IT)=0.
      DO 140 IT=1,ITMAX
      ITM=IT
C
C SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
      CALL TRANSC
      CALL SCATTR
      IF (DEBUG) PRINT 122,EX(ISCAT),DUM,P,DTAUS
122   FORMAT('0EX,SP1,SP2,SP3,P,DTAUS=',E10.4/(/10(1X,1P,10E12.4/)))
C
C CORRECTION TO THE SOURCE FUNCTION
      DO 120 K=1,JTAU1
      P(K)=ERROR(K)+P(K)*FJ(K)*S(K)/(X(K)+S(K))
      A(IT)=AMAX1(A(IT),ABS(P(K)/SOURCE(K)))
120   CONTINUE
C
C CHECK ERROR IN SOURCE FUNCTION
      IF (A(IT).LT.0.001) GO TO 141
      DO 130 K=1,JTAU1
130   SOURCE(K)=amax1(SOURCE(K)+P(K),source(k)/2.)
C
C SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
      CALL FORMAL
      NTAU=KIMPAC(ISCAT)
C
C NOTE THAT FJ() SHOULD ONLY BE PICKED UP ABOVE JTAU0.  THE ISCAT
C BECOMES TO INCLINED BELOW JTAU0.
      DO 131 K=1,JTAU0
131   FJ(K)=XJ(K)/PFEAU(ISCAT,K)
      IF (DEBUG) PRINT 132,XJ,SOURCE,ERROR,FJ
     & ,((PFEAU(I,K),K=1,NDP),I=1,NIMP1)
132   FORMAT('0XJ,SO,ERR,FJ,PF='/(/10(1X,1P,10E12.4/)))
      IF (IDEBUG.GT.1) GO TO 150
C
C END OF ITERATION LOOP
140   CONTINUE
C
C NOT CONVERGED
      IDEBUG=1
CCC      WRITE (13) JTAU,TAU,X,S,BPLAN,RADIUS,RR,RHO,ROSS
      WRITE(6,142) (A(IT),IT=1,ITM)
142   FORMAT(' MAXFEL =',12F9.5)
C
C CONVERGED, IF IN FIRST ITERATION, HAVE TO CALCULATE FJ().
141   IF (ITM.GT.1) GO TO 143
      NTAU=KIMPAC(ISCAT)
      DO 144 K=1,NTAU
144   FJ(K)=XJ(K)/PFEAU(ISCAT,K)
143   CONTINUE
C
C CALCULATE MOMENTS, AND CHECK DEBUG CONTROL
      CALL TRMOM
150   IF (DEBUG.AND.IDEBUG.GT.1) STOP
      IF (DEBUG.AND.IDEBUG.EQ.1) IDEBUG=0
      DEBUG=IDEBUG.GT.1
      IF (DEBUG) GO TO 109
C
      RETURN
      END
