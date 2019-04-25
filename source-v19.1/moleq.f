      SUBROUTINE MOLEQ(T,PE,G2,XIH,XKHM,XIHM,XNENH,F1,F2,F3,F4,F5,FE,
     *                 FSUM,EH)
C
C        THIS ROUTINE COMPUTES DISSOCIATION EQUILIBRIA FOR THE MOLECULES
C        H2 AND H2+ WITH H+, H AND H- CONSIDERED. IT MAINLY FOLLOWS MIHALAS,
C        METH. COMP. PHYS. 7, 1 (1967).
C
C        THE INNER ENERGY OF THE HYDROGEN GAS, EH, IS ALSO EVALUATED.
C
C        XIH=THE IONIZATION ENERGY OF HYDROGEN
C        XKHM=THE 'DISSOCIATION CONSTANT' OF H-
C        XIHM=THE 'DISSOCIATION ENERGY' OF H-.
C        XNENH=THE NUMBER OF ELECTRONS PER UNIT VOLUME FROM ELEMENTS OTHER THAN
C        HYDROGEN (Q IN MIHALAS'S ARTICLE)
C        G2,F1,F2 ETC. SEE REF.
C        DOUBLE PRECISION NECESSARY FOR RELATIVELY LOW PRESSURES.
C
C
      COMMON/UTPUT/IREAD,IWRIT
      DOUBLE PRECISION G3,G4,G5,A,E,B,C,D,C1,C2,C3,CAM,F1D,F2D,F3D,F4D,
     *F5D,FED,FSUMD,ROOT
C
C        CALL MOLFYS FOR PHYSICAL DATA
      CALL MOLFYS(T,XKH2,XKH2P,DEH2,DEH2P,deh2nodis,deh2pnodis)
C
C        CALCULATION OF THE EQUILIBRIUM
      G3=PE/XKHM
      G4=PE/XKH2P
      G5=PE/XKH2
      A=1.+G2+G3
      E=G2*G4/G5
      B=2.*(1.+E)
      C=G5
      D=G2-G3
      C1=C*B*B+A*D*B-E*A*A
      C2=2.*A*E-D*B+A*B*XNENH
      C3=-(E+B*XNENH)
      CAM=C2/(2.*C1)
      ROOT=DSQRT(CAM*CAM-C3/C1)
      F1D=-CAM+ROOT
      IF(F1D.GT.1.D0)F1D=-CAM-ROOT
      F5D=(1.D0-A*F1D)/B
      F4D=E*F5D
      F3D=G3*F1D
      F2D=G2*F1D
      FED=F2D-F3D+F4D+XNENH
      FSUMD=F1D+F2D+F3D+F4D+F5D
      F1=F1D
      F2=F2D
      F3=F3D
      F4=F4D
      F5=F5D
      FE=FED
      FSUM=FSUMD
C
C        CALCULATION OF THE ENERGIES
      EH2=(-2.*XIH+DEH2)*F5
      EH2P=(-XIH+DEH2P)*F4
      EHM=-(XIHM+XIH)*F3
      EHJ=-XIH*F1
      EH=EHJ+EHM+EH2+EH2P
    1 CONTINUE
      RETURN
      END
