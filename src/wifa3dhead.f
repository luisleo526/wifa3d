C
C     WIFA3D 5.6 PARAMETER HEAD CHANGED WITH CPU NUMBER
C                                                   CSW
C                                             2021/10/3
C
      PROGRAM WIFA3D_HEAD
C
      PARAMETER (MPRO=128,MWTD=10,MWTN=10000,MX=10000,MY=10000,MZ=1000)
C
      INTEGER I,J,K,ITMP,IT,IRA,IAZ,IWT
      INTEGER NPRO,NI,NJ,NK,NX,NY,NZ,NXT
      INTEGER MTY,MWT,MNP,MCPU,NTYD,NWT,MTT,NAZ,NRA,MRA,MAZ
      INTEGER NAJ,NAK,MAJ,MAK,MAT
      INTEGER NWID(MWTD),IWTY(MWTN)
      INTEGER IIX(MWTN),IJL(MWTN),IJU(MWTN),IKL(MWTN),IKU(MWTN)
      DOUBLE PRECISION TMP,ROTD(MWTD),HALF,HUBH(MWTD)
      DOUBLE PRECISION WBX(MWTN),WBY(MWTN),WBZ(MWTN)
      DOUBLE PRECISION YE,YW,ZT,ZB,DR,RO(3),ROTR
      DOUBLE PRECISION X(MX),Y(MY),Z(MZ),XC(MX),YC(MY),ZC(MZ)
      LOGICAL LTMP
      CHARACTER*128 FCASE,FILGRD,FILWTD,FILPRO,FILHED
      CHARACTER*256 CTMP
      CHARACTER*13  SNX,SNY,SNZ
C
      HALF=0.5D0
      OPEN (UNIT=9,FILE="CASE")
      READ(9,*) FCASE
      CLOSE(9)
C
      FILGRD=TRIM(FCASE)//'.grd'
      FILWTD=TRIM(FCASE)//'.wtd'
      FILPRO=TRIM(FCASE)//'.pro'
      FILHED='head.inc'
C
      OPEN (UNIT=1,FILE=TRIM(FILPRO))
      READ(1,*) NPRO
      IF(NPRO.LT.1) NPRO=1
      IF(NPRO.GT.MPRO) NPRO=MPRO
      CLOSE(1)
      OPEN (UNIT=1,FILE=TRIM(FILPRO))
      WRITE(1,*) NPRO
      CLOSE(1)
      MCPU=NPRO
      IF(MCPU.GT.MPRO) CALL ERROR("NPRO","MPRO",NPRO,MPRO)
C
      OPEN (UNIT=1,FILE=TRIM(FILGRD))
      READ(1,*) NI
      READ(1,*) NJ
      READ(1,*) NK
      IF(NX.GT.MX) CALL ERROR("NX","MX",NX,MX)
      IF(NY.GT.MY) CALL ERROR("NY","MY",NY,MY)
      IF(NZ.GT.MZ) CALL ERROR("NZ","MZ",NZ,MZ)
      I=(NI-2)/NPRO
      J=(NI-2)-I*NPRO
      IF(J.GT.0) J=1
      NX=I+J+2
      NXT=NI
      NY=NJ
      NZ=NK
      NXTM=NXT-1
      NYM=NJ-1
      NZM=NK-1
      READ(1,*) (X(I),I=1,NI)
      READ(1,*) (Y(J),J=1,NJ)
      READ(1,*) (Z(K),K=1,NK)
      CLOSE(1)
C
C.....X- COORDINATES OF CV-CENTERS
C
      DO I=2,NXTM
        XC(I)=HALF*(X(I)+X(I-1))
      END DO
      XC(1) =X(1)
      XC(NXT)=X(NXT)
C
C.....Y- COORDINATES OF CV-CENTERS
C
      DO J=2,NYM
        YC(J)=HALF*(Y(J)+Y(J-1))
      END DO
      YC(1) =Y(1)
      YC(NY)=Y(NY)
C
C.....Z- COORDINATES OF CV-CENTERS
C
      DO K=2,NZM
        ZC(K)=HALF*(Z(K)+Z(K-1))
      END DO
      ZC(1) =Z(1)
      ZC(NZ)=Z(NZ)
C
   11 FORMAT(A256)
      OPEN (UNIT=5,FILE=TRIM(FILWTD))
      READ(5,11) CTMP 
      READ(5,*)  LTMP,LTMP,LTMP,LTMP
      READ(5,11) CTMP 
      READ(5,*)  NWT,NWTD
      MTY=NWTD
      MWT=NWT
      IF(NWT.GT.MWTN) CALL ERROR("NWT","MWTN",NWT,MWTN)
      IF(NWT.GT.MWTN) CALL ERROR("NWTD","MWTD",NWTD,MWTD)
      READ(5,11) CTMP 
      READ(5,*)  (WBX(I),WBY(I),WBZ(I),IWTY(I),I=1,NWT)
C
      DO I=1,NWTD
      READ(5,11) CTMP 
      READ(5,*)  NWID(I)
      READ(5,11) CTMP 
      READ(5,*)  TMP,TMP,TMP,TMP
      READ(5,11) CTMP 
      READ(5,*)  ROTD(I),HUBH(I),TMP,TMP,TMP,TMP
      READ(5,11) CTMP 
      READ(5,*)  TMP,TMP
      READ(5,11) CTMP 
      READ(5,*)  ITMP,TMP,TMP
      READ(5,11) CTMP 
      READ(5,*) (TMP,TMP,K=1,ITMP)
      READ(5,11) CTMP 
      READ(5,*)  ITMP
      IF(I.eq.1)THEN
        MNP=ITMP 
      ELSE
        IF(ITMP.GT.MNP) MNP=ITMP
      ENDIF
      READ(5,11) CTMP 
      READ(5,*) (TMP,TMP,TMP,TMP,K=1,ITMP)
      END DO
      CLOSE(5)
C
      MRA=0
      MRA=0
      MAJ=0
      MAK=0
      DO IWT=1,NWT
C
        ITY=IWTY(IWT)
        ROTR=ROTD(ITY)*HALF
        RO(1)=WBX(IWT)
        RO(2)=WBY(IWT)
        RO(3)=WBZ(IWT)+HUBH(ITY)
        YW=RO(2)-ROTR
        YE=RO(2)+ROTR
        ZB=RO(3)-ROTR
        ZT=RO(3)+ROTR
        IL=2+ICZ*2
        IU=NXTM-ICZ*2
C
        IF(RO(1).LT.X(IL).OR.RO(1).GT.X(IU))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',I,' is outside X range! '
             STOP
          ENDIF
        ENDIF
        IF(YW.LT.Y(1).OR.YE.GT.Y(NY))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',I,' is outside Y range!'
             STOP
          ENDIF
        ENDIF
        IF(ZB.LT.Z(1).OR.ZT.GT.Z(NZ))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',I,' is outside Z range!'
             STOP
          ENDIF
        ENDIF
C
        IIX(IWT)=INDEXFIND(NXT,XC,RO(1))
        IJU(IWT)=INDEXFIND(NY ,YC,YE)
        IJL(IWT)=INDEXFIND(NY ,YC,YW)
        IKU(IWT)=INDEXFIND(NZ ,ZC,ZT)
        IKL(IWT)=INDEXFIND(NZ ,ZC,ZB)
        NAJ=IJU(IWT)-IJL(IWT)+2
        NAK=IKU(IWT)-IKL(IWT)+2
        IF(MAJ.LT.NAJ) MAJ=NAJ
        IF(MAK.LT.NAK) MAK=NAK
C
        IL=IXBN-1-2*ICZ
        IU=IXED+1+2*ICZ
C
        NRA=(IJU(IWT)-IJL(IWT)+IKU(IWT)-IKL(IWT)+3)/2
        NAZ=(IJU(IWT)-IJL(IWT)+IKU(IWT)-IKL(IWT)+2)*2
        write(*,*) IWT,NRA,NAZ,NAJ,NAK
        IF(MRA.LT.NRA) MRA=NRA
        IF(MAZ.LT.NAZ) MAZ=NAZ
C
  123 END DO
      MTT=MRA*MAZ
      MAT=MAJ*MAK
C 
      OPEN (UNIT=1,FILE=TRIM(FILHED))
      WRITE(SNX,15) NX
      WRITE(SNY,15) NY
      WRITE(SNZ,15) NZ
      CTMP="      PARAMETER (MX ="//SNX//",NY ="//SNY//",NZ ="//SNZ//","
      WRITE(1,13) TRIM(CTMP)
C     
      WRITE(SNX,15) MTY
      WRITE(SNY,15) MWT
      WRITE(SNZ,15) MNP
      CTMP="     &           MTY="//SNX//",MWT="//SNY//",MNP="//SNZ//","
      WRITE(1,13) TRIM(CTMP)
C     
      WRITE(SNX,15) MCPU
      WRITE(SNY,15) NI
      WRITE(SNZ,15) MTT
      CTMP="     &          MCPU="//SNX//",NXT="//SNY//",MTT="//SNZ//","
      WRITE(1,13) TRIM(CTMP)
C     
      WRITE(SNX,15) MAT 
      WRITE(SNY,15) MAT 
      WRITE(SNZ,15) MAT 
      CTMP="     &           MAT="//SNX//",MBT="//SNY//",MCT="//SNZ//")"
      WRITE(1,13) TRIM(CTMP)
C     
   13 FORMAT(A71)
   14 FORMAT(A53)
   15 FORMAT(I13)
      CLOSE(1)
C
      OPEN (UNIT=1,FILE='run')
      WRITE(1,*) 'mpirun -np ',MCPU,'./wifa3d'
      CLOSE(1)
C
      DO MYID=0,MCPU-1
      IF(MYID.EQ.0)THEN
        IXBN=MYID*(NI-2)/NPRO+2
        IXED=(MYID+1)*(NI-2)/NPRO+1
        NX=IXED-IXBN+3
      ELSE IF(MYID.EQ.(MCPU-1))THEN
        IXBN=MYID*(NI-2)/NPRO+1
        IXED=(MYID+1)*(NI-2)/NPRO
        NX=IXED-IXBN+3
      ELSE
        IXBN=MYID*(NI-2)/NPRO+1
        IXED=(MYID+1)*(NI-2)/NPRO+1
        NX=IXED-IXBN+2
      ENDIF
      write(*,*) NI,NI-2,MYID,IXBN,IXED,NX
      END DO
C
      RETURN
      END
C
      SUBROUTINE ERROR(SA,SB,VA,VB)
      CHARACTER*1 SA,SB
      INTEGER VA,VB
      write(*,*) SA,'=',VA,'>',SB,'=',VB
      STOP
C
      END
C############################################################
      FUNCTION INDEXFIND(NP,XP,TMP)
C############################################################
      INTEGER INDEXFIND,NPM
      DOUBLE PRECISION XP(NP),TMP
C
      NPM=NP-1
C
      DO I=2,NPM
        IF(TMP.LE.XP(I)) GOTO 123
      END DO
      I=NP
C
  123 INDEXFIND=I-1
C
      RETURN
      END

