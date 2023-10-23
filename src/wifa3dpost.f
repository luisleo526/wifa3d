C
C wifa3d FORTRAN BINARY OUTPUT -> TECPLOT ASCII FORMAT
C CSW at NTUESOE 
C
C#############################################################
      PROGRAM COLLECT_SOLUTION
C#############################################################
      INCLUDE 'param.inc'
      PARAMETER (NG=NXT*NY*NZ)
      CHARACTER*256 WORD,W1,W2,CTMP
      DIMENSION LI(NXT),LJ(NY)
      DIMENSION UAVE(NXT),VAVE(NXT),WAVE(NXT),VUFM(NXT),WUFM(NXT)
      DIMENSION FLUX(NXT),FLVX(NXT),FLWX(NXT),AREA(NXT)
      DIMENSION X(NXT),Y(NY),Z(NZ),F1(NG),F2(NXYZ),F3(NXYZ)
      DIMENSION XC(NXT),YC(NY),ZC(NZ),U(NXYZ),V(NXYZ),W(NXYZ)
     *         ,P(NXYZ),SHR(NXYZ),VIS(NXYZ),TE(NXYZ),ED(NXYZ),VIST(NXYZ)
     *         ,BDY1(NXYZ),BDY2(NXYZ),BDY3(NXYZ),VOL(NXYZ),WTE(NXYZ)
      DIMENSION ZXC(NG),ZYC(NG),ZZC(NG),ZU(NG),ZV(NG),ZW(NG)
     *         ,ZP(NG),ZSHR(NG),ZVIS(NG),ZTE(NG),ZED(NG),ZVIST(NG)
     *         ,ZBDY1(NG),ZBDY2(NG),ZBDY3(NG),ZWTE(NG)
     *         ,TX(MWT),TY(MWT),TZ(MWT)
      DOUBLE PRECISION X,Y,Z,XC,YC,ZC,U,V,W,P,SHR,VIS,TE,ED,VIST,WTE
     *                ,BDY1,BDY2,BDY3,VOL,F1,F2,F3,AREA,RC,PC,TMA
     *                ,XDIST,TDS,FLUX,FLVX,FLWX,UAVE,VAVE,WAVE,SMALL
     *                ,RYU,RYL,RZU,RZL,HH,VUFM,WUFM,RD,TMP,TX,TY,TZ
      DOUBLE PRECISION ZXC,ZYC,ZZC,ZU,ZV,ZW,ZP,ZSHR,ZVIS,ZTE,ZED,ZVIST
     *                ,ZBDY1,ZBDY2,ZBDY3,ZWTE
      LOGICAL LTMP
      INTEGER LI,LJ
C--------------------------------------------------------------
      HALF=0.5D0
      SMALL=1D-3
      WRITE(WORD,22) '2022/01/01 (Ver:2022R1)'
      WRITE(*,*) '(*) Built on ',TRIM(WORD)
C
      FWDIR="RESULT"
      FRDIR="RESTART"
      OPEN (UNIT=9,FILE="CASE")
      READ(9,*) FCASE
      CLOSE(9)
      WRITE(*,*) '(*) CASE = ',TRIM(FCASE)
      IF(LEN_TRIM(FCASE).GT.124)THEN
        PRINT *,'(E) CASE NAME TOO LONG!'
        STOP 
      ENDIF
C
      FILWTD=TRIM(FCASE)//'.wtd'
      FILOUT=TRIM(FCASE)//'a.plt'
      FILCHK=TRIM(FCASE)//'ave.plt'
C
      OPEN (UNIT=5,FILE=TRIM(FILWTD),ACTION='READ',IOSTAT=IERR)
c     write(*,*) IERR
      IF(IERR.NE.0)THEN
        write(*,*) "(E) ",TRIM(FILWTD),": file not found!"
        STOP
      ENDIF
      READ(5,11) CTMP
      READ(5,*)  LTMP,LTMP,LTMP,LTMP
      READ(5,11) CTMP
      READ(5,*)  NWT,TMP
      IF(NWT.EQ.1)THEN
        READ(5,11) CTMP
        READ(5,*)  (TX(I),TY(I),TZ(I),I=1,NWT),TMP
        READ(5,11) CTMP
        READ(5,11) TMP
        READ(5,11) CTMP
        READ(5,11) CTMP
        READ(5,11) CTMP
        READ(5,*)  RD,HH,TMP,TMP,TMP,TMP
        CLOSE(5) 
        WRITE(*,2) 'RD=',RD,' HH=',HH
      ELSE
        WRITE(*,*) '(*) Wind turbine array mode'
      ENDIF
   2  FORMAT(1X,A3,E10.2,A4,E10.2)
C
      FILRES=TRIM(FWDIR)//'/'//TRIM(FCASE)//'.1' 
      OPEN (UNIT=4,FILE=TRIM(FILRES),FORM='UNFORMATTED'
     *     ,ACTION='READ',IOSTAT=IERR)
c     write(*,*) IERR
      IF(IERR.NE.0)THEN
        write(*,*) "(E) ",TRIM(FILRES),": file not found!"
        STOP
      ENDIF
      READ(4) NPRO,IID
      CLOSE(4)
C
      DO I=1,NXT
        LI(I)=(I-1)*NY*NZ
      END DO
      DO J=1,NY
        LJ(J)=(J-1)*NZ
      END DO
C
      DO N=0,NPRO-1
        MYID=N
        FILRES=TRIM(FWDIR)//'/'//TRIM(FCASE)
C
        CALL APPEND_INDEX(MYID+1,FILRES)
        WRITE(*,*) 'READING FILE: ',TRIM(FILRES)
C
        OPEN (UNIT=4,FILE=TRIM(FILRES),FORM='UNFORMATTED'
     *       ,ACTION='READ',IOSTAT=IERR)
c       write(*,*) IERR
        IF(IERR.NE.0)THEN
          write(*,*) "(E) ",TRIM(FILRES),": file not found!"
          STOP
        ENDIF
        READ(4) IPRO,IID,NI,NJ,NK,NII,NIJ,NJK,NIK,NIJK
        IF(IPRO.NE.NPRO)THEN
          CLOSE(4)
          PRINT *,'NO. OF PARTITIONS DISMATCHED. (',IPRO,')'
          STOP
        ENDIF
        IF(IID.NE.MYID)THEN
          CLOSE(4)
          PRINT *,'PARTITION NUMBER DISMATCHED. (',IID,')'
          STOP
        ENDIF
        READ(4) (WTE(IJK),IJK=1,NIJK)
        READ(4) ( X(I),I=1,NI),( Y(J),J=1,NJ),( Z(K),K=1,NK),
     *          (XC(I),I=1,NI),(YC(J),J=1,NJ),(ZC(K),K=1,NK),
     *          ( F1(IJK),IJK=1,NIJK),( F2(IJK),IJK=1,NIJK),
     *          ( F3(IJK),IJK=1,NIJK),(  U(IJK),IJK=1,NIJK),
     *          (  V(IJK),IJK=1,NIJK),(  W(IJK),IJK=1,NIJK),
     *          (  P(IJK),IJK=1,NIJK),(VIS(IJK),IJK=1,NIJK),
     *          (SHR(IJK),IJK=1,NIJK),(VOL(IJK),IJK=1,NIJK),
     *          ( TE(IJK),IJK=1,NIJK),( ED(IJK),IJK=1,NIJK),
     *          (VIST(IJK),IJK=1,NIJK),(BDY1(IJK),IJK=1,NIJK),
     *          (BDY2(IJK),IJK=1,NIJK),(BDY3(IJK),IJK=1,NIJK)
        CLOSE(4)
C
        IF(NXT.NE.NI.OR.NJ.NE.NY.OR.NK.NE.NZ)THEN
          WRITE(*,*) "(E) NI,NJ,NK inconsistent!"
          STOP
        ENDIF
        NIM=NI-1
        NJM=NJ-1
        NKM=NK-1
        NIMM=NIM-1
        NJMM=NJM-1
        NKMM=NKM-1
C
        IZBN=MYID*(NI-2)/NPRO+2
        IZED=(MYID+1)*(NI-2)/NPRO+1
        IXBN=IZBN
        IXED=IZED
        IF(IXBN.EQ.2) IXBN=1 
        IF(IXED.EQ.NIM) IXED=NI
        write(*,13) "(",MYID,")=(",IXBN,",",IXED,")"
C
        DO I=IXBN,IXED
        DO J=1,NJ
        DO K=1,NK
           II=I-IZBN+2
           KKK=LI(I )+LJ(J)+K
           KK =LI(II)+LJ(J)+K
           ZXC(KKK)=XC(I)
           ZYC(KKK)=YC(J)
           ZZC(KKK)=ZC(K)
           ZU(KKK)=U(KK)
           ZV(KKK)=V(KK)
           ZW(KKK)=W(KK)
           ZP(KKK)=P(KK)
           ZSHR(KKK)=SHR(KK)
           ZVIS(KKK)=VIS(KK)
           ZTE(KKK)=TE(KK)
           ZED(KKK)=ED(KK)
           ZVIST(KKK)=VIST(KK)
           ZBDY1(KKK)=BDY1(KK)
           ZBDY2(KKK)=BDY2(KK)
           ZBDY3(KKK)=BDY3(KK)
           ZWTE(KKK)=WTE(KK)
 111    END DO
        END DO 
        END DO 
      END DO
C
CCCCCC TECPLOT OUTPUT
C
      WRITE(*,*) 'WRITING FILE: ',TRIM(FILOUT)
      OPEN (UNIT=11,FILE=TRIM(FILOUT),FORM='FORMATTED',ACTION='WRITE'
     *     ,IOSTAT=IERR)
c     write(*,*) IERR
      IF(IERR.NE.0)THEN
        write(*,*) "(E) ",TRIM(FILRES),": file open error!"
        STOP
      ENDIF
      WORD='VARIABLES="X","Y","Z","U","V","W",'//
     *'"P","SHR","VIS","TE","ED","VIST","BX","BY","BZ","WTE"'
      WRITE(11,*) TRIM(WORD)
      WRITE(WORD,10) 
     *'ZONE T="N=',IPRO,',ID=',MYID,'",F=POINT','I=',NI,'J=',NJ,'K=',NK
      WRITE(11,*) TRIM(WORD)
 10   FORMAT(A10,I3,A4,I3,A9,3(1X,A2,I4))
C
      DO K=1,NK
      DO J=1,NJ
      DO I=1,NI
         M=K+LJ(J)+LI(I)
         WRITE(11,12)
     *   ZXC(M),ZYC(M),ZZC(M),ZU(M),ZV(M),ZW(M)
     *  ,ZP(M),ZSHR(M),ZVIS(M),ZTE(M),ZED(M),ZVIST(M)
     *  ,ZBDY1(M),ZBDY2(M),ZBDY3(M),ZWTE(M)
      END DO
      END DO
      END DO
C
 123  CLOSE(11)
C
      IF(NWT.EQ.1)THEN                                                  !YC
C
CCCCCC TECPLOT OUTPUT: WIND AVERAGE
C
      WRITE(*,*) 'WRITING VAVE FILE: ',TRIM(FILCHK)
      OPEN (UNIT=12,FILE=TRIM(FILCHK),FORM='FORMATTED')
      WORD='VARIABLES="X","UAVE","VAVE","WAVE","AREA"'
      WRITE(12,*) WORD
      WRITE(12,*) 'ZONE T="AVERAGE WIND SPEED",F=POINT,I=',NI-2
      RYL=TY(1)-RD*HALF
      RYU=TY(1)+RD*HALF
      RZL=TZ(1)-RD*HALF+HH
      RZU=TZ(1)+RD*HALF+HH
      DO I=2,NI-1
       FLUX(I)=ZERO
       FLVX(I)=ZERO
       FLWX(I)=ZERO
       AREA(I)=ZERO
C
       JB=INDEXFIND(NY,Y,RYU)+1
       JA=INDEXFIND(NY,Y,RYL)+1
       KB=INDEXFIND(NZ,Z,RZU)+1
       KA=INDEXFIND(NZ,Z,RZL)+1
C
       DO J=JA,JB
       DO K=KA,KB
         IJK=LI(I)+LJ(J)+K
         DY=YC(J)-TY(1)
         DZ=ZC(K)-TZ(1)-HH
         TMP=2*SQRT(DY*DY+DZ*DZ)
         IF(TMP.LE.RD)THEN
           TMP=(Z(K)-Z(K-1))*(Y(J)-Y(J-1))
           TMA=TMP*ZU(IJK)
           AREA(I)=AREA(I)+TMP
           FLUX(I)=FLUX(I)+TMA
           FLVX(I)=FLVX(I)+TMA*ZU(IJK)
           FLWX(I)=FLWX(I)
     &          +TMA*(ZU(IJK)*ZU(IJK)+ZV(IJK)*ZV(IJK)+ZW(IJK)*ZW(IJK))
         ENDIF
       END DO
       END DO
       UAVE(I)=FLUX(I)/AREA(I)
       VAVE(I)=FLVX(I)/FLUX(I)
       WAVE(I)=SQRT(FLWX(I)/FLUX(I))
C
       IF(VAVE(I).LT.SMALL) VAVE(I)=UNITY
       IF(WAVE(I).LT.SMALL) WAVE(I)=UNITY
       WRITE(12,*) XC(I),UAVE(I),VAVE(I),WAVE(I)
     *            ,AREA(I)
      END DO
      CLOSE(12)
C
      ENDIF                                                             
C
 11   FORMAT(A256)
 12   FORMAT(1X,19(1X,E15.9))
 13   FORMAT(1X,A1,I3,A3,I4,A1,I4,A1)
 22   FORMAT(A23)
 23   FORMAT(A43,I6)
      END
C
C##########################################################
      FUNCTION XDIST(PA,PB)
C##########################################################
C
      DOUBLE PRECISION XDIST,TMP,PA(3),PB(3)
C
      XDIST=0.0
      DO I=1,3
        TMP=PA(I)-PB(I)
        XDIST=XDIST+TMP*TMP
      END DO 
      XDIST=SQRT(XDIST)
      RETURN
      END
C##########################################################
      SUBROUTINE APPEND_INDEX(ID,FNAME)
C##########################################################
      CHARACTER*64 FNAME
      DIMENSION IDX(3)
C
      N=ID
      J=LEN_TRIM(FNAME)+1
      FNAME(J:J)='.'
      DO I=3,1,-1
        IDX(I)=0
      END DO
      I=0
 222  IF(N.GT.0)then
        I=I+1
        IDX(I)=N-(N/10)*10
        N=(N-IDX(I))/10
        GOTO 222
      ENDIF
      DO K=1,I
        M=J+I-K+1
        FNAME(M:M)=ACHAR(IDX(K)+48)
      END DO
      RETURN
      END
C############################################################
      FUNCTION INDEXFIND(NP,XP,TMP) 
C############################################################
      INTEGER INDEXFIND,NP
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
