C#############################################################
      PROGRAM Wind_Farm_3D
C#############################################################
c     
C     This code solves the 3D wind farm flows using
C     Cartesian grid and colocated variable arrangement.
C     Actuator Disk Approach
C
C     I,X,U : Downstream Direction
C     J,Y,W : Horizontal Direction
C     K,Z,V : Vertical Direction
C
C                                 Chau,S.W. at NTUESOE
C-------------------------------------------------------------
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
      CHARACTER*16 PDATE,PVER
      DOUBLE PRECISION SOURCE,VISTMAX,VISTMIN
C
C--------------------------------------------------------------
C
      WRITE(PDATE,1) '2023/10/17'
      WRITE(PVER ,1) 'Ver.2023R1'
 1    FORMAT(10A)
C
      NYM=NY-1
      NZM=NZ-1
      NYMM=NYM-1
      NZMM=NZM-1
      NXTM=NXT-1
C
      ZERO   =0.0D0
      UNITY  =1.0D0
      HALF   =0.5D0
      TWO    =2.0D0
      THREE  =3.0D0
      TEN    =1.0D1
      SMALL  =1.0D-30
      GREAT  =1.0D+30
      PI     =3.14159265358979324D+0
      DEG2RAD=PI/180.D+0
      TENTH  =UNITY/TEN
      EIGHTH =0.125D0
      QUARTER=0.25D0
      THIRD  =UNITY/THREE
      THREEQR=0.75D0
      SEC2MIN=60.0D0
      ONETHIRD=UNITY/THREE
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPRO,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      IROOT=0
      IPRN=NPRO-1
C                          
      IF(NPRO.NE.MCPU)THEN
       IF(MYID.EQ.IROOT) THEN
         WRITE(*,*) '(E) Total Processor No. != ',MCPU
       ENDIF
       GOTO 530 
      ENDIF
      IF(MYID.EQ.IPRN)THEN
        WRITE(*,*) '(*) ',TRIM(PVER),'(',TRIM(PDATE),')' 
      ENDIF
C
C.....I/O FILE NAMES
C
      TITLE="CASE"
      OPEN (UNIT=9,FILE=TRIM(TITLE),ACTION='READ',IOSTAT=IERR)
      IF(IERR.NE.0) CALL FILEERROR(TITLE)
      READ(9,*) FCASE
      CLOSE(9)
      FWDIR="RESULT"
      FRDIR="RESTART"
C     
      IF(MYID.EQ.IPRN)THEN
        WRITE(*,*) '(*) CASE = ',TRIM(FCASE)
        WRITE(*,*) '(*) DIR  = ',TRIM(FWDIR)
        WRITE(*,*) '(*) DIR  = ',TRIM(FRDIR)
      ENDIF
      IF(LEN_TRIM(FCASE).GT.124)THEN
        IF(MYID.EQ.IPRN)THEN
          PRINT *,'(E) CASE NAME TOO LONG!'
        ENDIF
        GOTO 530
      ENDIF
C
      CALL FINDFILENAME
C
      CALL MODINP
C
      CALL INITMPI(FILRST,FILRES,NPRO)
C
      CALL INITFIELD
C
  116 FORMAT(1X,'(*) Processor Number =',I4)
  111 FORMAT(1X,'(*) Using ADM Module, NWT =',I4)
C
      IF(MYID.EQ.IROOT) WRITE(*,116) NPRO
      IF(LADM)THEN
        IF(MYID.EQ.IROOT) WRITE(*,111) NWT
        CALL INITADM
      ENDIF
c
      IF(MYID.EQ.IPRN)THEN
        OPEN (UNIT=2,FILE=TRIM(FILOUT),IOSTAT=IERR)
        IF(IERR.NE.0) CALL FILEERROR(TRIM(FILOUT))
        OPEN (UNIT=9,FILE=TRIM(FILPWR),IOSTAT=IERR)
        IF(IERR.NE.0) CALL FILEERROR(TRIM(FILPWR))
      ENDIF
C
      IF(NPRO.GT.1) CALL LOADCAL
      CALL INITFP
C
      IF(LCAL(ITE)) CALL MODVIS
C
C.....READ RESULTS OF PREVIOUS RESULTS
C
      IF(LREAD) THEN
          IF(MYID.EQ.IPRN) THEN
            WRITE(*,*) 'READ PREVIOUS RESULTS: ',FILRST
          ELSE IF(MYID.EQ.IROOT)THEN
            WRITE(2,*) 'READ PREVIOUS RESULTS: ',FILRST
          ENDIF
          OPEN (UNIT=3,FILE=TRIM(FILRST),FORM='UNFORMATTED',IOSTAT=IERR)
          write(*,*) MYID,IERR
          IF(IERR.NE.0) CALL FILEERROR(TRIM(FILRST))
          READ(3) INPRO,IMYID,NI,NJ,NK,NX,NIJ,NJK,NIK,NIJK
          IF(INPRO.NE.NPRO)THEN
            PRINT *,"NPRO MISMATCHED (",INPRO,",",NPRO,")"
            GOTO 530
          ENDIF
          IF(IMYID.NE.MYID)THEN
            PRINT *,"MYID MISMATCHED (",IMYID,",",MYID,")"
            GOTO 530
          ENDIF
          READ(3)(WTE(IJK),IJK=1,NXYZ)
          READ(3)( X(I),I=1,NXT),( Y(J),J=1,NY),( Z(K),K=1,NZ),
     *           (XC(I),I=1,NXT),(YC(J),J=1,NY),(ZC(K),K=1,NZ),
     *           ( F1(IJK),IJK=1,NXYZ),( F2(IJK),IJK=1,NXYZ),
     *           ( F3(IJK),IJK=1,NXYZ),(  U(IJK),IJK=1,NXYZ),
     *           (  V(IJK),IJK=1,NXYZ),(  W(IJK),IJK=1,NXYZ),
     *           (  P(IJK),IJK=1,NXYZ),(VIS(IJK),IJK=1,NXYZ),
     *           (SHR(IJK),IJK=1,NXYZ),(VOL(IJK),IJK=1,NXYZ),
     *           ( TE(IJK),IJK=1,NXYZ),( ED(IJK),IJK=1,NXYZ),
     *           (VIST(IJK),IJK=1,NXYZ),(BDY1(IJK),IJK=1,NXYZ),
     *           (BDY2(IJK),IJK=1,NXYZ),(BDY3(IJK),IJK=1,NXYZ)
           CLOSE(3)
      ENDIF
C
      CALL SETBC
      CALL INFLOW
C
C.....DEFINE MONITORING LOCATION (NODE WITH I=IMON, J=JMON, K=KMON)
C
C
      IF(LRATED)THEN
          IF(MYID.EQ.IPRN)
     $      WRITE(2,*) '(*) RATED Calculation'
          IF(MYID.EQ.IROOT)
     $      WRITE(*,*) '(*) RATED Calculation'
      ENDIF
C
      IF(MYID.EQ.IPRN)THEN
        WRITE(2,605) 'ITER','RES(IU)','RES(IV)','RES(IW)','RES(IP)'
     *              ,'RES(ITE)','RES(IED)'
      ENDIF
      IF(MYID.EQ.IROOT)THEN
        WRITE(*,605) 'ITER','RES(IU)','RES(IV)','RES(IW)','RES(IP)'
     *              ,'RES(ITE)','RES(IED)'
      ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C.....OUTER ITERATIONS (SIMPLE RELAXATIONS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMAX=MAXIT           
      DO ITER=1,MAXIT
        IF(LADM)  CALL ADM
        IF(LUVW)  CALL CALCUVW
        IF(LUVW)  CALL OUTBC
        IF(LPRE)  CALL CALCP
        IF(LVIST) CALL CALCKE(ITE,TE)
        IF(LVIST) CALL CALCKE(IED,ED)
        IF(LVIST) CALL MODVIS
        IF(LVIST) CALL FIND(VIST,VISTMAX,VISTMIN)
        CALL SETBC
C
C.....CHECK CONVERGENCE OF OUTER ITERATIONS
C
        DO I=1,NPHI
          RESOR(I)=RESOR(I)/(SNORM(I)+SMALL)
        END DO
C
        TPWR=ZERO
        DO I=1,NWT
          TPWR=TPWR+POWER(I)
        END DO
C
        IF(MYID.EQ.IPRN)THEN
          WRITE(2,606) ITER,RESOR(IU),RESOR(IV),RESOR(IW),RESOR(IP)
     *                 ,RESOR(ITE),RESOR(IED)
C    *                 ,YPMIN,YPMAX,VISTMAX,BFP,URFP,PWR,TPWR,FAC
          WRITE(9,608) ITER,(VSCA(I),POWER(I),RPM(I),I=1,NWT)
        ENDIF
        IF(MYID.EQ.IROOT)THEN
          WRITE(*,606) ITER,RESOR(IU),RESOR(IV),RESOR(IW),RESOR(IP)
     *                 ,RESOR(ITE),RESOR(IED)
C    *                 ,YPMIN,YPMAX,VISTMAX,BFP,URFP,PWR,TPWR,FAC
C    *                 ,(VSCA(I),I=1,NWT)
        ENDIF
C
        LSTOP=.TRUE.
C
        SOURCE=ZERO
        DO I=1,NPHI
          IF(LSORM(I))THEN
            LSTOP=(LSTOP.AND.(RESOR(I).LT.SORMAX)) 
            SOURCE=MAX(SOURCE,RESOR(I)) 
          ENDIF
        END DO
C
      IF(LSTOP)THEN
        K=INT(UNITY/BFINC)+1
        IF(ITER.LT.K) LSTOP=.FALSE.
      ENDIF
C
      IWF=ITER-(ITER/IFRE)*IFRE
      IF(LSTOP) MAXIT=ITER
      IF((IWF.EQ.0).OR.(ITER.EQ.MAXIT))THEN
C
        OPEN (UNIT=4,FILE=TRIM(FILRES),FORM='UNFORMATTED',IOSTAT=IERR)
        IF(IERR.NE.0) CALL FILEERROR(TRIM(FILRES))
        WRITE(4) NPRO,MYID,NXT,NY,NZ,NX,NXY,NYZ,NXZ,NXYZ
        WRITE(4) (WTE(IJK),IJK=1,NXYZ)
        WRITE(4) ( X(I),I=1,NXT),( Y(J),J=1,NY),( Z(K),K=1,NZ),
     *           (XC(I),I=1,NXT),(YC(J),J=1,NY),(ZC(K),K=1,NZ),
     *           ( F1(IJK),IJK=1,NXYZ),( F2(IJK),IJK=1,NXYZ),
     *           ( F3(IJK),IJK=1,NXYZ),(  U(IJK),IJK=1,NXYZ),
     *           (  V(IJK),IJK=1,NXYZ),(  W(IJK),IJK=1,NXYZ),
     *           (  P(IJK),IJK=1,NXYZ),(VIS(IJK),IJK=1,NXYZ),
     *           (SHR(IJK),IJK=1,NXYZ),(VOL(IJK),IJK=1,NXYZ),
     *           ( TE(IJK),IJK=1,NXYZ),( ED(IJK),IJK=1,NXYZ),
     *           (VIST(IJK),IJK=1,NXYZ),(BDY1(IJK),IJK=1,NXYZ),
     *           (BDY2(IJK),IJK=1,NXYZ),(BDY3(IJK),IJK=1,NXYZ)
        CLOSE(4)
C
      ENDIF
C
        IF(SOURCE.GT.SLARGE) GO TO 510
        IF(LSTOP) GOTO 520
C
      END DO
C
      GOTO 520
C
C==============================================================
C......MESSAGE FOR DIVERGENCE 
C==============================================================
C
  510 IF(MYID.EQ.IROOT)THEN 
        WRITE(*,*)' *** TERMINATED - OUTER ITERATIONS DIVERGING ***'
      ENDIF
      GOTO 520
C
C==============================================================
C......FORMAT SPECIFICATIONS
C==============================================================
C
c 600 FORMAT(1X,'ITER.',3X,
c    *'I---------ABSOLUTE RESIDUAL SOURCE SUMS--------I',3X,
c    *'I----FIELD VALUES AT MONITORING LOCATION (',I3,',',I3,',',I3,
c    *')----I',/,2X,'NO.',9X,'U',11X,'V',9X,'MASS',10X,'T',
c    *16X,'U',11X,'V',11X,'P',11X,'T',/)
C
  605 FORMAT(1X,A8,1X,20A11)
  606 FORMAT(1X,I8,1X,30E11.4)
C 607 FORMAT(1X,I8,1X,200E13.6)
  608 FORMAT(1X,I8,1X,2000E13.6)
C
  520 IF(MYID.EQ.IPRN)THEN
      CLOSE(2)
      CLOSE(9)
      ENDIF
  530 CALL MPI_FINALIZE(IERR)
      STOP
      END
C
C###########################################################
      SUBROUTINE CALCUVW
C###########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      DOUBLE PRECISION URFU,URFV,URFW,START,FINISH
      DOUBLE PRECISION S,VISC,D,CP,TMP,GX,GY,GZ
      DOUBLE PRECISION FXE,FXP,DXPE,CE,FLUXE,FUUDS,FUCDS
      DOUBLE PRECISION FYN,FYP,DYPN,CN,FLUXN,FVUDS,FVCDS
      DOUBLE PRECISION FZT,FZP,DZPT,CT,FLUXT,FWUDS,FWCDS
      DOUBLE PRECISION FXI,FXIM,PE,PW,DX
      DOUBLE PRECISION FYJ,FYJM,PN,PS,DY
      DOUBLE PRECISION FZK,FZKM,PT,PB,DZ
C
      GX=GRAVX*DENSIT
      GY=GRAVY*DENSIT
      GZ=GRAVZ*DENSIT
C
      IF(LADM.AND.BFP.LT.UNITY)THEN
        BFP =BFP+BFINC
        URFP=UNITY*URFUVW
      ELSE
        BFP =UNITY
        URFP=UNITY
      END IF
C
C
C.....RECIPROCAL VALUES OF UNDER-RELAXATION FACTORS FOR U, V AND W
C
C     URFU=UNITY/URF(IU)
C     URFV=UNITY/URF(IV)
C     URFW=UNITY/URF(IW)
C
      URFU=UNITY/(URF(IU)*URFP)
      URFV=UNITY/(URF(IV)*URFP)
      URFW=UNITY/(URF(IW)*URFP)
C
C.....SET BOUNDARY PRESSURE (LINEAR EXTRAPOLATION FROM INSIDE)
C
      CALL PBOUND(P)
C
C.....INITIALIZE TEMPORARILY STORED VARIABLES
C
      DO IJK=1,NXYZ
        SU(IJK)=ZERO
        SV(IJK)=ZERO
        SW(IJK)=ZERO
        APU(IJK)=ZERO
        APV(IJK)=ZERO
        APW(IJK)=ZERO
        AE(IJK)=ZERO
        AW(IJK)=ZERO
        AN(IJK)=ZERO
        AS(IJK)=ZERO
        AT(IJK)=ZERO
        AB(IJK)=ZERO
      END DO
C
C==========================================================
C.....FLUXES THROUGH INTERNAL EAST CV FACES 
C==========================================================
C
c     write(*,*) 'MYID=',MYID,IFBN,IFED
      DO I=IFBN,IFED
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO E
C
        II  =I-IREF
        FXE =FX(I)
        FXP =UNITY-FXE
        DXPE=XC(I+1)-XC(I)
C
        DO J=2,NYM 
        DO K=2,NZM 
C
          IJK=LI(II)+LJ(J)+K
          IJE=IJK+NYZ
C
C.....CELL FACE AREA S = DY*DZ
C
          S=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VIS(IJK)*FXP+VIS(IJE)*FXE                                
          D=VISC*S/DXPE
C
C.....EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
C
          FLUXE=F1(IJK)
C
          CE=MIN(FLUXE,ZERO)
          CP=MAX(FLUXE,ZERO)
C
          FUUDS=CP*U(IJK)+CE*U(IJE)
          FVUDS=CP*V(IJK)+CE*V(IJE)
          FWUDS=CP*W(IJK)+CE*W(IJE)
          FUCDS=FLUXE*(U(IJE)*FXE+U(IJK)*FXP)
          FVCDS=FLUXE*(V(IJE)*FXE+V(IJK)*FXP)
          FWCDS=FLUXE*(W(IJE)*FXE+W(IJK)*FXP)
C
C.....COEFFICIENTS AE(P) AND AW(E) DUE TO UDS
C
          AE(IJK)= CE-D
          AW(IJE)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND E DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(IU)*(FUUDS-FUCDS)
          SU(IJE)=SU(IJE)-GDS(IU)*(FUUDS-FUCDS)
          SV(IJK)=SV(IJK)+GDS(IU)*(FVUDS-FVCDS)
          SV(IJE)=SV(IJE)-GDS(IU)*(FVUDS-FVCDS)
          SW(IJK)=SW(IJK)+GDS(IU)*(FWUDS-FWCDS)
          SW(IJE)=SW(IJE)-GDS(IU)*(FWUDS-FWCDS)
C
        END DO
        END DO
C
      END DO
c     DO I =2,NX-1
c     II=I
c     J=2
c     JJ=NY-1
c     K=2
c     KK=NZ-1
c     write(*,*) "I=",I
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,SV,"SV.")
c     CALL PRINTOUT(II,II,J,JJ,K,KK,NXYZ,SW,"SWW")
c     END DO
c     STOP
C
C=========================================================
C.....FLUXES THROUGH INTERNAL NORTH CV FACES 
C=========================================================
C
      DO I=IXBN,IXED
        II=I-IREF
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO N
C
        DO J=2,NYMM
C
          FYN =FY(J)
          FYP =UNITY-FYN
          DYPN=YC(J+1)-YC(J)
C
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          IJN=IJK+NZ
C
C.....CELL FACE AREA S = DX*DZ
C
          S=(X(I)-X(I-1))*(Z(K)-Z(K-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VIS(IJK)*FYP+VIS(IJN)*FYN               
          D=VISC*S/DYPN
C
C.....EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
C
          FLUXN=F2(IJK)
C
          CN=MIN(FLUXN,ZERO)
          CP=MAX(FLUXN,ZERO)
C
          FUUDS=CP*U(IJK)+CN*U(IJN)
          FVUDS=CP*V(IJK)+CN*V(IJN)
          FWUDS=CP*W(IJK)+CN*W(IJN)
          FUCDS=FLUXN*(U(IJN)*FYN+U(IJK)*FYP)
          FVCDS=FLUXN*(V(IJN)*FYN+V(IJK)*FYP)
          FWCDS=FLUXN*(W(IJN)*FYN+W(IJK)*FYP)
C
C.....COEFFICIENTS AN(P) AND AS(N) DUE TO UDS
C
          AN(IJK)= CN-D
          AS(IJN)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND N DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(IV)*(FUUDS-FUCDS)
          SU(IJN)=SU(IJN)-GDS(IV)*(FUUDS-FUCDS)
          SV(IJK)=SV(IJK)+GDS(IV)*(FVUDS-FVCDS)
          SV(IJN)=SV(IJN)-GDS(IV)*(FVUDS-FVCDS)
          SW(IJK)=SW(IJK)+GDS(IV)*(FWUDS-FWCDS)
          SW(IJN)=SW(IJN)-GDS(IV)*(FWUDS-FWCDS)
C
        END DO
        END DO
C
      END DO
C
C=========================================================
C.....FLUXES THROUGH INTERNAL TOP CV FACES 
C=========================================================
C
      DO I=IXBN,IXED
        II=I-IREF
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO T
C
        DO J=2,NYM 
        DO K=2,NZMM
C
          FZT =FZ(K)
          FZP =UNITY-FZT
          DZPT=ZC(K+1)-ZC(K)
C
          IJK=LI(II)+LJ(J)+K
          IJT=IJK+1
C
C.....CELL FACE AREA S = DX*DY
C
          S=(X(I)-X(I-1))*(Y(J)-Y(J-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VIS(IJK)*FZP+VIS(IJT)*FZT               
          D=VISC*S/DZPT
C
C.....EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
C
          FLUXT=F3(IJK)
C
          CT=MIN(FLUXT,ZERO)
          CP=MAX(FLUXT,ZERO)
C
          FUUDS=CP*U(IJK)+CT*U(IJT)
          FVUDS=CP*V(IJK)+CT*V(IJT)
          FWUDS=CP*W(IJK)+CT*W(IJT)
          FUCDS=FLUXT*(U(IJT)*FZT+U(IJK)*FZP)
          FVCDS=FLUXT*(V(IJT)*FZT+V(IJK)*FZP)
          FWCDS=FLUXT*(W(IJT)*FZT+W(IJK)*FZP)
C
C.....COEFFICIENTS AT(P) AND AB(T) DUE TO UDS
C
          AT(IJK)= CT-D
          AB(IJT)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND T DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(IW)*(FUUDS-FUCDS)
          SU(IJT)=SU(IJT)-GDS(IW)*(FUUDS-FUCDS)
          SV(IJK)=SV(IJK)+GDS(IW)*(FVUDS-FVCDS)
          SV(IJT)=SV(IJT)-GDS(IW)*(FVUDS-FVCDS)
          SW(IJK)=SW(IJK)+GDS(IW)*(FWUDS-FWCDS)
          SW(IJT)=SW(IJT)-GDS(IW)*(FWUDS-FWCDS)
C
        END DO
        END DO
C
      END DO
C
C=============================================================
C.....VOLUME INTEGRALS (SOURCE TERMS)
C=============================================================
C
      TMP=BFP
      DO I=IXBN,IXED
        II=I-IREF
        DX=X(I)-X(I-1)
C
        DO J=2,NYM 
          DY=Y(J)-Y(J-1)
C
          DO K=2,NZM 
            DZ=Z(K)-Z(K-1)
            IJK=LI(II)+LJ(J)+K
C
C...... PRESSURE & GRAVITY SOURCE TERM 
C
            FXI =FX(I)
            FXIM=FX(I-1)                                                
            FYJ =FY(J)                                                  
            FYJM=FY(J-1)                                              
            FZK =FZ(K)                                                  
            FZKM=FZ(K-1)                                              
C
            PE=P(IJK+NYZ)*FXI +P(IJK)    *(UNITY-FXI)
            PW=P(IJK)    *FXIM+P(IJK-NYZ)*(UNITY-FXIM)
            PN=P(IJK+NZ) *FYJ +P(IJK)    *(UNITY-FYJ)
            PS=P(IJK)    *FYJM+P(IJK-NZ) *(UNITY-FYJM)
            PT=P(IJK+1)  *FZK +P(IJK)    *(UNITY-FZK)
            PB=P(IJK)    *FZKM+P(IJK-1)  *(UNITY-FZKM)
C
            DPX(IJK)=(PE-PW)/DX
            DPY(IJK)=(PN-PS)/DY
            DPZ(IJK)=(PT-PB)/DZ
C
            SU(IJK)=SU(IJK)+(GX-DPX(IJK))*VOL(IJK)
            SV(IJK)=SV(IJK)+(GY-DPY(IJK))*VOL(IJK)
            SW(IJK)=SW(IJK)+(GZ-DPZ(IJK))*VOL(IJK)
C
C...... BLADE BODY FORCE TERM
C
            IF(LADM)THEN
              SU(IJK)=SU(IJK)+BDY1(IJK)*TMP
              SV(IJK)=SV(IJK)+BDY2(IJK)*TMP
              SW(IJK)=SW(IJK)+BDY3(IJK)*TMP
            ENDIF
C
          END DO
        END DO
      END DO
C
C=============================================================
C.....PROBLEM MODIFICATIONS - BOUNDARY CONDITIONS
C=============================================================
C 
c     DO I =1,1
c     II=NX
c     J =1
c     JJ=NY
c     K=1
c     KK=NZ
c     write(*,*) "I=",I
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     END DO
c     ID=80+MYID
c     J=1
c     JJ=NY
c     K=1
c     KK=NZ
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ, U,"UUU")
c     END DO
c
      CALL BCUVW
C 
c     ID=70+MYID
c     DO I =1,NX
c     II=I
c     J=1
c     JJ=NY
c     K=1
c     KK=NZ
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,APU,"APU")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,U ,"F1U")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,F2,"F2U")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,F3,"F3U")
c     END DO
C
C=============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR U-VELOCITY
C=============================================================
C
      TMP=UNITY-URF(IU)*URFP
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        AP(IJK)=(-AE(IJK)-AW(IJK)-AN(IJK)-AS(IJK)-AT(IJK)-AB(IJK)
     &           +APU(IJK))*URFU
        SU(IJK)=SU(IJK)+TMP*AP(IJK)*U(IJK)
        APU(IJK)=UNITY/AP(IJK)
      END DO
      END DO
      END DO
C
c     ID=80+MYID
c     J=1
c     JJ=NY
c     K=1
c     KK=NZ
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     END DO
C
      CALL SIPSOL(U,IU)
C
c     ID=80+MYID
c     J=1
c     JJ=NY 
c     K=1
c     KK=NZ 
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ, U,"UUU")
c     END DO
c     STOP
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IF(U(IJK).LT.ZERO) U(IJK)=-U(IJK) 
      END DO
      END DO
      END DO
C
C=============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR V-VELOCITY
C=============================================================
C
      TMP=UNITY-URF(IV)*URFP
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        AP(IJK)=(-AE(IJK)-AW(IJK)-AN(IJK)-AS(IJK)-AT(IJK)-AB(IJK)
     &           +APV(IJK))*URFV
        SU(IJK)=SV(IJK)+TMP*AP(IJK)*V(IJK)
        APV(IJK)=UNITY/AP(IJK)
      END DO
      END DO
      END DO
C
      CALL SIPSOL(V,IV)
C
C=============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR W-VELOCITY
C=============================================================
C
      TMP=UNITY-URF(IW)*URFP
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        AP(IJK)=(-AE(IJK)-AW(IJK)-AN(IJK)-AS(IJK)-AT(IJK)-AB(IJK)
     &           +APW(IJK))*URFW
        SU(IJK)=SW(IJK)+TMP*AP(IJK)*W(IJK)
        APW(IJK)=UNITY/AP(IJK)
      END DO
      END DO
      END DO
C
      CALL SIPSOL(W,IW)
C
      IF(NPRO.GT.1)THEN
        CALL EXCHFI(ISIZE,APU,IU+100)
        CALL EXCHFI(ISIZE,APV,IV+100)
        CALL EXCHFI(ISIZE,APW,IW+100)
        CALL EXCHFI(ISIZE,DPX,IU+200)
        CALL EXCHFI(ISIZE,DPY,IV+200)
        CALL EXCHFI(ISIZE,DPZ,IW+200)
      ENDIF
C
c     ID=80+MYID
c     J =1
c     JJ=NY
c     K=1
c     KK=NZ
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,APU,"UUU")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,APV,"VVV")
c     END DO
c     DO I =1,NX
c     II=I
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,APW,"WWW")
c     END DO
c100  STOP
c
      RETURN
      END 
C
C############################################################## 
      SUBROUTINE CALCP 
C##############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION S,D,SUM,PPO,START,FINISH
      DOUBLE PRECISION FXE,FXP,VOLE,DENE,DPXEL,UEL,APUE,DPXE,UE
      DOUBLE PRECISION FYN,FYP,VOLN,DENN,DPYNL,VNL,APVN,DPYN,VN
      DOUBLE PRECISION FZT,FZP,VOLT,DENT,DPZTL,WTL,APWT,DPZT,WT
      DOUBLE PRECISION FXI,FXIM,PPE,PPW
      DOUBLE PRECISION FYJ,FYJM,PPN,PPS
      DOUBLE PRECISION FZK,FZKM,PPT,PPB
      DOUBLE PRECISION DXPE,DYPN,DZPT
      DOUBLE PRECISION DX,DY,DZ
C
C.....INITIALIZE TEMPORARILY STORED VARIABLES
C
      DO IJK=1,NXYZ
        SU(IJK)=ZERO
        AP(IJK)=ZERO
        AE(IJK)=ZERO
        AW(IJK)=ZERO
        AN(IJK)=ZERO
        AS(IJK)=ZERO
        AT(IJK)=ZERO
        AB(IJK)=ZERO
        PP(IJK)=ZERO
      END DO
C
C============================================================
C.....EAST CV FACES (S - AREA, VOLE - VOLUME BETWEEN P AND E)
C============================================================
C
      DO I=IFBN,IFED
C
        II  =I-IREF
        DXPE=XC(I+1)-XC(I)
        FXE =FX(I)
        FXP =UNITY-FXE
C
        DO J=2,NYM 
        DO K=2,NZM 
C
          IJK=LI(II)+LJ(J)+K
          IJE=IJK+NYZ
C
          S=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
          VOLE=(VOL(IJK)+VOL(IJK+NYZ))*HALF
          DENE=DENSIT
          D=DENE*S
C
C.....INTERPOLATED CELL FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
C
          DPXEL=DPX(IJE)*FXE+DPX(IJK)*FXP
          UEL  =  U(IJE)*FXE+  U(IJK)*FXP
          APUE =APU(IJE)*FXE+APU(IJK)*FXP
C
C.....CELL FACE GRADIENT, VELOCITY AND MASS FLUX
C
          DPXE=(P(IJE)-P(IJK))/DXPE
          UE=UEL-APUE*VOLE*(DPXE-DPXEL)
          F1(IJK)=D*UE
C
C.....COEFFICIENTS OF P' EQUATION, AE(P) AND AW(E)
C
          AE(IJK)=-D*APUE*S
          AW(IJE)=AE(IJK)
C
        END DO
        END DO
      END DO
C
C=============================================================
C.....NORTH CV FACES (S - AREA, VOLN - VOLUME BETWEEN P AND N)
C=============================================================
C
      DO I=IXBN,IXED
C
        II=I-IREF
        DO J=2,NYMM
C
          DYPN=YC(J+1)-YC(J)
          FYN=FY(J)
          FYP=UNITY-FYN
C
          DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          IJN=IJK+NZ
C
          S=(X(I)-X(I-1))*(Z(K)-Z(K-1))
          VOLN=(VOL(IJK)+VOL(IJK+NZ))*HALF
          DENN=DENSIT
          D=DENN*S
C
C.....INTERPOLATED CELL-FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
C
          DPYNL=DPY(IJN)*FYN+DPY(IJK)*FYP
          VNL  =  V(IJN)*FYN+  V(IJK)*FYP
          APVN =APV(IJN)*FYN+APV(IJK)*FYP
C
C.....CELL-FACE GRADIENT, VELOCITY AND MASS FLUX
C
          DPYN=(P(IJN)-P(IJK))/DYPN
          VN=VNL-APVN*VOLN*(DPYN-DPYNL)
          F2(IJK)=D*VN
C
C.....COEFFICIENTS OF P' EQUATION, AN(P) AND AS(N)
C
          AN(IJK)=-D*APVN*S
          AS(IJN)=AN(IJK)
C
        END DO
        END DO
      END DO
C
C=============================================================
C.....TOP CV FACES (S - AREA, VOLT - VOLUME BETWEEN P AND T)
C=============================================================
C
      DO I=IXBN,IXED
C
        II=I-IREF
        DO J=2,NYM 
C
          DO K=2,NZMM
C
          DZPT=ZC(K+1)-ZC(K)
          FZT=FZ(K)
          FZP=UNITY-FZT
          IJK=LI(II)+LJ(J)+K
          IJT=IJK+1
C
          S=(X(I)-X(I-1))*(Y(J)-Y(J-1))
          VOLT=(VOL(IJK)+VOL(IJT))*HALF
          DENT=DENSIT
          D=DENT*S
C
C.....INTERPOLATED CELL-FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
C
          DPZTL=DPZ(IJT)*FZT+DPZ(IJK)*FZP
          WTL  =  W(IJT)*FZT+  W(IJK)*FZP
          APWT =APW(IJT)*FZT+APW(IJK)*FZP
C
C.....CELL-FACE GRADIENT, VELOCITY AND MASS FLUX
C
          DPZT=(P(IJT)-P(IJK))/DZPT
          WT=WTL-APWT*VOLT*(DPZT-DPZTL)
          F3(IJK)=D*WT
C
C.....COEFFICIENTS OF P' EQUATION, AT(P) AND AB(T)
C
          AT(IJK)=-D*APWT*S
          AB(IJT)=AT(IJK)
C
        END DO
        END DO
      END DO
C
C===============================================================
C..... SORCE TERM AND COEFFICIENT OF NODE P
C===============================================================
C
      SUM=ZERO
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        SU(IJK)=F1(IJK-NYZ)-F1(IJK)+F2(IJK-NZ)-F2(IJK)+F3(IJK-1)-F3(IJK)
        AP(IJK)=-(AE(IJK)+AW(IJK)+AN(IJK)+AS(IJK)+AT(IJK)+AB(IJK))
        SUM=SUM+SU(IJK)
      END DO
      END DO
      END DO
c      
c     DO I =1,1
c     II=NX
c     J =1
c     JJ=NY
c     K=1
c     KK=NZ
c     write(*,*) "I=",I
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     END DO
c100  STOP
c
C.....SUM MUST BE ZERO IF GLOBAL MASS CONSERVATION IS ASSURED!
C
C     IF(LTEST.AND.(MYID.EQ.IPRN)) WRITE(2,*) '       SUM = ',SUM
C     IF(MYID.EQ.IPRN) WRITE(*,*) '       SUM = ',SUM
C
C===============================================================
C.....SOLVE EQUATIONS SYSTEM FOR P' AND APPLY CORRECTIONS
C===============================================================
C
c     ID=70+MYID
c     DO I =1,NX
c     II=I
c     J=1
c     JJ=NY
c     K=1
c     KK=NZ
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AW,"AWW")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AE,"AEE")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AS,"ASS")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AN,"ANN")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AB,"ABB")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,AT,"ATT")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,SU,"SUU")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,U ,"F1U")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,F2,"F2U")
c     CALL PRINTOUT(ID,I,II,J,JJ,K,KK,NXYZ,F3,"F3U")
c     END DO
c     STOP
C
      CALL SIPSOL(PP,IP)
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,PP,"PUU")
C
C.....CALCULATE PRESSURE CORRECTION AT BOUNDARIES
C
      CALL PBOUND(PP)
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,PP,"PUU")
C
C.....VALUE OF P' AT REFERENCE LOCATION TO BE SUBTRACTED FROM ALL P'
C
      IF(IPPR.EQ.MYID)THEN
        II=IPR-IREF
        IJKPREF=LI(II)+LJ(JPR)+KPR
        PPO=PP(IJKPREF)
      ELSE
        PPO=ZERO
      ENDIF
C
      IF(NPRO.GT.1)THEN
        CALL MPI_BCAST(PPO,1,MPI_DOUBLE_PRECISION,IPPR,MPI_COMM_WORLD,
     *                 IERR)
      ENDIF 
C
C.....CORRECT EAST MASS FLUXES 
C
      DO I=IFBN,IFED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        F1(IJK)=F1(IJK)+AE(IJK)*(PP(IJK+NYZ)-PP(IJK))
      END DO
      END DO
      END DO
C
C.....CORRECT NORTH MASS FLUXES 
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYMM
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        F2(IJK)=F2(IJK)+AN(IJK)*(PP(IJK+NZ)-PP(IJK))
      END DO
      END DO
      END DO
C
C.....CORRECT TOP MASS FLUXES 
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZMM
        IJK=LI(II)+LJ(J)+K
        F3(IJK)=F3(IJK)+AT(IJK)*(PP(IJK+1)-PP(IJK))
      END DO
      END DO
      END DO
C
C.....CORRECT PRESSURE AND VELOCITIES AT CELL CENTER
C
      DO I=IXBN,IXED
        II=I-IREF
        DX=X(I)-X(I-1)
C
        DO J=2,NYM 
          DY=Y(J)-Y(J-1)
        DO K=2,NZM 
          DZ=Z(K)-Z(K-1)
          IJK=LI(II)+LJ(J)+K
C
          FXI =FX(I)
          FXIM=FX(I-1)
          FYJ =FY(J)
          FYJM=FY(J-1)
          FZK =FZ(K)
          FZKM=FZ(K-1)
C
          PPE=PP(IJK+NYZ)*FXI +PP(IJK)    *(UNITY-FXI)
          PPW=PP(IJK)    *FXIM+PP(IJK-NYZ)*(UNITY-FXIM)
          PPN=PP(IJK+NZ) *FYJ +PP(IJK)    *(UNITY-FYJ)
          PPS=PP(IJK)    *FYJM+PP(IJK-NZ) *(UNITY-FYJM)
          PPT=PP(IJK+1)  *FZK +PP(IJK)    *(UNITY-FZK)
          PPB=PP(IJK)    *FZKM+PP(IJK-1)  *(UNITY-FZKM)
C
          U(IJK)=U(IJK)-(PPE-PPW)*DY*DZ*APU(IJK)
          V(IJK)=V(IJK)-(PPN-PPS)*DX*DZ*APV(IJK)
          W(IJK)=W(IJK)-(PPT-PPB)*DX*DY*APW(IJK)
          P(IJK)=P(IJK)+URF(IP)*(PP(IJK)-PPO)
C
        END DO
        END DO
      END DO
C
      CALL PBOUND(P)          
C
      IF(NPRO.GT.1)THEN
        CALL EXCHFI(ISIZE,U,IU+10)
        CALL EXCHFI(ISIZE,V,IV+10)
        CALL EXCHFI(ISIZE,W,IW+10)
        CALL EXCHFI(ISIZE,P,IP+10)
      ENDIF
C
      CALL GRADXYZ(U,DUX,DUY,DUZ)
      CALL GRADXYZ(V,DVX,DVY,DVZ)
      CALL GRADXYZ(W,DWX,DWY,DWZ)
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        SHR(IJK)=TWO*(DUX(IJK)*DUX(IJK)+DVY(IJK)*DVY(IJK)
     *               +DWZ(IJK)*DWZ(IJK))
        SHR(IJK)=SHR(IJK)+(DUY(IJK)+DVX(IJK))*(DUY(IJK)+DVX(IJK))
        SHR(IJK)=SHR(IJK)+(DUZ(IJK)+DWX(IJK))*(DUZ(IJK)+DWX(IJK))
        SHR(IJK)=SHR(IJK)+(DVZ(IJK)+DWY(IJK))*(DVZ(IJK)+DWY(IJK))
      END DO
      END DO
      END DO
C
c     DO I =1,1
c     II=NX
c     J =1
c     JJ=NY
c     K=1
c     KK=NZ
c     write(*,*) "I=",I
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,U ,"UUU")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,V ,"VVV")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,W ,"WWW")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,F1,"FUU")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,F2,"FVV")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,F3,"FWW")
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,P ,"PPP")
c     END DO
c100  STOP
C      
      RETURN
      END 
C
C#############################################################
      SUBROUTINE SIPSOL(FI,IFI)
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION FI(NXYZ),UE(NXYZ),UN(NXYZ),UT(NXYZ),RES(NXYZ)
      DOUBLE PRECISION LW(NXYZ),LS(NXYZ),LB(NXYZ),LPR(NXYZ)
      DOUBLE PRECISION P1,P2,P3,RESL,TMP(1),RSM,RTMP(1)
C     
C.....COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
C
      DO IJK=1,NXYZ
        LPR(IJK)=ZERO
        RES(IJK)=ZERO
        UE(IJK) =ZERO
        UN(IJK) =ZERO
        UT(IJK) =ZERO
        LW(IJK) =ZERO
        LS(IJK) =ZERO
        LB(IJK) =ZERO
      END DO
C
      DO I=2,NXM 
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(I)+LJ(J)+K
          LW(IJK)=AW(IJK)/(UNITY+ALFA*(UN(IJK-NYZ)+UT(IJK-NYZ)))
          LS(IJK)=AS(IJK)/(UNITY+ALFA*(UE(IJK-NZ) +UT(IJK-NZ)))
          LB(IJK)=AB(IJK)/(UNITY+ALFA*(UN(IJK-1)  +UE(IJK-1)))
          P1=ALFA*(LW(IJK)*UN(IJK-NYZ)+LB(IJK)*UN(IJK-1))
          P2=ALFA*(LS(IJK)*UE(IJK-NZ) +LB(IJK)*UE(IJK-1))
          P3=ALFA*(LW(IJK)*UT(IJK-NYZ)+LS(IJK)*UT(IJK-NZ))
          LPR(IJK)=UNITY/(AP(IJK)+P1+P2+p3
     *            -LW(IJK)*UE(IJK-NYZ)-LS(IJK)*UN(IJK-NZ)
     *            -LB(IJK)*UT(IJK-1)+SMALL)
          UN(IJK)=(AN(IJK)-P1)*LPR(IJK)
          UE(IJK)=(AE(IJK)-P2)*LPR(IJK)
          UT(IJK)=(AT(IJK)-P3)*LPR(IJK)
        END DO
        END DO
      END DO
C
C==============================================================
C.....INNER ITERATIONS LOOP
C==============================================================
C
      DO IT=1,NSW(IFI)
        RESL=ZERO
C      
C.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
C
        DO I=2,NXM 
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(I)+LJ(J)+K
          RES(IJK)=SU(IJK)            -AP(IJK)*FI(IJK)
     *            -AN(IJK)*FI(IJK+NZ) -AS(IJK)*FI(IJK-NZ)
     *            -AE(IJK)*FI(IJK+NYZ)-AW(IJK)*FI(IJK-NYZ)
     *            -AT(IJK)*FI(IJK+1)  -AB(IJK)*FI(IJK-1)
          RESL=RESL+ABS(RES(IJK))
C
          RES(IJK)=(RES(IJK)-LS(IJK)*RES(IJK-NZ)-LW(IJK)*RES(IJK-NYZ)
     *                      -LB(IJK)*RES(IJK-1))*LPR(IJK)
C
        END DO
        END DO
        END DO
C
C.....STORE INITIAL RESIDUAL SUM FOR CHECKING CONV. OF OUTER ITER.
C
        IF(NPRO.GT.1)THEN
          RTMP(1)=RESL
          CALL MPI_ALLREDUCE(RTMP,TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     *                       MPI_COMM_WORLD,IERR)
          RESL=TMP(1)
        ENDIF
C
        IF(IT.EQ.1) RESOR(IFI)=RESL
        RSM=RESL/(RESOR(IFI)+SMALL)
c       write(MYID+80,*) MYID,IFI,IT,RSM,RESL
C
C.....BACK SUBSTITUTION AND CORRECTION
C
        DO I=NXM,2,-1
        DO J=NYM,2,-1
        DO K=NZM,2,-1
          IJK=LI(I)+LJ(J)+K
          RES(IJK)=RES(IJK)-UN(IJK)*RES(IJK+NZ)-UE(IJK)*RES(IJK+NYZ)
     *                     -UT(IJK)*RES(IJK+1)
          FI(IJK)=FI(IJK)+RES(IJK)
        END DO
        END DO
        END DO
C
C       CALL CPU_TIME(START)
C
        IF(NPRO.GT.1) CALL EXCHFI(ISIZE,FI,IFI*10)
C
C       CALL CPU_TIME(FINISH)
C       IF(LPRN) WRITE(MYID+10,*) "SIPSOL =",FINISH-START
C
        IF(RSM.LT.SOR(IFI)) GOTO 999
C
      END DO
C
 999  RETURN
      END
C
C#############################################################
      SUBROUTINE PBOUND(FI)
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'logic.inc'
      DOUBLE PRECISION FI(NXYZ)
C
C.....TOP AND BOTTOM BOUNDARIES
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        K=1
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK+1)+(FI(IJK+1)-FI(IJK+2))*FZ(2)
        K=NZ
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK-1)
     *         +(FI(IJK-1)-FI(IJK-2))*(UNITY-FZ(NZMM)) 
      END DO 
      END DO 
C
C.....SOUTH AND NORTH BOUNDARIES
C
      NZ2=2*NZ
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        J=1
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK+NZ)+(FI(IJK+NZ)-FI(IJK+NZ2))*FY(2)
        J=NY
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK-NZ)
     *         +(FI(IJK-NZ)-FI(IJK-NZ2))*(UNITY-FY(NYMM))
      END DO 
      END DO 
C
C..... WEST ANS EAST BOUNDARY
C
      NYZ2=2*NYZ
      IF(MYID.EQ.IROOT)THEN
      DO J=2,NYM 
      DO K=2,NZM 
        II=1
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK+NYZ)+(FI(IJK+NYZ)-FI(IJK+NYZ2))*FX(2)
      END DO
      END DO
      ENDIF
C
      IF(MYID.EQ.IPRN)THEN
      DO J=2,NYM 
      DO K=2,NZM 
        II=NX
        I=II+IREF
        IJK=LI(II)+LJ(J)+K
        FI(IJK)=FI(IJK-NYZ)
     *         +(FI(IJK-NYZ)-FI(IJK-NYZ2))*(UNITY-FX(I-2))
      END DO
      END DO
      ENDIF
C
 111  RETURN
      END
C
C#############################################################
      SUBROUTINE VISWALL(ID,WVIS,YPLWD,VISWD) 
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      DOUBLE PRECISION WVIS,YPLWD,VISWD
C
      IF(ID.NE.IWAL)THEN
        WVIS=VISWD
      ELSE IF(LCAL(ITE).AND.(YPLWD.GT.CTRANS))THEN
        WVIS=VISWD
      ELSE 
        WVIS=VISCOS
      ENDIF
      RETURN
      END
C
C#############################################################
      SUBROUTINE BCUVW
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'wind.inc'
      DOUBLE PRECISION VISC,D,CP
      DOUBLE PRECISION CT,DXY
      DOUBLE PRECISION CN,DXZ
      DOUBLE PRECISION CE,DYZ
C
C.....BOTTOM BOUNDARY : WALL
C
      K=2
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCB(IJ),VISC,YPLB(IJ),VSWB(IJ))
        DXY=(X(I)-X(I-1))*(Y(J)-Y(J-1))
        D=VISC*DXY/(ZC(K)-ZC(K-1))
        AB(IJK)=-D
      END DO
      END DO
C
C.....TOP BOUNDARY : INLET
C
      K=NZM 
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCT(IJ),VISC,YPLT(IJ),VSWT(IJ))
        DXY=(X(I)-X(I-1))*(Y(J)-Y(J-1))
        D=VISC*DXY/(ZC(K+1)-ZC(K))
        CT=MIN(F3(IJK),ZERO)
        AT(IJK)=CT-D
      END DO
      END DO
C
C.....SOUTH BOUNDARY : INLET
C
      J=2
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCS(IK),VISC,YPLS(IK),VSWS(IK))
        DXZ=(X(I)-X(I-1))*(Z(K)-Z(K-1))
        D=VISC*DXZ/(YC(J)-YC(J-1))
        CP=MAX(F2(IJK-NZ),ZERO)
        AS(IJK)=-CP-D
      END DO
      END DO
C
C.....NORTH BOUNDARY : INLET
C
      J=NYM 
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCN(IK),VISC,YPLN(IK),VSWN(IK))
        DXZ=(X(I)-X(I-1))*(Z(K)-Z(K-1))
        D=VISC*DXZ/(YC(J+1)-YC(J))
        CN=MIN(F2(IJK),ZERO)
        AN(IJK)= CN-D
      END DO
      END DO
C
C.....WEST BOUNDARY : INLET
C
      IF(MYID.EQ.IROOT)THEN
        I=2
        II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          JK=(J-1)*NZ+K
          CALL VISWALL(IBCW(JK),VISC,YPLW(JK),VSWW(JK))
          DYZ=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
          D=VISC*DYZ/(XC(I)-XC(I-1))
          CP=MAX(F1(IJK-NYZ),ZERO)
          AW(IJK)=-CP-D
        END DO
        END DO
      ENDIF
C
C.....EAST BOUNDARY : OUTLET
C
      IF(MYID.EQ.IPRN)THEN
        I=NXTM 
        II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          JK=(J-1)*NZ+K
          CALL VISWALL(IBCE(JK),VISC,YPLE(JK),VSWE(JK))
          DYZ=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
          D=VISC*DYZ/(XC(I+1)-XC(I))
          CE=MIN(F1(IJK),ZERO)
          AE(IJK)= CE-D
        END DO
        END DO
      ENDIF
C
      RETURN
      END
C
C#############################################################
      SUBROUTINE OUTBC
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION AR
C
C..... OUTLET: CALCULATE MASS FLUXES AT EAST 
C
      FMOUT=ZERO
      IF(MYID.EQ.IPRN)THEN
        I=NXT
        II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          AR=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
          F1(IJK-NYZ)=DENSIT*AR*U(IJK-NYZ)
          FMOUT=FMOUT+F1(IJK-NYZ)
        END DO
        END DO
C
        FAC=FLOMAS/(FMOUT+SMALL)
c
c       FAC=(FLOMAS-FMOUT)/(DENSIT*(ARTOP+SMALL))
C
C.....CORRECT VELOCITY AND MASS FLUXES AT EAST
C
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          F1(IJK-NYZ)=  F1(IJK-NYZ)*FAC
               U(IJK)=   U(IJK-NYZ)*FAC
               V(IJK)=   V(IJK-NYZ)
               W(IJK)=   W(IJK-NYZ)
              TE(IJK)=  TE(IJK-NYZ)
              ED(IJK)=  ED(IJK-NYZ)
             VIS(IJK)= VIS(IJK-NYZ)
            VIST(IJK)=VIST(IJK-NYZ)
        END DO
        END DO
      ENDIF
C      
c     IF(NPRO.GT.1)THEN
c       CALL MPI_BCAST(FAC,1,MPI_DOUBLE_PRECISION,IPRN,
c    *                 MPI_COMM_WORLD,IERR)
c     ENDIF
C
C.....CORRECT VELOCITY AND MASS FLUXES AT TOP
C
c       K=NZ
c       DO I=IXBN,IXED
c         II=I-IREF
c       DO J=2,NYM 
c         IJK=LI(II)+LJ(J)+K
c         AR=(X(I)-X(I-1))*(Y(J)-Y(J-1))
c         W(IJK)=-FAC
c         F3(IJK-1)=W(IJK)*DENSIT*AR
c         F3(IJK)  =F3(IJK-1)
c       END DO
c       END DO
C
      RETURN
      END
C
C############################################################
      SUBROUTINE SETBC
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
C
C---------------------------------------------------
C..... SET BOUNDARY CONDITIONS
C---------------------------------------------------
C
      DO I=2,NXM 
      DO K=2,NZM 
C
        IJK=LI(I)+LJ(1)+K
        U(IJK)   =UIN*((ZC(K)/REFH)**EXPON)
        V(IJK)   =VIN
        W(IJK)   =WIN
        TE(IJK)  =TEIN
        ED(IJK)  =EDIN
        VIST(IJK)=CMU*DENSIT*TEIN*TEIN/(EDIN+SMALL)
        VIS(IJK) =VISCOS+VIST(IJK)
C
        IJK=LI(I)+LJ(NY)+K
        U(IJK)   =UIN*((ZC(K)/REFH)**EXPON)
        V(IJK)   =VIN
        W(IJK)   =WIN
        TE(IJK)  =TEIN
        ED(IJK)  =EDIN
        VIST(IJK)=CMU*DENSIT*TEIN*TEIN/(EDIN+SMALL)
        VIS(IJK) =VISCOS+VIST(IJK)
      END DO
      END DO
C
C     WEST & EAST 
C
      IF(MYID.EQ.IROOT)THEN
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(1)+LJ(J)+K
        U(IJK)   =UIN*((ZC(K)/REFH)**EXPON)
        V(IJK)   =VIN
        W(IJK)   =WIN
        TE(IJK)  =TEIN
        ED(IJK)  =EDIN
        VIST(IJK)=CMU*DENSIT*TEIN*TEIN/(EDIN+SMALL)
        VIS(IJK) =VISCOS+VIST(IJK)
      END DO
      END DO
      ENDIF
C
c     IF(MYID.EQ.IPRN)THEN
c     DO J=2,NYM 
c     DO K=2,NZM 
c       IJK=LI(NX)+LJ(J)+K
c       U(IJK)   = U(IJK-NYZ)
c       V(IJK)   = V(IJK-NYZ)
c       W(IJK)   = W(IJK-NYZ)
c       TE(IJK)  =TE(IJK-NYZ)
c       ED(IJK)  =ED(IJK-NYZ)
c     END DO
c     END DO
c     ENDIF
C
C     BOTTOM & TOP
C
      DO J=2,NYM 
      DO I=2,NXM 
C
        IJK=LI(I)+LJ(J)+1
        U(IJK)   =ZERO
        V(IJK)   =ZERO
        W(IJK)   =ZERO
C
        IJK=LI(I)+LJ(J)+NZ
        U(IJK)   =UIN*((ZC(K)/REFH)**EXPON)
        V(IJK)   =VIN
        W(IJK)   =WIN
        TE(IJK)  =TEIN
        ED(IJK)  =EDIN
        VIST(IJK)=CMU*DENSIT*TEIN*TEIN/(EDIN+SMALL)
        VIS(IJK) =VISCOS+VIST(IJK)
C
      END DO
      END DO
C
      RETURN
      END
C
C############################################################
      SUBROUTINE INFLOW
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION AREA,AR,TMP,TT
C
C---------------------------------------------------
C
C.....CALCULATE INFLOW MASS FLUX 
C
      AREA=ZERO
      FLOMAS=ZERO
C
C     I=1
C     write(*,*) 'MYID=',MYID,IROOT,IREF,I-IREF
      IF(MYID.EQ.IROOT)THEN
      I=1
      II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        AR=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
        AREA=AREA+AR
        F1(IJK)=DENSIT*AR*U(IJK)
        FLOMAS=FLOMAS+F1(IJK)
      END DO
      END DO
      ENDIF
C
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(1)+K
        F2(IJK)=ZERO
        IJK=LI(II)+LJ(NY)+K
        F2(IJK)=ZERO
        F2(IJK-NZ)=F2(IJK)
C
      END DO
      END DO
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+1
        F3(IJK)=ZERO
        IJK=LI(II)+LJ(J)+NZ
        AR=(X(I)-X(I-1))*(Y(J)-Y(J-1))
        F3(IJK)=DENSIT*AR*W(IJK)
        F3(IJK-1)=F3(IJK)
      END DO
      END DO
C
      IF(NPRO.GT.1)THEN
        CALL MPI_BCAST(FLOMAS,1,MPI_DOUBLE_PRECISION,IROOT,
     *                 MPI_COMM_WORLD,IERR)
      ENDIF
      TMP=FLOMAS
      TT=UREF
      SNORM(IU)=TMP*UREF
      SNORM(IV)=TMP*UREF
      SNORM(IW)=TMP*UREF
      SNORM(IP)=TMP
      SNORM(ITE)=TMP*TT*TT
      SNORM(IED)=TMP*(TT**(UNITY+HALF))/XLREF
C
      RETURN
      END
C
C############################################################
      SUBROUTINE MODINP
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
      CHARACTER CTMP*512
C
C-----------------------------------------------------------
C.....READ CONTROL DATA FROM UNIT 2
C-----------------------------------------------------------
C
C.....SET SOME CONTROL VARIABLES
C
      IU =1
      IV =2
      IW =3
      IP =4
      ITE=5
      IED=6
C
      IVIS=0
      IVIST=1
C
      IWAL=0
      IINL=1
      IOUT=2
      ISYM=3
C
C     OPEN OUTPUT FILE
C
  116 FORMAT(1X,'(*) Processor Number =',I4)
  111 FORMAT(1X,'(*) Using ADM Module, NWT =',I4)
      IF(MYID.EQ.IPRN)THEN
        OPEN (UNIT=2,FILE=TRIM(FILCHK),IOSTAT=IERR)
        IF(IERR.NE.0) CALL FILEERROR(TRIM(FILCHK))
        WRITE(2,'(80A1)') ('=',I=1,80) 
        WRITE(2,*) 'CONTROL DATA:' 
        WRITE(2,'(80A1)') ('=',I=1,80) 
        WRITE(2,116) NPRO 
        WRITE(2,111) NWT 
      ENDIF
C
C     LHUB   : Hub considered (T/F)
C     LTOWER : Tower considered (T/F)
C     LBLADE : Blade considered (T/F)
C     LRATED : Power Calibration (T/F)
C     EXPCT  : Exponent of CT of the Fitted Curve
C     EXPCP  : Exponent of CP of the Fitted Curve 
C     EXPLA  : Exponent of Lamda of the Fitted Curve
C     CGSCA  : Scaling constant of CG of the Fitted Curve
C     NWT    : Number of wind turbines
C     NWTD   : Number of wind turbine type
C     WBX    : X coordinate of the tower base center (m)
C     WBY    : Y coordinate of the tower base center (m)
C     WBZ    : Z coordinate of the tower base center (m)
C     ROTD   : Rotor diameter (m)
C     HUBH   : Hub height (m)
C     TTANG  : Tilt angle (deg)
C     VRATED : Rated wind speed (m/s)
C     WRATED : Rated Rotor speed (RPM)
C     PRATED : Rated Power (W)
C
   11 FORMAT(A256)
      OPEN (UNIT=5,FILE=TRIM(FILWTD),IOSTAT=IERR)
      IF(IERR.NE.0) CALL FILEERROR(TRIM(FILWTD))
      READ(5,11) CTMP 
      READ(5,*)  LHUB,LTOWER,LBLADE,LRATED
      READ(5,11) CTMP 
      READ(5,*)  NWT,NWTD
      READ(5,11) CTMP 
      READ(5,*)  (WBX(I),WBY(I),WBZ(I),IWTY(I),I=1,NWT)
C     write(*,*)  (I,WBX(I),WBY(I),WBZ(I),IWTY(I),I=1,NWT)
C
      DO I=1,NWTD
      READ(5,11) CTMP 
      READ(5,*)  J
      READ(5,11) CTMP 
      READ(5,*)  EXPCT(J),EXPCP(J),EXPLA(J),CGSCA
      READ(5,11) CTMP 
      READ(5,*)  ROTD(J),HUBH(J),TTANG(J),VRATED(J),PRATED(J),WRATED(J)
     *           ,TERATED(J)
      READ(5,11) CTMP 
      READ(5,*)  HUBD(J),HCD(J)
      READ(5,11) CTMP 
      READ(5,*)  NTR(J),TCD(J),XTC(J)
      READ(5,11) CTMP 
      READ(5,*) (HTW(J,K),RTW(J,K),K=1,NTR(J))
C
      HUBR(J)=HUBD(J)*HALF
      TTANG(J)=TTANG(J)*DEG2RAD
C
C     NAZ   : Number of azimuthal divisions
C     NRA   : Number of radial divisions
C     NPC   : Number of points describing power curve
C     PCS   : The average wind speed arcoss rotor (m/s) 
C     WDS   : The wind speed (m/s)
C     PCR   : The rotor speed (rpm) 
C     PCP   : The power (MW)
C     HUBR  : Hub radius (m)
C     HCD   : Drag coefficient of disk
C     TCD   : Drag coefficient of cylinder
C     NTR   : Number of points describing tower radius curve 
C     XTC   : Tower center offset
C     NTP   : Number of tower nodes
C     HTW   : The height offset of the Ith point on radius curve (m)
C     RTW   : The radius offset of the Ith point on radius curve (m)
C     TTE   : Total turbulent kinetic energu production at rotor
C
      READ(5,11) CTMP 
      READ(5,*)  NPC(J)
      READ(5,11) CTMP 
      READ(5,*) (WDS(J,K),PCS(J,K),PCP(J,K),PCR(J,K),TTE(J,K),
     *           K=1,NPC(J))
      END DO
C
      CLOSE(5)
C
C     OPEN CONTROL FILE
C
      OPEN (UNIT=5,FILE=TRIM(FILCTL),IOSTAT=IERR)
      IF(IERR.NE.0) CALL FILEERROR(TRIM(FILCTL))
C
      READ(5,11) CTMP 
      READ(5,6) TITLE
    6 FORMAT(A80)
C
C     LREAD:  read results from previous run
C     LWRITE: write results onto a file
C     LTEST:  print additional information 
C     LOUTS:  print convergence history on the screen          
C     LOUTE:  print convergence history onto a file          
C
      READ(5,11) CTMP 
      READ(5,*) XMF,YMF,ZMF
C
      READ(5,11) CTMP 
      READ(5,*) LREAD,LWRITE,IFRE,LADM,BFINC,URFUVW
      FILRES=TRIM(FWDIR)//'/'//TRIM(FCASE)
      FILRST=TRIM(FRDIR)//'/'//TRIM(FCASE)
C
C     MAXIT:  the maximum number of outer iterations
C     IMON:   the I index of the monitoring location
C     JMON:   the J index of the monitoring location
C     KMON:   the K index of the monitoring location
C     IPR:    the I index of the pressure reference point (p=0)
C     JPR:    the J index of the pressure reference point (p=0)
C     KPR:    the K index of the pressure reference point (p=0)
C     SORMAX: the level of  residual norm 
C     SLARGE: the large value 
C     ALFA:   the parameter in the SIP solver
C
      READ(5,11) CTMP 
      READ(5,*) MAXIT,SORMAX
      SLARGE=1D30
      ALFA=0.92D0
      IPR=2
      JPR=2
      KPR=2
      IMON=2
      JMON=2
      KMON=2
C
C     DENSIT: the initial fluid density (kg/m3)
C     VISCOS: the initial viscosity (Pa-s) 
C     GRAVX:  the X component of the gravity vector (m/s2)
C     GRAVY:  the Y component of the gravity vector (m/s2)
C     GRAVZ:  the Z component of the gravity vector (m/s2)
C
      READ(5,11) CTMP 
      READ(5,*) DENSIT,VISCOS
      GRAVX=ZERO
      GRAVY=ZERO
      GRAVZ=ZERO
C
C     UIN, VIN, WIN, PIN, TEIN, EDIN are the values of 
C     U, V, W, P, TE, ED employed to initialize fields 
C     (usually zero, or some mean values).
C     EXPON:   the log-law expont of wind profile 
C     REFH:    reference height for UIN
C     REFD:    reference diameter for TKE 
C
      READ(5,11) CTMP 
      READ(5,*) UIN,TEIN,XLREF,EXPON,GDS(IU),REFH,REFD
      PIN=ZERO
      VIN=ZERO
      WIN=ZERO
      IF(LRATED) UIN=VRATED(1)
C
C     IF(MYID.EQ.IPRN) WRITE(2,*) UIN,VIN,WIN,PIN,TEIN,EDIN,EXPON
C
C     XLREF:   the reference length of turbulence (m)
C     IFRE:    the output frequency of binary result
C
      UREF=SQRT(UIN*UIN+VIN*VIN+WIN*WIN)
      TEIN=TEIN*UREF*UREF
      EDIN=TEIN**1.5/XLREF
C
C     LCAL(I) defines which equations are to be solved 
C     (I defines variable as follows: 
C     1->U,2->V,3->W,4->P,5->TE,6->ED
C
      DO I=1,NPHI
        LCAL(I)=.TRUE.
      END DO
C
C     URF(I) is the under-relaxation factor for the Ith variable.
C
      READ(5,11) CTMP 
      READ(5,*) URF(1),URF(4),URF(5)
      URF(6)=URF(5)
      URF(2)=URF(1)
      URF(3)=URF(1)
C
C     SOR(I) is the required ratio of reduction of the residual norm
C     during inner iterations for Ith variable before they are stoped 
C
      READ(5,11) CTMP 
      READ(5,*) SOR(1),SOR(4),SOR(5)
      SOR(6)=SOR(5)
      SOR(2)=SOR(1)
      SOR(3)=SOR(1)
C
C     NSW(I) is the maximum allowed number of inner iterations for the
C     Ith variable.
C
      READ(5,11) CTMP 
      READ(5,*) NSW(1),NSW(4),NSW(5)
      NSW(6)=NSW(5)
      NSW(2)=NSW(1)
      NSW(3)=NSW(1)
C
C     GDS(I) is the blending factor for UDS and CDS in the equation for
C     the Ith variable.
C
      DO I=2,NPHI
        GDS(I)=ZERO
      END DO
      GDS(IV)=GDS(IU)
      GDS(IW)=GDS(IU)
C
C     LSORM(I) is for checking the residuum of the Ith variable.
C
      DO I=1,NPHI
        LSORM(I)=.TRUE.
      END DO
C
C     URFM(I) is the under-relaxation factor for the Ith material property.
C     IVIST=1
C
C     Constants for k-e turbulence model                                       
C
      READ(5,11) CTMP 
      READ(5,*) CMU,CAPPA,STE,SED,CED1,CED2,ELOG,CTRANS,SFVIST,URFM(1)
      SFVIST=MAX(1.0,SFVIST)
      READ(5,11) CTMP
      READ(5,*) LPRN
C
C     Constants for earth rotation and latitude angle
C     
C     EROT in rad per second
C     ALAT in deg      
C
      READ(5,11) CTMP
      READ(5,*) EROT,ALAT
      ALAT=ALAT*DEG2RAD
      CLOSE(5)
C
      LUVW=LCAL(IU).OR.LCAL(IV).OR.LCAL(IW)
      LPRE=LCAL(IP)
      LVIST=LCAL(ITE).OR.LCAL(IED)
C
      DO I=1,NPHI
        RESOR(I)=ZERO
      END DO
      CMU75=CMU**THREEQR
      CMU25=CMU**QUARTER
C
C.....READ GRID DATA 
C
      OPEN (UNIT=1,FILE=TRIM(FILGRD),IOSTAT=IERR)
      IF(IERR.NE.0) CALL FILEERROR(TRIM(FILGRD))
C
      READ(1,*) NI
      READ(1,*) NJ
      READ(1,*) NK
C
C.....CHECK DIMENSIONS
C
      IF(MYID.EQ.IPRN)THEN
        I=0
        IF(NI.NE.NXT)THEN
          WRITE(*,*) '(E) NI != ',NXT
          I=I+1
        ENDIF
        IF(NJ.NE.NY)THEN
          WRITE(*,*) '(E) NJ != ',NY
          I=I+1
        ENDIF
        IF(NK.NE.NZ)THEN 
          WRITE(*,*) '(E) NK != ',NZ
          I=I+1
        ENDIF
        IF(I.NE.0)THEN
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
      ENDIF
C
      READ(1,*) (X(I),I=1,NXT)
      READ(1,*) (Y(J),J=1,NY)
      READ(1,*) (Z(K),K=1,NZ)
C
      CLOSE(1)
C
      END
C
C############################################################
      SUBROUTINE INITFIELD
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
C
      DO I=1,NX
        LI(I)=(I-1)*NYZ
      END DO
      DO J=1,NY
        LJ(J)=(J-1)*NZ
      END DO
C
      DO I=1,NXT
       X(I)=X(I)*XMF 
      END DO
      DO J=1,NY                             
       Y(J)=Y(J)*YMF 
      END DO
      DO K=1,NZ                             
       Z(K)=Z(K)*ZMF 
      END DO
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
C.....INTERPOLATION FACTORS
C
      DO I=1,NXTM 
        FX(I)=(X(I)-XC(I))/(XC(I+1)-XC(I))
      END DO
      FX(NXT)=UNITY
C
      DO J=1,NYM 
        FY(J)=(Y(J)-YC(J))/(YC(J+1)-YC(J))
      END DO
      FY(NY)=UNITY
C
      DO K=1,NZM 
        FZ(K)=(Z(K)-ZC(K))/(ZC(K+1)-ZC(K))
      END DO
      FZ(NZ)=UNITY
C
C......SET BOUNDARY CONDITION INDEX
C
      DO J=1,NY
      DO K=1,NZ
        JK=(J-1)*NZ+K
        IBCW(JK)=IINL
        IBCE(JK)=IOUT 
      END DO
      END DO
C
      DO I=1,NX
      DO K=1,NZ
        IK=(I-1)*NZ+K
        IBCS(IK)=ISYM
        IBCN(IK)=ISYM
      END DO
      END DO
C
      DO I=1,NX
      DO J=1,NY
        IJ=(I-1)*NY+J
        IBCB(IJ)=IWAL
        IBCT(IJ)=IINL
      END DO
      END DO
C
C.....INITIALIZE VARIBLE VALUES
c
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
          IJK=LI(I)+LJ(J)+K
          U(IJK)=UIN*((ZC(K)/REFH)**EXPON)
          V(IJK)=VIN 
          W(IJK)=WIN 
          P(IJK)=ZERO
          TE(IJK)=TEIN
          ED(IJK)=EDIN
          VOL(IJK)=ZERO
          F1(IJK)=ZERO
          F2(IJK)=ZERO
          F3(IJK)=ZERO
          VIS(IJK)=VISCOS
          VIST(IJK)=ZERO
      END DO
      END DO
      END DO
C
      DO I=1,NYZ
         YPLW(I)=ZERO
         YPLE(I)=ZERO
         VSWW(I)=VISCOS
         VSWE(I)=VISCOS
      END DO
      DO J=1,NXZ
         YPLS(J)=ZERO
         YPLN(J)=ZERO
         VSWS(J)=VISCOS
         VSWN(J)=VISCOS
      END DO
      DO K=1,NXY
         YPLB(K)=ZERO
         YPLT(K)=ZERO
         VSWB(K)=VISCOS
         VSWT(K)=VISCOS
      END DO
C
C.....CALCULATE VOLUME
C
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          VOL(IJK)=(X(I)-X(I-1))*(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
      END DO
      END DO
      END DO
C
      IF(LREAD)THEN
       BFP=UNITY
      ELSE
       BFP=ZERO
      ENDIF
C
      ARTOP=ZERO
      DO I=2,NXTM 
      DO J=2,NYM 
        ARTOP=ARTOP+(X(I)-X(I-1))*(Y(J)-Y(J-1))
      END DO
      END DO
C
C------------------------------------------------------
C.....INITIAL OUTPUT - PRINTOUT OF FLOW PARAMETERS
C------------------------------------------------------
C
      IF((MYID.EQ.IPRN).AND.LWRITE)THEN
      WRITE(2,'(80A1)') ('=',I=1,80) 
      WRITE(2,*) 'SUMMARY OF SYSTEM PARAMETERS:',TITLE
      WRITE(2,'(80A1)') ('=',I=1,80) 
      WRITE(2,599) LREAD,LWRITE
  599 FORMAT(1X,//,10X,50('*'),/,10X,
     *        ' LREAD :  ',L1,/,10X,
     *        ' LWRITE:  ',L1,/,10X)
      WRITE(2,600) UIN,VIN,WIN,TEIN,EDIN
  600 FORMAT(1X,//,10X,50('*'),/,10X,
     *        ' UIN   :  ',1P1E12.4,/,10X,
     *        ' VIN   :  ',1P1E12.4,/,10X,
     *        ' WIN   :  ',1P1E12.4,/,10X,
     *        ' TEIN  :  ',1P1E12.4,/,10X,
     *        ' EDIN  :  ',1P1E12.4,/,10X)
      WRITE(2,601) DENSIT,VISCOS
  601 FORMAT(1X,//,10X,50('*'),/,10X,
     *        ' FLUID DENSITY    :  ',1P1E10.4,/,10X,
     *        ' DYNAMIC VISCOSITY:  ',1P1E10.4,/,10X)
   11 FORMAT(A256)
    7 FORMAT(1X,A25,L1,1X,E10.4,1X,I2)
    8 FORMAT(1X,A25,I6)
    9 FORMAT(1X,A25,E10.4,1X,L1,1X,I2)
      WRITE(2,*) '  '
      WRITE(2,*) '          ALFA  PARAMETER  :  ',ALFA
      WRITE(2,*) '  '
      WRITE(2,*) '          UNDERRELAXATION  FACTORS, ACTIVE AND ID' 
      WRITE(2,*) '          ==========================================='
      WRITE(2,9) '          U-VELOCITY  :  ',URF(IU),LCAL(IU),IU
      WRITE(2,9) '          V-VELOCITY  :  ',URF(IV),LCAL(IV),IV
      WRITE(2,9) '          W-VELOCITY  :  ',URF(IW),LCAL(IW),IW
      WRITE(2,9) '          PRESSURE    :  ',URF(IP),LCAL(IP),IP
      WRITE(2,9) '          TUR. K.E.   :  ',URF(ITE),LCAL(ITE),ITE
      WRITE(2,9) '          DISSP. RATE :  ',URF(IED),LCAL(IED),IED
      WRITE(2,*) '  '
      WRITE(2,*) '          SPATIAL BLENDING FACTORS (CDS-UDS)'
      WRITE(2,*) '          =================================='
      WRITE(2,9) '          U-VELOCITY  :  ',GDS(IU)
      WRITE(2,9) '          V-VELOCITY  :  ',GDS(IV)
      WRITE(2,9) '          W-VELOCITY  :  ',GDS(IW)
      WRITE(2,9) '          PRESSURE    :  ',GDS(IP)
      WRITE(2,9) '          TUR. K.E.   :  ',GDS(ITE)
      WRITE(2,9) '          DISSP. RATE :  ',GDS(IED)
      WRITE(2,*) '  '
      WRITE(2,*) '          SOR'
      WRITE(2,*) '          =================================='
      WRITE(2,9) '          U-VELOCITY  :  ',SOR(IU)
      WRITE(2,9) '          V-VELOCITY  :  ',SOR(IV)
      WRITE(2,9) '          W-VELOCITY  :  ',SOR(IW)
      WRITE(2,9) '          PRESSURE    :  ',SOR(IP)
      WRITE(2,9) '          TUR. K.E.   :  ',SOR(ITE)
      WRITE(2,9) '          DISSP. RATE :  ',SOR(IED)
      WRITE(2,*) '  '
      WRITE(2,*) '          NSW'
      WRITE(2,*) '          =================================='
      WRITE(2,8) '          U-VELOCITY  :  ',NSW(IU)
      WRITE(2,8) '          V-VELOCITY  :  ',NSW(IV)
      WRITE(2,8) '          W-VELOCITY  :  ',NSW(IW)
      WRITE(2,8) '          PRESSURE    :  ',NSW(IP)
      WRITE(2,8) '          TUR. K.E.   :  ',NSW(ITE)
      WRITE(2,8) '          DISSP. RATE :  ',NSW(IED)
      WRITE(2,*) '  '
      WRITE(2,*) '          LSORM'
      WRITE(2,*) '          =================================='
      WRITE(2,7) '          U-VELOCITY  :  ',LSORM(IU)
      WRITE(2,7) '          V-VELOCITY  :  ',LSORM(IV)
      WRITE(2,7) '          W-VELOCITY  :  ',LSORM(IW)
      WRITE(2,7) '          PRESSURE    :  ',LSORM(IP)
      WRITE(2,7) '          TUR. K.E.   :  ',LSORM(ITE)
      WRITE(2,7) '          DISSP. RATE :  ',LSORM(IED)
      WRITE(2,*) '  '
      WRITE(2,*) '          FIELD PROPERTY UPDATE: (T/F), URFM, ID'
      WRITE(2,*) '          ==========================================='
      WRITE(2,7) '          TURB. VISCOS:  ',LVIST,URFM(IVIST),IVIST
      WRITE(2,7) '  '
      WRITE(2,603) XMF,YMF,ZMF
  603 FORMAT(1X,/,10x,
     *        ' XMF    :  ',1P1E10.4,/,10X,
     *        ' YMF    :  ',1P1E10.4,/,10X,
     *        ' ZMF    :  ',1P1E10.4)
      WRITE(2,604) IRCFP
  604 FORMAT(1X,/,10x,
     *        ' IRCFP  :  ',I10,/,10X)
      WRITE(2,605) CMU,CMU25,CMU75,CAPPA,STE,SED,CED1,CED2,ELOG,CTRANS
  605 FORMAT(1X,/,10x,
     *        ' CMU    :  ',1P1E10.4,/,10X,
     *        ' CMU25  :  ',1P1E10.4,/,10X,
     *        ' CMUR75 :  ',1P1E10.4,/,10X,
     *        ' CAPPA  :  ',1P1E10.4,/,10X,
     *        ' STE    :  ',1P1E10.4,/,10X,
     *        ' SED    :  ',1P1E10.4,/,10X,
     *        ' CED1   :  ',1P1E10.4,/,10X,
     *        ' CED2   :  ',1P1E10.4,/,10X,
     *        ' ELOG   :  ',1P1E10.4,/,10X,
     *        ' CTRANS :  ',1P1E10.4)
      WRITE(2,607) MAXIT,SORMAX,SLARGE,ALFA,IFRE
  607 FORMAT(1X,/,10x,
     *        ' MAXIT  :  ',I10,/,10X,
     *        ' SORMAX :  ',E10.4,/,10X,
     *        ' SLARGE :  ',E10.4,/,10X,
     *        ' ALFA   :  ',E10.4,/,10X,
     *        ' IFRE   :  ',I10)
      WRITE(2,608) UREF,XLREF
  608 FORMAT(1X,/,10x,
     *        ' UREF   :  ',E10.4,/,10X,
     *        ' XLREF  :  ',E10.4)
      ENDIF
C
      CLOSE(2)
      RETURN
      END
C
C##########################################################
      SUBROUTINE INITFP                                       
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
C
      DO IJK=1,NXYZ
        VIS(IJK)=VISCOS
        VIST(IJK)=ZERO
        SHR(IJK)=ZERO
        DPX(IJK)=ZERO
        DPY(IJK)=ZERO
        DPZ(IJK)=ZERO
        DUX(IJK)=ZERO
        DUY(IJK)=ZERO
        DUZ(IJK)=ZERO
        DVX(IJK)=ZERO
        DVY(IJK)=ZERO
        DVZ(IJK)=ZERO
        DWX(IJK)=ZERO
        DWY(IJK)=ZERO
        DWZ(IJK)=ZERO
      END DO 
C
      RETURN
      END
C
C##########################################################
      SUBROUTINE FIND(XYZ,VALMAX,VALMIN)
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION XYZ(NXYZ),VALMAX,VALMIN,TMPA,TMPB,TMP
      DOUBLE PRECISION TTA(1),TTB(1),TT(1)
C
      TMPA=-GREAT
      TMPB= GREAT
      DO I=2,NXM 
      DO J=2,NYM 
      DO K=2,NZM 
         IJK=LI(I)+LJ(J)+K
         IF(XYZ(IJK).GT.TMPA)THEN
           TMPA=XYZ(IJK)
         ENDIF
         IF(XYZ(IJK).LT.TMPB)THEN
           TMPB=XYZ(IJK)
         ENDIF
      END DO
      END DO
      END DO
      IF(NPRO.GT.1)THEN
        TTA(1)=TMPA
        TTB(1)=TMPB
        CALL MPI_ALLREDUCE(TTA,TT,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     *                     MPI_COMM_WORLD,IERR)
        CALL MPI_ALLREDUCE(TTB,TT,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     *                     MPI_COMM_WORLD,IERR)
        TMPA=TT(1)
        TMPB=TT(1) 
      ENDIF
      VALMAX=TMPA
      VALMIN=TMPB
C
      RETURN
      END
C
C##########################################################
      SUBROUTINE APPEND_INDEX(ID,FNAME)
C##########################################################
      CHARACTER*128 FNAME
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
      END
C
C##########################################################
      SUBROUTINE INITMPI(FNRST,FNRES,IPRO) 
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'mpif.h'
      CHARACTER*128 FNRST,FNRES
C     DIMENSION ITMP(2*IPRO),IX(2)
C
      IPROM=IPRO-1
C
C     IPRN: OUTPUT PROCESSOR
C
      ILFT=MYID-1
      IRHT=MYID+1
C     
      IXBN=MYID*(NXT-2)/IPRO+2
      IXED=(MYID+1)*(NXT-2)/IPRO+1
      NX=IXED-IXBN+3
      NXM=NX-1
      NXMM=NXM-1
C
      IF(NX.GT.MX)THEN
        CALL MPI_FINALIZE(IERR)
        IF(MYID.EQ.IPRN) WRITE(*,*) "(E) Array Dimension."
        STOP
      ENDIF
C
      IFBN=IXBN-1
      IF(IFBN.EQ.1) IFBN=2
      IFED=IXED
      IF(IFED.EQ.NXTM) IFED=NXT-2
C
      IREF=IXBN-2
C
      DO I=1,IPRO
        IZBN(I)=(I-1)*(NXT-2)/IPRO+2
        IZED(I)=(I  )*(NXT-2)/IPRO+1
      END DO
C
C     First: SEND DATA TO LEFT AND RECEIVE DATA FROM RIGHT
C     
      IF(ILFT.LT.0)THEN
        INS(1)=MPI_PROC_NULL
      ELSE
        INS(1)=ILFT
      ENDIF
      IF(IRHT.GT.IPROM)THEN
        INR(1)=MPI_PROC_NULL
      ELSE
        INR(1)=IRHT
      ENDIF
C        
      ISR(1)=2
      IRR(1)=NX
C
C     Second: SEND DATA TO RIGHT AND RECEIVE DATA FROM LEFT 
C     
      IF(IRHT.GT.IPROM)THEN
        INS(2)=MPI_PROC_NULL
      ELSE
        INS(2)=IRHT
      ENDIF
      IF(ILFT.LT.0)THEN
        INR(2)=MPI_PROC_NULL
      ELSE
        INR(2)=ILFT
      ENDIF
C     
      ISR(2)=NXM 
      IRR(2)=1
C
      CALL APPEND_INDEX(MYID+1,FNRST)
      CALL APPEND_INDEX(MYID+1,FNRES)
C
      NP=NYMM*NZMM
      CALL MPI_PACK_SIZE(NP,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ISIZE,
     *                   IERR)
C
      IF(IPR.LT.2) IPR=2
      IF(IPR.GE.NXT) IPR=NXTM
      IF(JPR.LT.2) JPR=2
      IF(JPR.GE.NY) JPR=NYM 
      IF(KPR.LT.2) KPR=2
      IF(KPR.GE.NZ) KPR=NZM 
C
      IPPR=0
      DO I=1,IPRO
        IF((IPR.GE.IZBN(I)).AND.(IPR.LE.IZED(I))) IPPR=I-1
      ENDDO  
C
      END
C
C##########################################################
      SUBROUTINE EXCHFI(NSIZE,FI,ID)
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION BUFS(NYZ),FI(NXYZ)
      INTEGER ISTATUS(MPI_STATUS_SIZE),ID
      CHARACTER BCS(NSIZE),BCR(NSIZE)
C
      NP=NYMM*NZMM
      DO N=1,2
C
        I=ISR(N)
        JK=0
C
        JSTMP=NYM 
        JRTMP=NYM 
        KSTMP=NZM 
        KRTMP=NZM 
        IF(INS(N).NE.MPI_PROC_NULL)THEN
          DO J=2,JSTMP
          DO K=2,KSTMP
            IJK=LI(I)+LJ(J)+K
            JK=JK+1
            BUFS(JK)=FI(IJK)
          ENDDO
          ENDDO
        ELSE
          DO J=2,JSTMP
          DO K=2,KSTMP
            IJK=LI(I)+LJ(J)+K
            JK=JK+1
            BUFS(JK)=ZERO
          ENDDO
          ENDDO
        ENDIF
C
        IPOS=0
        CALL MPI_PACK(BUFS,NP,MPI_DOUBLE_PRECISION,BCS,NSIZE,IPOS,
     &                MPI_COMM_WORLD,IERR)
        ITAG=ID+N
        CALL MPI_SENDRECV(BCS,NSIZE,MPI_CHARACTER,INS(N),ITAG,
     &                    BCR,NSIZE,MPI_CHARACTER,INR(N),ITAG,
     &                    MPI_COMM_WORLD,ISTATUS,IERR)
        IPOS=0
        CALL MPI_UNPACK(BCR,NSIZE,IPOS,BUFS,NP,MPI_DOUBLE_PRECISION,
     &                  MPI_COMM_WORLD,IERR)
C
        I=IRR(N)
        JK=0
        IF(INR(N).NE.MPI_PROC_NULL)THEN
          DO J=2,JRTMP
          DO K=2,KRTMP
            IJK=LI(I)+LJ(J)+K
            JK=JK+1
            FI(IJK)=BUFS(JK)
          ENDDO
          ENDDO
        ELSE
          DO J=2,JRTMP
          DO K=2,KRTMP
            IJK=LI(I)+LJ(J)+K
            JK=JK+1
            BUFS(JK)=ZERO
          ENDDO
          ENDDO
        ENDIF
C
      ENDDO
C
      END
C
C##########################################################
      SUBROUTINE LOADCAL
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION RT(NPRO),TMP(1),ALL(1)
C
      TMP(1)=NXMM*NYMM*NZMM
      CALL MPI_ALLREDUCE(TMP,ALL,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     *                       MPI_COMM_WORLD,IERR)
      TMP(1)=TMP(1)/ALL(1)
C
      CALL MPI_ALLGATHER(TMP,1,MPI_DOUBLE_PRECISION,RT,1,
     *                   MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
      IF(MYID.EQ.IPRN)THEN
        WRITE(2,*) ' LOADING of CPUS:' 
        WRITE(2,*) ' =================================='
        DO I=1,NPRO
          WRITE(2,*) ' CPU ID=',I,' LOADING=',RT(I)
        END DO
        WRITE(2,*) ' ' 
      ENDIF
      END
C
C############################################################
      SUBROUTINE CALCKE(ID,WW)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'wind.inc'
      DOUBLE PRECISION WW(NXYZ),TMPA,START,FINISH
      DOUBLE PRECISION URFI,SF,WWIN,S,VISC,VISD,D,CP,FUDS,FCDS,TMP
      DOUBLE PRECISION FXE,FXP,DXPE,CE
      DOUBLE PRECISION FYN,FYP,DYPN,CN
      DOUBLE PRECISION FZT,FZP,DZPT,CT
C
      IF(LADM.AND.BFP.LT.UNITY)THEN
        BFP =BFP+BFINC
        URFP=UNITY*URFUVW
      ELSE
        BFP =UNITY
        URFP=UNITY
      END IF
C
C.....INITIALIZATION OF TEMPORARILY STORED VARIABLES
C
      DO IJK=1,NXYZ
        SU(IJK)=ZERO
        AP(IJK)=ZERO
        AE(IJK)=ZERO
        AW(IJK)=ZERO
        AN(IJK)=ZERO
        AS(IJK)=ZERO
        AT(IJK)=ZERO
        AB(IJK)=ZERO
      END DO
C      
      IF(ID.EQ.ITE)THEN
        SF=UNITY/STE
        WWIN=TEIN
        URFI=UNITY/(URF(ID)*URFP)
      ELSE IF(ID.EQ.IED)THEN
        SF=UNITY/SED
        WWIN=EDIN
        URFI=UNITY/URF(ID)
      ELSE 
        IF(MYID.EQ.IPRN) WRITE(2,*) 'CALCKE: ID Error!'
        STOP
      ENDIF
C
C==========================================================
C.....FLUXES THROUGH INTERNAL EAST CV-FACES 
C==========================================================
C
      DO I=IFBN,IFED
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO E
C
        II  =I-IREF
        FXE =FX(I)
        FXP =UNITY-FXE
        DXPE=XC(I+1)-XC(I)
C
        DO J=2,NYM 
        DO K=2,NZM 
C
          IJK=LI(II)+LJ(J)+K
          IJE=IJK+NYZ
C
C.....CELL FACE AREA S = DY*RE*1
C
          S=(Y(J)-Y(J-1))*(Z(K)-Z(K-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VISCOS
          VISD=VIST(IJK)*FXP+VIST(IJE)*FXE
          VISC=VISC+VISD*SF
          D=VISC*S/DXPE
C
C.....EXPLICIT CONVECTIVE FLUX FOR UDS AND CDS
C
          CE=MIN(F1(IJK),ZERO)
          CP=MAX(F1(IJK),ZERO)
C
          FUDS=CP*WW(IJK)+CE*WW(IJE)
          FCDS=F1(IJK)*(WW(IJE)*FXE+WW(IJK)*FXP)
C
C.....COEFFICIENTS AE(P) AND AW(E) DUE TO UDS
C
          AE(IJK)= CE-D
          AW(IJE)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND E DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(ID)*(FUDS-FCDS)
          SU(IJE)=SU(IJE)-GDS(ID)*(FUDS-FCDS)
C
        END DO
        END DO
      END DO
C
C=========================================================
C.....FLUXES THROUGH INTERNAL NORTH CV FACES 
C=========================================================
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO N
C
      DO I=IXBN,IXED
C
        II=I-IREF
        DO J=2,NYMM
C
          FYN =FY(J)
          FYP =UNITY-FYN
          DYPN=YC(J+1)-YC(J)
C
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          IJN=IJK+NZ
C
C.....CELL FACE AREA S = DX*RN*1
C
          S=(X(I)-X(I-1))*(Z(K)-Z(K-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VISCOS
          VISD=VIST(IJK)*FYP+VIST(IJN)*FYN
          VISC=VISC+VISD*SF
          D=VISC*S/DYPN
C
C.....EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
C
          CN=MIN(F2(IJK),ZERO)
          CP=MAX(F2(IJK),ZERO)
C
          FUDS=CP*WW(IJK)+CN*WW(IJN)
          FCDS=F2(IJK)*(WW(IJN)*FYN+WW(IJK)*FYP)
C
C.....COEFFICIENTS AN(P) AND AS(N) DUE TO UDS
C
          AN(IJK)= CN-D
          AS(IJN)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND N DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(ID)*(FUDS-FCDS)
          SU(IJN)=SU(IJN)-GDS(ID)*(FUDS-FCDS)
C
       END DO
       END DO
      END DO
C
C=========================================================
C.....FLUXES THROUGH INTERNAL TOP CV FACES 
C=========================================================
C
C.....INTERPOLATION FACTORS, DISTANCE FROM P TO T
C
      DO I=IXBN,IXED
C
        II=I-IREF
        DO J=2,NYM 
        DO K=2,NZMM
C
          FZT =FZ(K)
          FZP =UNITY-FZT
          DZPT=ZC(K+1)-ZC(K)
C
          IJK=LI(II)+LJ(J)+K
          IJT=IJK+1
C
C.....CELL FACE AREA S = DX*RN*1
C
          S=(X(I)-X(I-1))*(Y(J)-Y(J-1))
C
C.....COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
C
          VISC=VISCOS
          VISD=VIST(IJK)*FZP+VIST(IJT)*FZT
          VISC=VISC+VISD*SF
          D=VISC*S/DZPT
C
C.....EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
C
          CT=MIN(F3(IJK),ZERO)
          CP=MAX(F3(IJK),ZERO)
C
          FUDS=CP*WW(IJK)+CT*WW(IJT)
          FCDS=F3(IJK)*(WW(IJT)*FZT+WW(IJK)*FZP)
C
C.....COEFFICIENTS AT(P) AND AB(T) DUE TO UDS
C
          AT(IJK)= CT-D
          AB(IJT)=-CP-D
C
C.....SOURCE TERM CONTRIBUTIONS AT P AND T DUE TO DEFERRED CORRECTION
C
          SU(IJK)=SU(IJK)+GDS(ID)*(FUDS-FCDS)
          SU(IJT)=SU(IJT)-GDS(ID)*(FUDS-FCDS)
C
       END DO
       END DO
      END DO
C
C=============================================================
C.....VOLUME INTEGRALS (SOURCE TERMS)
C=============================================================
      IF(ID.EQ.ITE)THEN
       DO I=IXBN,IXED
         II=I-IREF
       DO J=2,NYM 
       DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          GEN(IJK)=VIST(IJK)*SHR(IJK)
          SU(IJK)=SU(IJK)+GEN(IJK)*VOL(IJK)
          TMPA=DENSIT*VOL(IJK)*ED(IJK)/(TE(IJK)+SMALL)
          AP(IJK)=AP(IJK)+TMPA
       END DO
       END DO
       END DO
C
       IF(LADM)THEN
        DO I=IXBN,IXED
          II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          SU(IJK)=SU(IJK)+WTE(IJK)*BFP
        END DO
        END DO
        END DO
       END IF
C
      ELSE IF(ID.EQ.IED) THEN
       DO I=IXBN,IXED
         II=I-IREF
       DO J=2,NYM 
       DO K=2,NZM 
            IJK=LI(II)+LJ(J)+K
            TMP=ED(IJK)*VOL(IJK)/(TE(IJK)+SMALL)
            SU(IJK)=SU(IJK)+CED1*TMP*GEN(IJK)*DENSIT
            AP(IJK)=AP(IJK)+CED2*TMP*DENSIT
       END DO
       END DO
       END DO
      END IF
C
C=============================================================
C.....PROBLEM MODIFICATIONS - BOUNDARY CONDITIONS
C=============================================================
C
      IF(ID.EQ.ITE)THEN
        CALL BCTE
      ELSE          
        CALL BCED
      ENDIF
C
C==============================================================
C.....UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR TEMPERATURE
C==============================================================
C
      IF(ID.EQ.ITE)THEN
        TMP=UNITY-URF(ID)*URFP
      ELSE 
        TMP=UNITY-URF(ID)
      ENDIF
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          AP(IJK)=(AP(IJK)-AW(IJK)-AE(IJK)-AS(IJK)-AN(IJK)
     &                    -AB(IJK)-AT(IJK))*URFI
          SU(IJK)=SU(IJK)+TMP*AP(IJK)*WW(IJK)
      END DO
      END DO
      END DO
C
      CALL SIPSOL(WW,ID)
C     
      IF(ID.EQ.ITE)THEN
        DO I=IXBN,IXED
          II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
          IJK=LI(II)+LJ(J)+K
          WW(IJK)=MAX(ZERO,WW(IJK))
        END DO
        END DO
        END DO
      ENDIF
C     
      RETURN
      END
C     
C#############################################################
      SUBROUTINE EDWALL(IJK,DN)
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      DOUBLE PRECISION DN,TMP,TEXP
C     
      TMP=GREAT
      TEXP=UNITY+HALF
      ED(IJK)=CMU75*(MAX(ZERO,TE(IJK))**TEXP)/(CAPPA*DN)
      SU(IJK)=ED(IJK)*TMP
      AP(IJK)=TMP
      AS(IJK)=ZERO
      AN(IJK)=ZERO
      AW(IJK)=ZERO
      AE(IJK)=ZERO
      AT(IJK)=ZERO
      AB(IJK)=ZERO
C     
      RETURN
      END
C     
C#############################################################
      SUBROUTINE TEWALL(IJK,IWD,JWD,KWD,IJD,IXYZ,YPLWD,VISWD)
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      DOUBLE PRECISION YPLWD,VISWD
      DOUBLE PRECISION VISS,DN,TMP,TMPU,TMPV,TMPW,TAU
C
      IF(IJD.GT.0)THEN
        IJDN= 1  
      ELSE    
        IJDN=-1 
      ENDIF   
      IJP=IJK
      IJB=IJK+IJD
      TE(IJB)=ZERO
      SU(IJP)=SU(IJP)-GEN(IJP)*VOL(IJP)*DENSIT
      IF(LCAL(ITE).AND.(YPLWD.GT.CTRANS))THEN
        VISS=VISWD
      ELSE    
        VISS=VISCOS
      ENDIF   
C
      IF(IXYZ.EQ.1)THEN
        DN  =ABS(XC(IWD)-XC(IWD+IJDN))
        TMPV=ABS(V(IJB)-V(IJP))
        TMPW=ABS(W(IJB)-W(IJP))
        TMP =SQRT(TMPV*TMPV+TMPW*TMPW)
      ELSE IF(IXYZ.EQ.2)THEN
        DN  =ABS(YC(JWD)-YC(JWD+IJDN))
        TMPU=ABS(U(IJB)-U(IJP))
        TMPW=ABS(W(IJB)-W(IJP))
        TMP =SQRT(TMPU*TMPU+TMPW*TMPW)
      ELSE
        DN  =ABS(ZC(KWD)-ZC(KWD+IJDN))
        TMPU=ABS(U(IJB)-U(IJP))
        TMPV=ABS(V(IJB)-V(IJP))
        TMP =SQRT(TMPU*TMPU+TMPV*TMPV)
      ENDIF   
      TAU=VISS*TMP/DN
      GEN(IJP)=TAU*CMU25*SQRT(TE(IJP))/(DN*CAPPA)
      SU(IJP)=SU(IJP)+GEN(IJP)*VOL(IJP)*DENSIT
C
      RETURN  
      END     
C
C#############################################################
      SUBROUTINE BCED
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'logic.inc'
      DOUBLE PRECISION VISC,D,DN,CP,CE,CT
C
C.....SOUTH BOUNDARY
C
      J=2
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCS(IK),VISC,YPLS(IK),VSWS(IK))
        D=VISC*(X(I)-X(I-1))*(Z(K)-Z(K-1))/(YC(J)-YC(J-1))
        AS(IJK)=-D
        ED(IJK-NZ)=ED(IJK)
      END DO
      END DO
C
C.....NORTH BOUNDARY
C
      J=NYM 
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCN(IK),VISC,YPLN(IK),VSWN(IK))
        D=VISC*(X(I)-X(I-1))*(Z(K)-Z(K-1))/(YC(J+1)-YC(J))
        AN(IJK)=-D
        ED(IJK+NZ)=ED(IJK)
      END DO
      END DO
C
C.....WEST BOUNDARY
C
      IF(MYID.EQ.IROOT)THEN
      I=2
      II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        JK=(J-1)*NZ+K
        CALL VISWALL(IBCW(JK),VISC,YPLW(JK),VSWW(JK))
        D=VISC*(Y(J)-Y(J-1))*(Z(K)-Z(K-1))/(XC(I)-XC(I-1))
        CP=MAX(F1(IJK-NYZ),ZERO)
        AW(IJK)=-CP-D
      END DO
      END DO
      ENDIF
C
C.....EAST BOUNDARY
C
      IF(MYID.EQ.IPRN)THEN
      I=NXTM
      II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        JK=(J-1)*NZ+K
        CALL VISWALL(IBCE(JK),VISC,YPLE(JK),VSWE(JK))
        D=VISC*(Y(J)-Y(J-1))*(Z(K)+Z(K-1))/(XC(I+1)-XC(I))
        CE=MIN(F1(IJK),ZERO)
        AE(IJK)= CE-D
      END DO
      END DO
      ENDIF
C
C.....BOTTOM BOUNDARY
C
      K=2
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCB(IJ),VISC,YPLB(IJ),VSWB(IJ))
        D=VISC*(X(I)-X(I-1))*(Y(J)-Y(J-1))/(ZC(K)-ZC(K-1))
        DN=ABS(ZC(K)-Z(K-1))
        CALL EDWALL(IJK,DN)
      END DO
      END DO
C
C.....TOP BOUNDARY
C
      K=NZM 
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCT(IJ),VISC,YPLT(IJ),VSWT(IJ))
        D=VISC*(X(I)-X(I-1))*(Y(J)-Y(J-1))/(Z(K)-ZC(K))
        CT=MIN(F3(IJK),ZERO)
        AT(IJK)= CT-D
      END DO
      END DO
C
      RETURN
      END
C
C#############################################################
      SUBROUTINE BCTE
C#############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'logic.inc'
      DOUBLE PRECISION VISC,D,CP,CE,CT
C
C.....SOUTH BOUNDARY
C
      J=2
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCS(IK),VISC,YPLS(IK),VSWS(IK))
        D=VISC*(X(I)-X(I-1))*(Z(K)-Z(K-1))/(YC(J)-YC(J-1))
        AS(IJK)=-D
        TE(IJK-NZ)=TE(IJK)
      END DO
      END DO
C
C.....NORTH BOUNDARY
C
      J=NYM 
      DO I=IXBN,IXED
        II=I-IREF
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        IK=(II-1)*NZ+K
        CALL VISWALL(IBCN(IK),VISC,YPLN(IK),VSWN(IK))
        D=VISC*(X(I)-X(I-1))*(Z(K)-Z(K-1))/(YC(J+1)-YC(J))
        AN(IJK)=-D
        TE(IJK+NZ)=TE(IJK)
      END DO
      END DO
C       
C.....WEST BOUNDARY
C       
      IF(MYID.EQ.ROOT)THEN
      I=2
      II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        JK=(J-1)*NZ+K
        CALL VISWALL(IBCW(JK),VISC,YPLW(JK),VSWW(JK))
        D=VISC*(Y(J)-Y(J-1))*(Z(K)-Z(K-1))/(XC(I)-XC(I-1))
        CP=MAX(F1(IJK-NYZ),ZERO)
        AW(IJK)=-CP-D
      END DO
      END DO
      ENDIF
C 
C.....EAST BOUNDARY
C
      IF(MYID.EQ.IPRN)THEN
      I=NXTM 
      II=I-IREF
      DO J=2,NYM 
      DO K=2,NZM 
        IJK=LI(II)+LJ(J)+K
        JK=(J-1)*NZ+K
        CALL VISWALL(IBCE(JK),VISC,YPLE(JK),VSWE(JK))
        D=VISC*(Y(J)-Y(J-1))*(Z(K)-Z(K-1))/(XC(I+1)-XC(I))
        CE=MIN(F1(IJK),ZERO)
        AE(IJK)= CE-D
      END DO
      END DO
      ENDIF
C
C.....BOTTOM BOUNDARY
C
      K=2
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCB(IJ),VISC,YPLB(IJ),VSWB(IJ))
        D=VISC*(X(I)-X(I-1))*(Y(J)-Y(J-1))/(ZC(K)-ZC(K-1))
        CALL TEWALL(IJK,I,J,K,-1,3,YPLB(IJ),VSWB(IJ))
      END DO
      END DO
C
C.....TOP BOUNDARY
C
      K=NZM 
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        CALL VISWALL(IBCT(IJ),VISC,YPLT(IJ),VSWT(IJ))
        D=VISC*(X(I)-X(I-1))*(Y(J)-Y(J-1))/(ZC(K+1)-ZC(K))
        CT=MIN(F3(IJK),ZERO)
        AT(IJK)= CT-D
      END DO
      END DO
      RETURN
      END
C
C############################################################
      SUBROUTINE MODVIS
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'coef.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION TMP,TVISMAX,DN,YMIN(6),YMAX(6)
      DOUBLE PRECISION TA(1),TI(1),TAA(1),TII(1) 
C
      DO I=IXBN,IXED
        II=I-IREF
        DO J=2,NYM 
        DO K=2,NZM 
            IJK=LI(II)+LJ(J)+K
            VIS(IJK)=VISCOS
            IF(LVIST)THEN
              TMP=CMU*DENSIT*TE(IJK)*TE(IJK)/(ED(IJK)+SMALL)
              TVISMAX=VISCOS*SFVIST
              TMP=MIN(TMP,TVISMAX)
              VIST(IJK)=TMP*URFM(IVIST)+VIST(IJK)*(UNITY-URFM(IVIST))
              VIS(IJK)=VIS(IJK)+VIST(IJK) 
            ENDIF
        END DO
        END DO
      END DO
C
      IF(NPRO.GT.1)THEN
        IF(LVIST) CALL EXCHFI(ISIZE,VIST,IVIST*10)
        CALL EXCHFI(ISIZE,VIS,IVIS*10)
      ENDIF
C
c     I =1
c     II=NX
c     J=1
c     JJ=NY
c     K=1
c     KK=NZ
c     CALL PRINTOUT(I,II,J,JJ,K,KK,NXYZ,VIS,"VIS")
c
      IF(LVIST)THEN
C
C.....BOTTOM BOUNDARY
C
      K=2
      DO I=IXBN,IXED
        II=I-IREF
      DO J=2,NYM 
        IJK=LI(II)+LJ(J)+K
        IJ=(II-1)*NY+J
        DN=ABS(ZC(K)-Z(K-1))
        CALL WALLVIS(IJK,-1,DN,YPLB(IJ),VSWB(IJ))
      END DO
      END DO
C       
C.....FIND YPL MAX AND MIN
C
 111  CALL FINDYPL(NXY,YPLB,YMIN(1),YMAX(1))
      TA(1)=YMAX(1)
      TI(1)=YMIN(1)
C
      IF(NPRO.GT.1)THEN
        CALL MPI_ALLREDUCE(TA,TAA,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     *                     MPI_COMM_WORLD,IERR)
        CALL MPI_ALLREDUCE(TI,TII,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     *                     MPI_COMM_WORLD,IERR)
      ELSE
        TAA(1)=TA(1)
        TII(1)=TI(1)
      ENDIF
C
      YPMAX=TAA(1)
      YPMIN=TII(1)
C
      ENDIF
C
      RETURN
      END
C
C############################################################
      SUBROUTINE FINDYPL(IDIM,YPL,YMIN,YMAX)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'rcont.inc'
      DOUBLE PRECISION YPL(IDIM),TI,TA,YMIN,YMAX
C
      TI=GREAT
      TA=ZERO
      DO I=1,IDIM
        IF(YPL(I).GT.TA) TA=YPL(I) 
        IF(YPL(I).LT.TI.AND.YPL(I).GT.ZERO) TI=YPL(I) 
      END DO 
      YMAX=TA
      YMIN=TI
C
      RETURN
      END
C
C############################################################
      SUBROUTINE WALLVIS(IJK,IJD,DN,YPLWD,VISWD)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      DOUBLE PRECISION DN,VISC,CK,VISCW,YPLWD,VISWD
C
C.....UPDATE THE DIMENSIONLESS DISTANCE FROM THE WALL
C 
      IJP=IJK
      IJB=IJK+IJD
      VISC=VISCOS
      CK=CMU25*SQRT(MAX(ZERO,TE(IJP)))
      YPLWD=DENSIT*CK*DN/VISC
      IF(YPLWD.LT.CTRANS)THEN
        VISCW=VISC
      ELSE
        VISCW=YPLWD*VISC*CAPPA/LOG(ELOG*YPLWD)
      ENDIF
      VISWD=VISCW
      VIS(IJB)=VISCW
C     write(*,*) "WWW",IJB,VIS(IJB)
C
      RETURN
      END
C
C############################################################
      SUBROUTINE GRADXYZ(FT,DFX,DFY,DFZ)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'var.inc'
      DOUBLE PRECISION DFX(NXYZ),DFY(NXYZ),DFZ(NXYZ),FT(NXYZ)
      DOUBLE PRECISION FXI,FXIM,PE,PW
      DOUBLE PRECISION FYJ,FYJM,PN,PS
      DOUBLE PRECISION FZK,FZKM,PT,PB
      DOUBLE PRECISION DX,DY,DZ
C
      DO I=IXBN,IXED
        II=I-IREF
        DX=X(I)-X(I-1)
C
        DO J=2,NYM 
          DY=Y(J)-Y(J-1)
C
        DO K=2,NZM 
          DZ=Z(K)-Z(K-1)
C
          IJK=LI(II)+LJ(J)+K
C
          FXI =FX(I)      
          FXIM=FX(I-1) 
          FYJ =FY(J)     
          FYJM=FY(J-1)   
          FZK =FZ(K)     
          FZKM=FZ(K-1)   
C
          PE=FT(IJK+NYZ)*FXI +FT(IJK)    *(UNITY-FXI)
          PW=FT(IJK)    *FXIM+FT(IJK-NYZ)*(UNITY-FXIM)
          PN=FT(IJK+NZ) *FYJ +FT(IJK)    *(UNITY-FYJ)
          PS=FT(IJK)    *FYJM+FT(IJK-NZ) *(UNITY-FYJM)
          PT=FT(IJK+1)  *FZK +FT(IJK)    *(UNITY-FZK)
          PB=FT(IJK)    *FZKM+FT(IJK-1)  *(UNITY-FZKM)
C
          DFX(IJK)=(PE-PW)/DX
          DFY(IJK)=(PN-PS)/DY
          DFZ(IJK)=(PT-PB)/DZ
C
        END DO
        END DO
      END DO
C
      RETURN
      END
C
C############################################################
      SUBROUTINE GETVEL(IDX,SP,VP,ISZ,IDM)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'var.inc'
      DIMENSION IDX(4),IC(IDM)
      DOUBLE PRECISION SP(4),VP(4),TP(4),DT(IDM),TMP
C
      IX=IDX(1)
      JX=IDX(2)
      KX=IDX(3)
C
      IT=0
      DO I=-ISZ,ISZ
      DO J=-ISZ,ISZ
      DO K=-ISZ,ISZ
         IT=IT+1
         IC(IT)=LI(IX+I)+LJ(JX+J)+(KX+K)
         TP(1)=XC(IX+I)
         TP(2)=YC(JX+J)
         TP(3)=ZC(KX+K)
         CALL DISTFIND(SP,TP,DT(IT))
      END DO
      END DO
      END DO
C
      DO I=1,4
        TP(I)=ZERO
      END DO
C
      DO I=1,IDM
        IT=IC(I)
        TMP=UNITY/DT(I)
        TP(1)=TP(1)+U(IT)*TMP
        TP(2)=TP(2)+V(IT)*TMP
        TP(3)=TP(3)+W(IT)*TMP
        TP(4)=TP(4)+TMP
      END DO
C
      DO I=1,3
        VP(I)=TP(I)/TP(4)
      END DO
C
      RETURN
      END
C 
C#######################################################################
      SUBROUTINE FINDWINDSPEED
C#######################################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
C
      DOUBLE PRECISION WSD(2),TWS(MWT),TMP
C
      DO IWT=1,NWT
C
       IF(LWS(IWT))THEN
         DO I=1,2
           II=IIX(IWT)-IREF+I-1
C
           ID=0
           WSD(I)=ZERO
           DO J=IJL(IWT),IJU(IWT)+1
           DO K=IKL(IWT),IKU(IWT)+1
             ID=ID+1
             IJK=LI(II)+LJ(J)+K
             WSD(I)=WSD(I)+U(IJK)*STWS(IWT,ID)
           END DO
           END DO
C
         END DO
         IX=IIX(IWT)
         TMP=(WBX(IWT)-XC(IX))/(XC(IX+1)-XC(IX))
         VSCA(IWT)=WSD(2)*TMP+WSD(1)*(UNITY-TMP)
       ELSE
         VSCA(IWT)=ZERO
       ENDIF
C 
      END DO
C
      IF(NPRO.GT.1)THEN
        CALL MPI_ALLREDUCE(VSCA,TWS,NWT,MPI_DOUBLE_PRECISION,MPI_SUM,
     *                     MPI_COMM_WORLD,IERR)
        DO IWT=1,NWT
          VSCA(IWT)=TWS(IWT)/XISUM(IWT)
        END DO
      ENDIF
C
      RETURN
      END
C 
C#######################################################################
      SUBROUTINE ADM
C#######################################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
C
      DOUBLE PRECISION TMPAR,ROTA,R3,THRUST,ROTR,PI2,START,FINISH
      DOUBLE PRECISION VP(4),TMPAF,TMPTF,TMPRR,DR
      DOUBLE PRECISION THETA,COSTH,SINTH,DTH,TMPWTE,TMPTE
C     DOUBLE PRECISION WTCP,WTCT,WTLA,WTCG,U3,U2
C     DOUBLE PRECISION RPMSTAR,VSTAR 
      DIMENSION IJK(4)
C
      PI2=PI*PI
      DO I=1,NXYZ
        BDY1(I)=ZERO
        BDY2(I)=ZERO
        BDY3(I)=ZERO
        WTE(I)=ZERO
      END DO
C
      CALL FINDWINDSPEED
C
      DO IWT=1,NWT
C      IF(.NOT.LWT(IWT)) GOTO 123
C
       ITY=IWTY(IWT)
       ROTR=ROTD(ITY)*HALF
       ROTA=PI*ROTR*ROTR
       R3=ROTR*ROTR*ROTR
C
       IF(LRATED)THEN
         POWER(IWT)=PRATED(ITY)
         RPM(IWT)=WRATED(ITY)
         VEL(IWT)=UIN                                                   
         TMPWTE=TERATED(ITY)
       ELSE
         IF((VSCA(IWT).LT.PCS(ITY,1)).OR.
     *      ((VSCA(IWT).GT.PCS(ITY,NPC(ITY)))))THEN
           POWER(IWT)=ZERO
           RPM(IWT)=ZERO
           TMPWTE=ZERO
           GOTO 123
         ENDIF
         CALL GETPOINT(IWT,NPC,PCS,PCP,VSCA(IWT),POWER(IWT))
         CALL GETPOINT(IWT,NPC,PCS,PCR,VSCA(IWT),RPM(IWT))
         CALL GETPOINT(IWT,NPC,PCS,WDS,VSCA(IWT),VEL(IWT))
         CALL GETPOINT(IWT,NPC,PCS,TTE,VSCA(IWT),TMPWTE)
       ENDIF
       THRUST=POWER(IWT)/VSCA(IWT)
C      WRITE(*,*) (TTE(IWT,I),I=1,NPC(IWT))
C      WRITE(*,*) IWT,TMPWTE,VSCA(IWT)
C
C      U2=VSCA(IWT)*VSCA(IWT)
C      U3=U2*VSCA(IWT)
C      RPMSTAR=ABS(RPM(IWT)/WRATED(ITY))
C      VSTAR=ABS(VSCA(IWT)/VRATED(ITY))
C
C      TMP=HALF*DENSIT*ROTA
C      WTCP=POWER(IWT)/(TMP*U3)
C      WTCT=THRUST/(TMP*U2)
C      WTLA=RPM(IWT)*ROTR/VSCA(IWT)*PI*THIRD*TENTH
C      WTCG=(WTCP**EXPCP(ITY))*(WTCT**EXPCT(ITY))*(WTLA**EXPLA(ITY))
C      TMPWTE=WTCG*CGSCA(ITY)*TMP*U3
C      TMPWTE=EXPCP(ITY)*(VSTAR**EXPCT(ITY))*(RPMSTAR**EXPLA(ITY))
C
CCC    RADIAL DISTRIBUTION OF BLADE LOADING
C
       IF(LBLADE)THEN
C
         TMPAF=THRUST/ROTA
         TMPTF=45D0*POWER(IWT)/(PI2*RPM(IWT)*R3)
         TMPTE=TMPWTE/ROTA
         NRA=NWTRA(IWT)
         NAZ=NWTAZ(IWT)
         DTH=TWO*PI/NAZ
         DR=ROTR/NRA
c      write(*,*) IWT,TMPTE,TMPAF,TMPTF,VSCA(IWT)
C
         IAD=0                                                          
         DO IAZ=1,NAZ
           THETA=(2*IAZ-1)*DTH*HALF
           SINTH=SIN(THETA)
           COSTH=COS(THETA)
           DO IRA=1,NRA
             IAD=IAD+1                                                  
             TMPAR= HALF*DTH*DR*DR*(2*IRA-1)
             VP(1)=-TMPAR*TMPAF
             TMPRR= TMPAR*TMPTF
             VP(2)= TMPRR*SINTH
             VP(3)=-TMPRR*COSTH
             VP(4)= TMPAR*TMPTE
             IJK(1)=ICADP(IWT,IAD,1)                                    
             IJK(2)=ICADP(IWT,IAD,2)                                    
             IJK(3)=ICADP(IWT,IAD,3)                                    
             CALL DISTBDF(IWT,IAD,IJK,VP,ICZ)
           END DO
         END DO
C
       ENDIF
C
 123  END DO  
C
      END
C
C############################################################
      SUBROUTINE DISTFIND(SP,TP,DIST) 
C############################################################
      INCLUDE 'param.inc'
      DOUBLE PRECISION SP(4),TP(4),DIST
      DIST=ZERO
      DO I=1,3
        TP(4)=SP(I)-TP(I)
        DIST=DIST+TP(4)*TP(4)
      END DO
      DIST=SQRT(DIST)
      RETURN
      END
C
C############################################################
      SUBROUTINE DISTBDF(IWT,IAD,IDX,VP,ISZ)
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'var.inc'
      INCLUDE 'wind.inc'
      DIMENSION IDX(4)
      DOUBLE PRECISION VP(4)
      LOGICAL LFLAG
C
      IX=IDX(1)
      JX=IDX(2)
      KX=IDX(3)
C
      IST=0
      DO I=-ISZ,ISZ+1
       II=IX+I
       LFLAG=(II.GE.IXBN).AND.(II.LE.IXED)
       IF(LFLAG)THEN
         II=IX-IREF
         DO J=-ISZ,ISZ+1
         DO K=-ISZ,ISZ+1
            IST=IST+1
            IT=LI(II+I)+LJ(JX+J)+(KX+K)
            TMP=STADP(IWT,IAD,IST)
c      write(*,*) IST,II+I,JX+J,KX+K,IT,TMP,(VP(III),III=1,4)
            BDY1(IT)=BDY1(IT)+TMP*VP(1)
            BDY2(IT)=BDY2(IT)+TMP*VP(2)
            BDY3(IT)=BDY3(IT)+TMP*VP(3)
            WTE(IT) =WTE(IT) +TMP*VP(4)
         END DO
         END DO
       ELSE
         IST=IST+(ISZ+2)*(ISZ+2)
       ENDIF
      END DO
C
      RETURN
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
C############################################################
      SUBROUTINE GETPOINT(ID,NT,XR,XF,XMYR,XMYF) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'wind.inc'
      INTEGER NT(MTY)
      DOUBLE PRECISION XR(MTY,MNP),XF(MTY,MNP),XMYR,XMYF,TMP
C
      ITY=IWTY(ID)
      IF(XMYR.LE.XR(ITY,1))THEN
       XMYF=XF(ITY,1)
      ELSE
       DO I=2,NT(ITY)-1
        IF(XMYR.LE.XR(ITY,I)) GOTO 678 
       END DO
       I=NT(ITY)
  678  TMP=(XMYR-XR(ITY,I-1))/(XR(ITY,I)-XR(ITY,I-1)) 
       XMYF=XF(ITY,I-1)+TMP*(XF(ITY,I)-XF(ITY,I-1))
      ENDIF
C
      RETURN
      END
C############################################################
      SUBROUTINE UNITMX(A) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION A(4,4) 
C
      DO I=1,4
      DO J=1,4
         IF(I.EQ.J) THEN
           A(I,J)=UNITY
         ELSE
           A(I,J)=ZERO
         ENDIF
      END DO
      END DO
C
      RETURN
      END
C############################################################
      SUBROUTINE MULTPM(A,B,C) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION A(4,4),B(4,4),C(4,4) 
      DO I=1,4
      DO J=1,4
        C(I,J)=ZERO
        DO K=1,4
          C(I,J)=C(I,J)+A(I,K)*B(K,J)
        END DO
      END DO
      END DO
C
      RETURN
      END
C############################################################
      SUBROUTINE MULTPV(A,B,C) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION A(4,4),B(4),C(4) 
C
      DO I=1,4
        C(I)=ZERO
        DO K=1,4
          C(I)=C(I)+A(I,K)*B(K)
        END DO
      END DO
C
      RETURN
      END
C############################################################
      SUBROUTINE TRANSL(T,XMX) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION T(3),TMX(4,4),XMX(4,4),YMX(4,4)
C
      CALL UNITMX(TMX)
      TMX(1,4)=T(1)
      TMX(2,4)=T(2)
      TMX(3,4)=T(3)
      CALL MULTPM(TMX,XMX,YMX)
      CALL COPYMX(YMX,XMX)
C
      RETURN
      END
C############################################################
      SUBROUTINE COPYMX(A,B) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION A(4,4),B(4,4)
      DO I=1,4
      DO J=1,4
        B(I,J)=A(I,J)
      END DO
      END DO
      RETURN
      END
C############################################################
      SUBROUTINE PRINTM(A) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION A(4,4)
      IF(MYID.EQ.IPRN)THEN
        WRITE(*,*) ('*',J=1,40)
        DO I=1,4
          WRITE(*,*) (A(I,J),J=1,4)
        END DO
        WRITE(*,*) ('*',J=1,40)
      ENDIF
      RETURN
      END
C############################################################
      SUBROUTINE ROTATE(IAXIS,RAD,XMX) 
C############################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      DOUBLE PRECISION XMX(4,4),RAD
      DOUBLE PRECISION TMX(4,4),YMX(4,4),COSA,SINA
C
      CALL UNITMX(TMX)
      IF(IAXIS.LT.1.OR.IAXIS.GT.3)THEN
        IF(MYID.EQ.IPRN)THEN
          WRITE(*,*) "(E) IAXIS<0 or IAXIS>3."
        ENDIF
        STOP
      ENDIF
      COSA=COS(RAD)
      SINA=SIN(RAD)
      IF(IAXIS.EQ.1)THEN
        TMX(2,2)= COSA 
        TMX(2,3)=-SINA 
        TMX(3,2)= SINA 
        TMX(3,3)= COSA 
      ELSE IF(IAXIS.EQ.2)THEN
        TMX(1,1)= COSA 
        TMX(1,3)= SINA 
        TMX(3,1)=-SINA 
        TMX(3,3)= COSA 
      ELSE IF(IAXIS.EQ.3)THEN
        TMX(1,1)= COSA 
        TMX(1,2)=-SINA 
        TMX(2,1)= SINA 
        TMX(2,2)= COSA 
      ELSE
        WRITE(*,*) '(E) Unknown Rotation Axis.'
        STOP
      ENDIF
      CALL MULTPM(TMX,XMX,YMX)
      CALL COPYMX(YMX,XMX)
C
      RETURN
      END
C############################################################
      SUBROUTINE FINDFILENAME
C############################################################
      INCLUDE 'param.inc'
C
      FILCTL=TRIM(FCASE)//'.ctl'
      FILCHK=TRIM(FCASE)//'.chk'
      FILOUT=TRIM(FCASE)//'.out'
      FILGRD=TRIM(FCASE)//'.grd'
      FILPWR=TRIM(FCASE)//'.pwr'
      FILWTD=TRIM(FCASE)//'.wtd'
C
      RETURN
      END
C############################################################
      SUBROUTINE FINDMINMAX(NARRAY,XARRAY,XMIN,XMAX)
C############################################################
      INCLUDE 'param.inc'
      DOUBLE PRECISION XARRAY(NARRAY),XMIN,XMAX
C
      XMIN=XARRAY(1)
      XMAX=XARRAY(1)
      DO I=2,NARRAY
         IF(XARRAY(I).GT.XMAX)THEN
          XMAX=XARRAY(I) 
         ELSE 
          IF(XARRAY(I).LT.XMAX) XMIN=XARRAY(I) 
         ENDIF
      END DO
C
      RETURN
      END
C##########################################################
      SUBROUTINE INITADM                                   
C##########################################################
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
      INCLUDE 'wind.inc'
      INCLUDE 'mpif.h'
C  
      DOUBLE PRECISION ROTR,RC(MAT),YY(5),ZZ(5),Y1,Z1,Y2,Z2
      DOUBLE PRECISION SP(4),TMPRR,RO(4),TY,TZ,AREA,XYZ(4)
      DOUBLE PRECISION YE,YW,ZT,ZB,DR,TP(4),DT(ICE),STA(4)
      DOUBLE PRECISION THETA,COSTH,SINTH,DTH,TMP,TMPA,ROTA,TAA
      DOUBLE PRECISION XWS(MWT),XSUM(MWT)
      DIMENSION IJK(4),JAP(5),KAP(5),IXAP(5)
      LOGICAL LINSIDE(MAT),LAND,LOR
C
      DO IWT=1,MWT
        XWS(IWT)=ZERO
        XSUM(IWT)=ZERO
        XISUM(IWT)=ZERO
      END DO
C
      IF(LRATED)THEN
        IF(NPRO.NE.1)THEN
          IF(MYID.EQ.IPRN) 
     *      WRITE(*,*) "(E) RATED Condition is allowed for NPRO=1."
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
      ENDIF
C
      IF(MYID.EQ.IPRN) WRITE(*,*) '(*) Checking WT Positions ...'
C
      DO IWT=1,NWT
        IXL=IXBN-1-2*ICZ
        IXU=IXED+1+2*ICZ
        IF(IXL.LT.1) IXL=1
        IF(IXU.GT.NXT) IXU=NXT
        IF(WBX(IWT).GE.XC(IXL).AND.WBX(IWT).LE.XC(IXU))THEN
          LWT(IWT)=.TRUE.
        ELSE
          LWT(IWT)=.FALSE.
        ENDIF
        IXL=IXBN-1
        IXU=IXED+1
        IF(IXL.LT.1) IXL=1
        IF(IXU.GT.NXT) IXU=NXT
        IF(WBX(IWT).GE.XC(IXL).AND.WBX(IWT).LE.XC(IXU))THEN
          LWS(IWT)=.TRUE.
          XWS(IWT)=UNITY
        ELSE
          LWS(IWT)=.FALSE.
          XWS(IWT)=ZERO
        ENDIF
        IF(IWTY(IWT).LT.1.OR.IWTY(IWT).GT.NWTD)THEN
          IF(MYID.EQ.IPRN) WRITE(*,*) "(*) NWTD Error!"
          STOP
        ENDIF
      END DO
C
      DO IWT=1,NWT
C
        ITY=IWTY(IWT)
        ROTR=ROTD(ITY)*HALF
        ROTA=PI*ROTR*ROTR
        RO(1)=WBX(IWT)
        RO(2)=WBY(IWT)
        RO(3)=WBZ(IWT)+HUBH(ITY)
        RO(4)=ZERO
        YW=RO(2)-ROTR
        YE=RO(2)+ROTR
        ZB=RO(3)-ROTR
        ZT=RO(3)+ROTR
        IXL=2+ICZ*2
        IXU=NXTM-ICZ*2
C
        IF(RO(1).LT.X(IXL).OR.RO(1).GT.X(IXU))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',IWT,' is outside X range! '
             CALL MPI_FINALIZE(IERR)
             STOP
          ENDIF
        ENDIF
        IF(YW.LT.Y(1).OR.YE.GT.Y(NY))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',IWT,' is outside Y range!' 
             CALL MPI_FINALIZE(IERR)
             STOP
          ENDIF
        ENDIF
        IF(ZB.LT.Z(1).OR.ZT.GT.Z(NZ))THEN
          IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',IWT,' is outside Z range!' 
             CALL MPI_FINALIZE(IERR)
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
        NAT=NAJ*NAK
C
        IF(NAT.GT.MAT)THEN
          IF(MYID.EQ.IPRN) WRITE(*,*) "(*) MAT Error!"
          STOP
        ENDIF
C
        DO IAP=1,NAT
          STWS(IWT,IAP)=ZERO
        END DO 
        IAP=0
        SP(1)=RO(1)
        DO J=IJL(IWT),IJU(IWT)+1
        DO K=IKL(IWT),IKU(IWT)+1
          IAP=IAP+1
          SP(2)=YC(J)
          SP(3)=ZC(K)
          CALL DISTFIND(RO,SP,RC(IAP))
          IF(RC(IAP).GT.ROTR)THEN
            LINSIDE(IAP)=.FALSE.
          ELSE
            LINSIDE(IAP)=.TRUE.
          ENDIF
        END DO
        END DO
c
        TMPA=ZERO
        TAA=ZERO
        DO J=IJL(IWT),IJU(IWT)
        DO K=IKL(IWT),IKU(IWT)
          TAA=TAA+(YC(J+1)-YC(J))*(ZC(K+1)-ZC(K)) 
          IAP=(J-IJL(IWT))*NAK+K-IKL(IWT)+1
          IXAP(1)=IAP
          IXAP(2)=IAP+NAK
          IXAP(3)=IAP+1+NAK
          IXAP(4)=IAP+1
          IXAP(5)=IAP
          LAND=LINSIDE(IXAP(1)).AND.LINSIDE(IXAP(2)).AND.
     *         LINSIDE(IXAP(3)).AND.LINSIDE(IXAP(4))
          LOR =LINSIDE(IXAP(1)).OR.LINSIDE(IXAP(2)).OR.
     *         LINSIDE(IXAP(3)).OR.LINSIDE(IXAP(4))
          IF(LAND)THEN
            AREA=(YC(J+1)-YC(J))*(ZC(K+1)-ZC(K))
            TMP=AREA*QUARTER
            DO ID=1,4
              STWS(IWT,IXAP(ID))=STWS(IWT,IXAP(ID))+TMP
            END DO
          ELSE IF(LOR)THEN
            JAP(1)=J
            KAP(1)=K
            JAP(2)=J+1
            KAP(2)=K
            JAP(3)=J+1
            KAP(3)=K+1
            JAP(4)=J
            KAP(4)=K+1
            JAP(5)=J
            KAP(5)=K
            II=0
            DO ID=1,4 
              IF(LINSIDE(IXAP(ID)))THEN
                II=II+1
                YY(II)=YC(JAP(ID))  
                ZZ(II)=ZC(KAP(ID))  
              ENDIF
              IF(LINSIDE(IXAP(ID)).NEQV.LINSIDE(IXAP(ID+1)))THEN
                II=II+1
                IF(ID.EQ.1.OR.ID.EQ.3)THEN
                  TZ=ZC(KAP(ID))
                  TMP=TZ-RO(3)
                  TMP=SQRT(ROTR*ROTR-TMP*TMP)
                  TY=RO(2)-TMP
                  IF(((TY-YC(JAP(ID)))*(TY-YC(JAP(ID+1)))).GT.ZERO)THEN
                    TY=RO(2)+TMP
                  ENDIF 
                ELSE
                  TY=YC(JAP(ID))
                  TMP=TY-RO(2)
                  TMP=SQRT(ROTR*ROTR-TMP*TMP)
                  TZ=RO(3)-TMP
                  IF(((TZ-ZC(KAP(ID)))*(TZ-ZC(KAP(ID+1)))).GT.ZERO)THEN
                    TZ=RO(3)+TMP
                  ENDIF 
                ENDIF
                YY(II)=TY
                ZZ(II)=TZ
              ENDIF
            END DO 
            AREA=ZERO         
            DO ID=2,II-1 
              Y1=YY(ID)-YY(1)
              Z1=ZZ(ID)-ZZ(1)
              Y2=YY(ID+1)-YY(1)
              Z2=ZZ(ID+1)-ZZ(1)
              AREA=AREA+HALF*ABS(Y1*Z2-Y2*Z1) 
            END DO
            Y1=ZERO
            Z1=ZERO
            DO ID=1,II 
              Y1=Y1+YY(ID)
              Z1=Z1+ZZ(ID)
            END DO
            Y1=Y1/II
            Z1=Z1/II
            XYZ(1)=ZERO
            SP(1)=ZERO
            SP(2)=Y1
            SP(3)=Z1
            TMP=ZERO
            DO ID=1,4
              XYZ(2)=YC(JAP(ID))
              XYZ(3)=ZC(KAP(ID))
              CALL DISTFIND(XYZ,SP,STA(ID))
              STA(ID)=UNITY/(STA(ID)+SMALL)
              TMP=TMP+STA(ID)
            END DO
            DO ID=1,4
              STA(ID)=STA(ID)/TMP
            END DO
            TMP=AREA
            DO ID=1,4
              STWS(IWT,IXAP(ID))=STWS(IWT,IXAP(ID))+TMP*STA(ID)
            END DO
          ELSE
            AREA=ZERO
          ENDIF
          TMPA=TMPA+AREA
C
        END DO
        END DO
C
        DO IAP=1,NAT
          STWS(IWT,IAP)=STWS(IWT,IAP)/TMPA
        END DO
C
      END DO
C
      IF(NPRO.GT.1)THEN
        CALL MPI_ALLREDUCE(XWS,XSUM,MWT,MPI_DOUBLE_PRECISION,MPI_SUM,
     *                     MPI_COMM_WORLD,IERR)
        DO IWT=1,NWT
          IF(XSUM(IWT).LT.HALF)THEN
           IF(MYID.EQ.IPRN) THEN
             WRITE(*,*) '(E) WT ',IWT,': Wind Speed Location Error!' 
             CALL MPI_FINALIZE(IERR)
             STOP
           ENDIF
          ENDIF
          XISUM(IWT)=XSUM(IWT)
        END DO
      ELSE
        DO IWT=1,NWT
          XISUM(IWT)=UNITY
        END DO
      ENDIF
C
C     Compute the closet grid near the actuator point                       
C    
      SP(4)=ZERO
      DO IWT=1,NWT
C
        ITY=IWTY(IWT)
        ROTR=ROTD(ITY)*HALF
        RO(1)=WBX(IWT)
        RO(2)=WBY(IWT)
        RO(3)=WBZ(IWT)+HUBH(ITY)
C
        SP(1)=RO(1)
        NRA=(IJU(IWT)-IJL(IWT)+IKU(IWT)-IKL(IWT)+3)/2
        NAZ=(IJU(IWT)-IJL(IWT)+IKU(IWT)-IKL(IWT)+2)*2
        NWTRA(IWT)=NRA
        NWTAZ(IWT)=NAZ
        DTH=TWO*PI/NAZ                                                 
        DR=ROTR/NRA  
        IF(MYID.EQ.IPRN) WRITE(*,17) "(",IWT,") NRA =",NRA,", NAZ =",NAZ
  17    FORMAT(1X,A2,I4,A7,I3,A7,I3)
C
        IJK(1)=INDEXFIND(NXT,XC,SP(1))
        IAD=0
        DO IAZ=1,NAZ
          THETA=(2*IAZ-1)*DTH*HALF   
          SINTH=SIN(THETA)
          COSTH=COS(THETA)
          DO IRA=1,NRA
            IAD=IAD+1
            TMPRR=HALF*DR*(2*IRA-1)
            SP(2)= RO(2)+TMPRR*COSTH
            SP(3)= RO(3)+TMPRR*SINTH
            IJK(2)=INDEXFIND(NY ,YC,SP(2))
            IJK(3)=INDEXFIND(NZ ,ZC,SP(3))
            DO ID=1,3
              ICADP(IWT,IAD,ID)=IJK(ID)
            END DO
C
            IX=IJK(1)
            JX=IJK(2)
            KX=IJK(3)
            ISP=ICZ+1
            EPS=(XC(IX+ISP)-XC(IX-ICZ)
     *          +YC(JX+ISP)-YC(JX-ICZ)
     *          +ZC(KX+ISP)-ZC(KX-ICZ))/THREE
C
            IST=0
            DO I=-ICZ,ICZ+1
            DO J=-ICZ,ICZ+1
            DO K=-ICZ,ICZ+1
              IST=IST+1
              TP(1)=XC(IX+I)
              TP(2)=YC(JX+J)
              TP(3)=ZC(KX+K)
              CALL DISTFIND(SP,TP,DT(IST))
              DT(IST)=DT(IST)/EPS
              DT(IST)=EXP(-DT(IST)*DT(IST))
C             write(*,*) IWT,I,J,K,IST,DT(IST),EPS,ICZ
            END DO
            END DO
            END DO
C
            TP(4)=ZERO
            DO IST=1,ICE
C             write(*,*) IWT,IAD,IST,DT(IST)
              TP(4)=TP(4)+DT(IST)
            END DO
C           write(*,*) IWT,IAD,TP(4),ICE
            DO IST=1,ICE
              STADP(IWT,IAD,IST)=DT(IST)/TP(4)
c             IF(MYID.EQ.IPRN) write(*,*) IWT,IAD,IST,STADP(IWT,IAD,IST) 
            END DO
            TMP=ZERO
            DO IST=1,ICE
              TMP=TMP+STADP(IWT,IAD,IST)
            END DO
c           IF(MYID.EQ.IPRN) write(*,*) IWT,IAZ,IRA,IAD,TMP,UNITY 
C
          END DO
        END DO
C
  123 END DO
      RETURN
C 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C234567890
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRINTOUT(ID,I1,I2,J1,J2,K1,K2,NN,XX,WORD)
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
C
      INTEGER I1,I2,J1,J2,K1,K2,I,J,K,IC
      DOUBLE PRECISION XX(NN),YY
      CHARACTER*3 WORD
c     write(*,*) I1,I2,J1,J2,K1,K2,NN
      YY=ZERO
      IC=0
      DO I=I1,I2
      DO J=J1,J2
      DO K=K1,K2
        IC=IC+1
        IJK=LI(I)+LJ(J)+K
        YY=YY+XX(IJK)
c       write(ID,*) I,J,K,IJK,XX(IJK),YY
      END DO
      END DO
      END DO
      write(ID,*) WORD,I2,J2,K2,IC,YY
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C234567890
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRINTALL(ID)
      INCLUDE 'param.inc'
      INCLUDE 'index.inc'
      INCLUDE 'logic.inc'
      INCLUDE 'rcont.inc'
      INCLUDE 'var.inc'
      INCLUDE 'geo.inc'
      INCLUDE 'iter.inc'
      INCLUDE 'init.inc'
      INCLUDE 'turb.inc'
C
      INTEGER ID
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,U,"UUU")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,V,"VVV")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,W,"WWW")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,P,"PPP")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,F1,"FF1")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,F2,"FF2")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,F3,"FF3")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,TE,"TEE")
      CALL PRINTOUT(ID,1,NX,1,NY,1,NZ,NXYZ,ED,"EDD")
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C234567890
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FILEERROR(FNAME)
      CHARACTER*128 FNAME
      WRITE(*,*) "(E) ",TRIM(FNAME),": file open error!"
      STOP
      RETURN
      END
