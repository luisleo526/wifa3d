      COMMON /TUR1/ TE(NXYZ),ED(NXYZ),GEN(NXYZ),VIST(NXYZ),VIS(NXYZ),
     *   YPLW(NYZ),YPLE(NYZ),YPLS(NXZ),YPLN(NXZ),YPLT(NXY),YPLB(NXY),
     *   VSWW(NYZ),VSWE(NYZ),VSWS(NXZ),VSWN(NXZ),VSWT(NXY),VSWB(NXY) 
      COMMON /TUR2/ CMU,CMU75,CAPPA,STEMP,STE,SED,CED1,CED2,ELOG,
     *              CTRANS,LVIST,IVIST,CMU25,YPMAX,YPMIN,SFVIST,
     *              EROT,ALAN
      LOGICAL LVIST
      DOUBLE PRECISION TE,ED,GEN,VIST,VIS,
     *   YPLW,YPLE,YPLS,YPLN,YPLT,YPLB,VSWW,VSWE,VSWS,VSWN,VSWT,VSWB
      DOUBLE PRECISION CMU,CMU75,CAPPA,STEMP,STE,SED,CED1,CED2,
     *                 ELOG,CTRANS,CMU25,YPMAX,YPMIN,SFVIST,
     *                 EROT,ALAN