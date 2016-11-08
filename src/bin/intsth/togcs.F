      SUBROUTINE TOGCS(IGC,LA,NIR,NIRL,ND,ITYP,LAMBDA,SLAB,ICT,
     1  KLASS,KATOM,KTYPE,KMIN,NPERSH,ISC,IPC,LOC,LOC2,LAB,CC)
C
C  THIS ROUTINE CONVERTS HONDO SYMMETRY ARRAYS TO PITZER'S AO-TO-SO
C  MATRICES (IGC) AND FORMS OTHER USEFUL ARRAYS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      integer mctu,mcsu,mgcsu,mxcf,mru,kaords,mstu,mxsym,mxat,mxtrn
      integer mxsh,mxpsym
      COMMON /OUT/    IOUT,ICHECK,IFLAG3,NSECT,NINTS
      COMMON /PSIZE/  MSU,KAORDS,MCONSU,MGCSU,MSTU,MRU,MCSU,MCTU,MCONU,
     1                MCU,MSFU,MGU,MSFRU,MNSFU,MPRU,MCXU,MCCU,MBLU,MRCRU
      COMMON /MAXDIM/ MXAT,MXTRN,MXPR,MXSH,MXSYM,MXCF,MXT,MXISO,MXPSYM
      COMMON /NSHEL/  NSHELL
      COMMON /SALCS/  NIRRED
      COMMON /SYMTRY/ NT
      COMMON /INFOA/  NAT,NUM,NQ,NX,NE,NA,NB
      COMMON /PACKQ/  I1,I2,I3,I4,I5
      COMMON /ISOPAC/ INDIN(120),INDOUT(30)
C
      DIMENSION II(30),CA(30,18),KCS(6)
C
      REAL*8 IGC(MCTU,MCSU,MGCSU),CC(MXCF)
      INTEGER LA(MRU,KAORDS),NIR(KAORDS),NIRL(KAORDS),NPERSH(5)
      INTEGER ND(MSTU),LAMBDA(MXSYM),ICT(MXAT,MXTRN),KLASS(MXAT)
      INTEGER KATOM(MXSH),KTYPE(MXSH),KMIN(MXSH),ISC(MXSH,MXPSYM)
      INTEGER IPC(MXSH,MXPSYM),LOC(MXAT,5),LOC2(MXAT,5),LAB(MXCF)
      INTEGER ITYP(MSTU),SLAB(MXSYM)
      ITM=NT
      NTSAVE=NT
      NT=NIRRED
      NTD=(NT+3)/4
      IGCSU=0
      ISFRU=0
C
      DO 120 I=1,KAORDS
        NIR(I)=0
        NIRL(I)=0
        DO 110 J=1,MRU
          LA(J,I)=0
  110   CONTINUE
  120 CONTINUE
C
      DO 150 I=1,MGCSU
        DO 140 J=1,MCSU
          DO 130 K=1,MCTU
            IGC(K,J,I)=0.0D+00
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
C
      DO 160 I=1,MSTU
        ND(I)=LAMBDA(I)
        ITYP(I)=SLAB(I)
  160 CONTINUE
C
      DO 500 NATOM=1,NAT
C
        DO 170 IT=1,ITM
          IF(ICT(NATOM,IT).GT.NATOM) GO TO 500
  170   CONTINUE
C
cges    write(6,*)'atom',natom
        ICLASS=KLASS(NATOM)
C
        DO 180 I=1,NAT
  180     II(I)=0
C
        DO 190 IT=1,ITM
          JCEN=ICT(NATOM,IT)
  190     II(JCEN)=1
C
        NP=0
        DO 400 KT=1,5
          NP=NP+KT
C
          DO 210 I=1,NSHELL
            ISHELL=I
            J=KATOM(I)
            IF(J.NE.NATOM) GO TO 210
            IF(KTYPE(I).EQ.KT) GO TO 240
            IF(KT.EQ.1.AND.KTYPE(I).EQ.2.AND.KMIN(I).EQ.1) GO TO 220
  210     CONTINUE
          GO TO 400
  220     CONTINUE
          DO 230 IQ=1,NTD
            INDOUT(IQ)=ISC(I,IQ)
  230     CONTINUE
          GO TO 260
  240     CONTINUE
          DO 250 IQ=1,NTD
            INDOUT(IQ)=IPC(I,IQ)
  250     CONTINUE
C
  260     CALL ISOOUT                                 
          KA=LOC(ICLASS,KT)
          KAA=KA
          KB=LOC2(ICLASS,KT)
C
          IGCSU=IGCSU+1
          NSOL=0
          KSO=0
C
  270     I5=LAB(KA+1)
          CALL UNPAK
          IRP=I3
          IPA=I2+INDIN(IRP)-1
          LL=LAMBDA(IRP)
          MAX=NP*LL
C
          KSO=KSO+1
          LA(KSO,IGCSU)=IRP
cges      write(6,*)'kso,igcsu,la(kso,igcsu)',kso,igcsu,la(kso,igcsu)
C
          DO 280 I=1,MAX
  280       CA(NATOM,I)=CC(KA+I)
C
          DO 300 JCEN=1,NAT
            IF(II(JCEN).EQ.0) GO TO 300
            IF(JCEN.EQ.NATOM) GO TO 300
            JCLASS=KLASS(JCEN)
            KC=LOC(JCLASS,KT)+KA-KAA
            DO 290 I=1,MAX
  290         CA(JCEN,I)=CC(KC+I)
  300     CONTINUE
C
          K=0
          DO 330 L=1,LL
               NSOL=NSOL+1
               NAO=0
               DO 320 N=1,NAT
                 IF (II(N).EQ.0) GO TO 320
                 DO 310 I=1,NP
                   NAO=NAO+1
                   IGC(NAO,NSOL,IGCSU)=CA(N,K+I)
  310            CONTINUE
  320          CONTINUE
  330     K=K+NP
C
          KA=KA+MAX
          IF(KA.LT.KB) GO TO 270
C
          NIR(IGCSU)=KSO
          NIRL(IGCSU)=NSOL
cges
c    vtol is a realistic zero to decide if igc is zero or not ...
            vtol=1.0d-10  
cges      write(6,*) 'ka  kb  kc ', ka,kb,kc
cges      write(6,*)'i5',i5
C
C  REMOVE 3S SYMMETRY ORBITALS
C
          IF(KT.EQ.3.AND.NPERSH(3).EQ.1) THEN
cges        write(6,*)' entering NPERSH'
            ICTU=NIR(IGCSU)
            ICS=1
            ICSU=ICTU
C  LOOP OVER ICS
c    
cges 082189
c following has been changed to something more realistic than a hard zero
c 340       IF(IGC(1,ICS,IGCSU).EQ.0.0D0) GO TO 348
c
  340       continue
            val=dabs(igc(1,ics,igcsu))
            if(val.lt.vtol)go to 348
cges        write(6,*)' entering 348 if, ics,igcsu',ics,igcsu
cges        write(6,*)'             igc(1,ics,igcsu) ',igc(1,ics,igcsu)
            DO 342 JCT=1,ICTU
              IGC(JCT,ICS+2,IGCSU)=2*IGC(JCT,ICS+2,IGCSU)-
     1          IGC(JCT,ICS,IGCSU)-IGC(JCT,ICS+1,IGCSU)
              IGC(JCT,ICS,IGCSU)=IGC(JCT,ICS,IGCSU)-IGC(JCT,ICS+1,IGCSU)
  342       CONTINUE
            JCSL=ICS+2
            DO 346 JCS=JCSL,ICSU
            DO 344 JCT=1,ICTU
              IGC(JCT,JCS-1,IGCSU)=IGC(JCT,JCS,IGCSU)
  344       CONTINUE
            LA(JCS-1,IGCSU)=LA(JCS,IGCSU)
  346       CONTINUE
            NIR(IGCSU)=NIR(IGCSU)-1
cges        write(6,*)' substracted 1 from nir'
            ICSU=ICSU-1
            ICS=ICS+1
  348       ICS=ICS+1
            IF(ICS.LE.ICSU) GO TO 340
          END IF
C
C  REMOVE 4P SYMMETRY ORBITALS
C
          IF(KT.EQ.4.AND.NPERSH(4).EQ.1) THEN
            ICTU=NIR(IGCSU)
            ICS=1
            ICSU=ICTU
            ICS7=0
            ICS8=0
            ICS9=0
C  LOOP OVER ICS
  350       continue
            val1=dabs(igc(1,ics,igcsu)) 
            val2=dabs(igc(2,ics,igcsu)) 
            val3=dabs(igc(3,ics,igcsu)) 
            val4=dabs(igc(4,ics,igcsu)) 
            val5=dabs(igc(5,ics,igcsu)) 
            val6=dabs(igc(6,ics,igcsu)) 
            val7=dabs(igc(7,ics,igcsu)) 
            val8=dabs(igc(8,ics,igcsu)) 
            val9=dabs(igc(9,ics,igcsu)) 
            if(val1.gt.vtol)ics1=ics
            if(val2.gt.vtol)ics2=ics
            if(val3.gt.vtol)ics3=ics
            if(val4.gt.vtol)ics4=ics
            if(val5.gt.vtol)ics5=ics
            if(val6.gt.vtol)ics6=ics
            if(val7.gt.vtol)ics7=ics
            if(val8.gt.vtol)ics8=ics
            if(val9.gt.vtol)ics9=ics
c 350       IF(IGC(1,ICS,IGCSU).NE.0.0D0) ICS1=ICS
c           IF(IGC(2,ICS,IGCSU).NE.0.0D0) ICS2=ICS
c           IF(IGC(3,ICS,IGCSU).NE.0.0D0) ICS3=ICS
c           IF(IGC(4,ICS,IGCSU).NE.0.0D0) ICS4=ICS
c           IF(IGC(5,ICS,IGCSU).NE.0.0D0) ICS5=ICS
c           IF(IGC(6,ICS,IGCSU).NE.0.0D0) ICS6=ICS
c           IF(IGC(7,ICS,IGCSU).NE.0.0D0) ICS7=ICS
c           IF(IGC(8,ICS,IGCSU).NE.0.0D0) ICS8=ICS
c           IF(IGC(9,ICS,IGCSU).NE.0.0D0) ICS9=ICS
cold        IF(ICS.LT.10) GO TO 360
            IF(ICS7.EQ.0.OR.ICS8.EQ.0.OR.ICS9.EQ.0) GO TO 360
            DO 352 JCT=1,ICTU
             IGC(JCT,ICS1,IGCSU)=2*IGC(JCT,ICS1,IGCSU)-
     1         3*(IGC(JCT,ICS6,IGCSU)+IGC(JCT,ICS8,IGCSU))
             IGC(JCT,ICS6,IGCSU)=IGC(JCT,ICS6,IGCSU)-IGC(JCT,ICS8,IGCSU)
             IGC(JCT,ICS2,IGCSU)=2*IGC(JCT,ICS2,IGCSU)-
     1         3*(IGC(JCT,ICS4,IGCSU)+IGC(JCT,ICS9,IGCSU))
             IGC(JCT,ICS4,IGCSU)=IGC(JCT,ICS9,IGCSU)-IGC(JCT,ICS4,IGCSU)
             IGC(JCT,ICS3,IGCSU)=2*IGC(JCT,ICS3,IGCSU)-
     1         3*(IGC(JCT,ICS5,IGCSU)+IGC(JCT,ICS7,IGCSU))
             IGC(JCT,ICS5,IGCSU)=IGC(JCT,ICS5,IGCSU)-IGC(JCT,ICS7,IGCSU)
  352       CONTINUE
            KCS(1)=ICS
            KCS(3)=MIN(ICS7,ICS8,ICS9)
            KCS(2)=ICS7+ICS8+ICS9-KCS(1)-KCS(3)
            DO 358 I=1,3
            JCSL=KCS(I)+1
            DO 356 JCS=JCSL,ICSU
            DO 354 JCT=1,ICTU
              IGC(JCT,JCS-1,IGCSU)=IGC(JCT,JCS,IGCSU)
  354       CONTINUE
            LA(JCS-1,IGCSU)=LA(JCS,IGCSU)
  356       CONTINUE
            ICSU=ICSU-1
  358       CONTINUE
            NIR(IGCSU)=NIR(IGCSU)-3
            ICS=ICS-3
            ICS7=0
            ICS8=0
            ICS9=0
  360       ICS=ICS+1
            IF(ICS.LE.ICSU) GO TO 350
          END IF
C
C  REMOVE 5S AND 5D SYMMETRY ORBITALS
C
          IF(KT.EQ.5.AND.NPERSH(5).EQ.1) THEN
            ICTU=NIR(IGCSU)
            ICS=1
            ICSU=ICTU
            ICS12=0
            ICS13=0
            ICS14=0
            ICS15=0
C  LOOP OVER ICS
c 370       IF(IGC(1,ICS,IGCSU).NE.0.0D0) ICS1=ICS
c           IF(IGC(2,ICS,IGCSU).NE.0.0D0) ICS2=ICS
c           IF(IGC(3,ICS,IGCSU).NE.0.0D0) ICS3=ICS
c           IF(IGC(4,ICS,IGCSU).NE.0.0D0) ICS4=ICS
c           IF(IGC(5,ICS,IGCSU).NE.0.0D0) ICS5=ICS
c           IF(IGC(6,ICS,IGCSU).NE.0.0D0) ICS6=ICS
c           IF(IGC(7,ICS,IGCSU).NE.0.0D0) ICS7=ICS
c           IF(IGC(8,ICS,IGCSU).NE.0.0D0) ICS8=ICS
c           IF(IGC(9,ICS,IGCSU).NE.0.0D0) ICS9=ICS
c           IF(IGC(10,ICS,IGCSU).NE.0.0D0) ICS10=ICS
c           IF(IGC(11,ICS,IGCSU).NE.0.0D0) ICS11=ICS
c           IF(IGC(12,ICS,IGCSU).NE.0.0D0) ICS12=ICS
c           IF(IGC(13,ICS,IGCSU).NE.0.0D0) ICS13=ICS
c           IF(IGC(14,ICS,IGCSU).NE.0.0D0) ICS14=ICS
c           IF(IGC(15,ICS,IGCSU).NE.0.0D0) ICS15=ICS
  370       continue
            val1=dabs(igc(1,ics,igcsu)) 
            val2=dabs(igc(2,ics,igcsu)) 
            val3=dabs(igc(3,ics,igcsu)) 
            val4=dabs(igc(4,ics,igcsu)) 
            val5=dabs(igc(5,ics,igcsu)) 
            val6=dabs(igc(6,ics,igcsu)) 
            val7=dabs(igc(7,ics,igcsu)) 
            val8=dabs(igc(8,ics,igcsu)) 
            val9=dabs(igc(9,ics,igcsu)) 
            val10=dabs(igc(10,ics,igcsu)) 
            val11=dabs(igc(11,ics,igcsu)) 
            val12=dabs(igc(12,ics,igcsu)) 
            val13=dabs(igc(13,ics,igcsu)) 
            val14=dabs(igc(14,ics,igcsu)) 
            val15=dabs(igc(15,ics,igcsu)) 
            if(val1.gt.vtol)ics1=ics
            if(val2.gt.vtol)ics2=ics
            if(val3.gt.vtol)ics3=ics
            if(val4.gt.vtol)ics4=ics
            if(val5.gt.vtol)ics5=ics
            if(val6.gt.vtol)ics6=ics
            if(val7.gt.vtol)ics7=ics
            if(val8.gt.vtol)ics8=ics
            if(val9.gt.vtol)ics9=ics
            if(val10.gt.vtol)ics10=ics
            if(val11.gt.vtol)ics11=ics
            if(val12.gt.vtol)ics12=ics
            if(val13.gt.vtol)ics13=ics
            if(val14.gt.vtol)ics14=ics
            if(val15.gt.vtol)ics15=ics
            IF(ICS12.EQ.0.OR.ICS13.EQ.0.OR.ICS14.EQ.0.OR.ICS15.EQ.0)
     1         GO TO 390
C           IF(ICS.LT.15) GO TO 390
            DO 371 JCT=1,ICTU
             IGC(JCT,ICS15,IGCSU)=6*IGC(JCT,ICS15,IGCSU)-
     1         IGC(JCT,ICS4,IGCSU)-IGC(JCT,ICS6,IGCSU)
             IGC(JCT,ICS14,IGCSU)=6*IGC(JCT,ICS14,IGCSU)-
     1         IGC(JCT,ICS5,IGCSU)-IGC(JCT,ICS8,IGCSU)
             IGC(JCT,ICS13,IGCSU)=6*IGC(JCT,ICS13,IGCSU)-
     1         IGC(JCT,ICS7,IGCSU)-IGC(JCT,ICS9,IGCSU)
             IGC(JCT,ICS9,IGCSU)=-IGC(JCT,ICS9,IGCSU)+
     1         IGC(JCT,ICS7,IGCSU)
             IGC(JCT,ICS8,IGCSU)=IGC(JCT,ICS8,IGCSU)-
     1         IGC(JCT,ICS5,IGCSU)
             IGC(JCT,ICS6,IGCSU)=IGC(JCT,ICS4,IGCSU)-
     1         IGC(JCT,ICS6,IGCSU)
             TEMP1=IGC(JCT,ICS1,IGCSU)
             TEMP2=IGC(JCT,ICS2,IGCSU)
             TEMP3=IGC(JCT,ICS3,IGCSU)
             IGC(JCT,ICS3,IGCSU)=TEMP1-TEMP2
     1        -6*IGC(JCT,ICS11,IGCSU)+6*IGC(JCT,ICS12,IGCSU)
             IGC(JCT,ICS2,IGCSU)=-TEMP1-TEMP2+2*TEMP3
     1        +12*IGC(JCT,ICS10,IGCSU)
     2        -6*IGC(JCT,ICS11,IGCSU)-6*IGC(JCT,ICS12,IGCSU)
             IGC(JCT,ICS1,IGCSU)=TEMP1+TEMP2+TEMP3
     1        -3*IGC(JCT,ICS10,IGCSU)
     2        -3*IGC(JCT,ICS11,IGCSU)-3*IGC(JCT,ICS12,IGCSU)
  371       CONTINUE
***********************************************************************
*    NOW SORT OUT THE ORDER OF ICS 4,5,7,10,11,12                     *
***********************************************************************
            KCS(1)=ICS4
            KCS(2)=ICS5
            KCS(3)=ICS7
            KCS(4)=ICS10
            KCS(5)=ICS11
            KCS(6)=ICS12
            DO 373 I=1,6
               DO 372 J=1,I-1
                  IF(KCS(I).GT.KCS(J)) THEN
                     ITEMP=KCS(I)
                     KCS(I)=KCS(J)
                     KCS(J)=ITEMP
                  ENDIF
  372          CONTINUE
  373       CONTINUE
            DO 376 I=1,6
            JCSL=KCS(I)+1
            DO 375 JCS=JCSL,ICSU
            DO 374 JCT=1,ICTU
              IGC(JCT,JCS-1,IGCSU)=IGC(JCT,JCS,IGCSU)
  374       CONTINUE
            LA(JCS-1,IGCSU)=LA(JCS,IGCSU)
c           write(6,*) 'la(i) , i ', la(jcs-1,igcsu),jcs-1
  375       CONTINUE
            ICSU=ICSU-1
  376       CONTINUE
            NIR(IGCSU)=NIR(IGCSU)-6
            ICS=ICS-6
            ICS12=0
            ICS13=0
            ICS14=0
            ICS15=0
  390       ICS=ICS+1
            IF(ICS.LE.ICSU) GO TO 370
          END IF
C
  400   CONTINUE
  500 CONTINUE
      NT=NTSAVE
C
C     ----- PRINT ARRAYS -----
C
      WRITE(ICHECK,900)
  900 FORMAT (//,'    SYMMETRY INFORMATION',//,'   NIR    LA')
      DO 610 JGC=1,MGCSU
        MXR=NIR(JGC)
        WRITE(ICHECK,901) NIR(JGC),(LA(I,JGC),I=1,MXR)
  610 CONTINUE
  901 FORMAT (1X,I3,(T10,24I3))
      WRITE(ICHECK,902)
  902 FORMAT (//,'      IGC MATRICES ',/)
      DO 630 JGC=1,MGCSU
        IRL=NIRL(JGC)
        ICSU=NIR(JGC)
        DO 620 JSO=1,ICSU
          WRITE(ICHECK,903) (IGC(JAO,JSO,JGC),JAO=1,IRL)
  620   CONTINUE
  630 CONTINUE
  903 FORMAT (1X,15F8.3)
      WRITE(ICHECK,904) ND
  904 FORMAT (//,'    ND:',(20I5))
      WRITE(ICHECK,905) ITYP
  905 FORMAT (/,'  ITYP:',(20I5))
CC
      RETURN
      END
