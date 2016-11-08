      SUBROUTINE SIMULT(GMO,RMO,CC,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM20120
      REAL*16 DLIM,DETERM           
      DIMENSION GMO(NIND),RMO(NIND)
      DIMENSION CC(NIND,NIND+1)
      dimension lbli(maxbf4),bufi(maxbf2)
      COMMON/BASIS/NBASIS,NTRI,NTOT,NSYM
      COMMON/CLIMT/cmlim,dlamb
      common/limit/maxbas,maxbuf,maxbf2,maxbf4
      COMMON/SIGNS/IOFF(32897),IPRNT 
      common/works/jtape1,jtape2
      COMMON/CI102/NIND,NINTRI
      COMMON/CI103/MXMIC,IDLIM
      data a00 / 0.0d+00 /
    2 FORMAT(//,2X,' THE CC MATRIX IN SIMULT'/)
    3 FORMAT(//,2X,' DETERM IN SIMULT = ',Q20.10,5x,' IP100 = ',i5,
     1             5x,' IM100 = ',i5/)
    4 format(//,2x,' The solutions after flinq'/) 
C                                                                       SMM20320
c   Set the smallest positive determinant value 
c      allowed before going to the level shift
      dlim= 1.0q+00 / ( 10.0q+00**(idlim) )
C   FORM THE CC MATRIX                                                  SMM20430
      MICRO=-1                                                          SMM20470
  300 CONTINUE                                                          SMM20480
      MICRO=MICRO+1                                                     SMM20490
      DINCR=dlamb*dfloat(MICRO)
      write(3,*) '  micro = ',micro
      write(3,*) '  dincr = ',dincr
      CALL ZERO(CC,nind*(nind+1))
C                                                                       SMM20520
C   Read in the A11 matrix elements
      call rfile(jtape1)
  500 continue
      call sread(jtape1,lbli,maxbf4)
      do 101 m=1,maxbuf
      ii=lbli(m+m-1)
      if(ii.eq.0) go to 600
      jj=lbli(m+m)
      cc(ii,jj)=bufi(maxbuf+m)
      if(ii.eq.jj) go to 301
      cc(jj,ii)=bufi(maxbuf+m)
      GO TO 101                                                         SMM20610
  301 CONTINUE                                                          SMM20620
      cc(ii,ii)=cc(ii,ii)+dincr
  101 CONTINUE                                                          SMM20640
      go to 500
  600 continue
      call rclose(jtape1,3)
C                                                                       SMM20650
C   THE GMO TERMS                                                       SMM20860
      DO 102 I=1,NIND                                                   SMM20870
      CC(I,nind+1)= - GMO(I)   
  102 CONTINUE                                                          SMM20900
C                                                                       SMM20910
      if(iprnt.gt.3) then
        WRITE(6,2)                                           
        CALL MATOUT(CC,nind,nind+1,nind,nind+1,6)
      end if
C                                                                       SMM21000
C   SOLVE THE SIMULTANEOUS EQUATIONS                                    SMM21010
cyy   CALL FLINQ(CC,nind,nind,1,DETERM)
      call qflin(cc,nind,nind,1,determ,ip100,im100)
C                                                                       SMM21040
C   TRANSFER ANSWERES TO THE RMO ARRAY                                  SMM21050
      DO 103 I=1,nind
      RMO(I)=CC(I,nind+1)
  103 CONTINUE                                                          SMM21080
      AMOMAX=A00                                                        SMM21090
      DO 104 I=1,NIND                                                   SMM21100
      TEMP=DABS(RMO(I))                                            
      IF(TEMP.LE.AMOMAX) GO TO 104                                      SMM21120
      AMOMAX=TEMP                                                       SMM21130
  104 CONTINUE                                                          SMM21140
C*    WRITE(6,*) ' DETERM,AMOMAX = ',DETERM,AMOMAX
      WRITE(4,*) ' DETERM,AMOMAX = ',DETERM,AMOMAX
c
c   Determ is a ``negative'' value
      IF(DETERM.LT.0.0Q+00) THEN                                        SMM21220
        WRITE(6,*) '  MICRO = ',MICRO
        write(6,3) determ,ip100,im100
        IF(MICRO.GE.MXMIC) GO TO 302                                    SMM21250
        GO TO 300                                                       SMM21260
      END IF                                                            SMM21270
c
c   Determ is a small positive value
      if(qabs(determ).lt.Dlim) then
        WRITE(6,*) '  MICRO = ',MICRO
        write(6,3) determ,ip100,im100
        IF(MICRO.GE.MXMIC) GO TO 302                                    SMM21250
        GO TO 300                                                       SMM21260
      end if
c
        IF(AMOMAX.LE.cmlim) GO TO 302     
C                                                                       SMM21290
C   SCALE THE MO CHANGE IF REQUIRED                                     SMM21300
      CALL ANORM(DMO,RMO,NIND,NIND)                                     SMM21310
      SCALE=cmlim/DMO      
      DO 105 I=1,NIND                                                   SMM21330
      RMO(I)=RMO(I)*SCALE                                               SMM21340
  105 CONTINUE                                                          SMM21350
C                                                                       SMM21360
  302 CONTINUE                                                          SMM21380
      write(6,3) determ,ip100,im100
      if(Iprnt.gt.3) then
        WRITE(6,4)                                     
        CALL MATOUT(rmo,nind,1,nind,1,6)
      end if
C                                                                       SMM21520
      IF(DETERM.LT.0.0Q+00) THEN
        write(6,*) '  '
        write(6,*) '  '
        write(6,*) '  ERROR:  You have a negative determinant.'
        write(6,*) '  '
        write(6,*) '  To remedy the situation increase the number'
        write(6,*) '    of MICROITER and/or increase the HESSHIFT.'
        write(6,*) '  '
        write(6,*) '  If such a remedy does not work, then please'
        write(6,*) '    reexamine for CASSCF space and/or find'
        write(6,*) '    someone who can help you.'
        write(6,*) '    (e.g. preferably someone who has a clue)'
        write(6,*) '  '
        write(6,*) '  '
        call pabort
      endif

      RETURN                                                            SMM21540
      END                                                               SMM21550
