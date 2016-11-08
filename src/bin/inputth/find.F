      SUBROUTINE FIND(GET,ZNUC,IR,LINENO)
C
      REAL*8 ZNUC
      CHARACTER*80 GET,LINE
      CHARACTER*12 ATNAME(92)
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
C
      DATA ATNAME /'HYDROGEN    ','HELIUM      ','LITHIUM     '
     #,            'BERYLLIUM   ','BORON       ','CARBON      '
     #,            'NITROGEN    ','OXYGEN      ','FLUORINE    '
     #,            'NEON        ','SODIUM      ','MAGNESIUM   '
     #,            'ALUMINUM    ','SILICON     ','PHOSPHORUS  '
     #,            'SULFUR      ','CHLORINE    ','ARGON       '
     #,            'POTASSIUM   ','CALCIUM     ','SCANDIUM    '
     #,            'TITANIUM    ','VANADIUM    ','CHROMIUM    '
     #,            'MANGANESE   ','IRON        ','COBALT      '
     #,            'NICKEL      ','NICKEL      ','ZINC        '
     #,            'GALLIUM     ','GERMANIUM   ','ARSENIC     '
     #,            'SELENIUM    ','BROMINE     ','KRYPTON     '
     #,            56*'123456789012'                            /
C
      write(6,*)get(5:64)
      LINENO=0
      REWIND IR
      DO 10 I=1,10000
         READ (IR,900,END=9000,ERR=9000) LINE
         LINENO=LINENO+1
  900    FORMAT (A80)
         IF (LINE(1:60).EQ.GET(5:64)) GO TO 20
   10 CONTINUE
      GO TO 9000
C
   20 CONTINUE
      NUM=ZNUC+0.5D+00
      DO 30 I=1,10000
         READ (IR,900,END=9000,ERR=9000) LINE
         LINENO=LINENO+1
         IF (LINE(1:12).EQ.ATNAME(NUM)) RETURN
   30 CONTINUE
C
 9000 CONTINUE
      WRITE(IOUT,901) GET
  901 FORMAT (//,' INPUT (FIND) : CANNOT FIND LINE:',1X,A80,//)
      STOP
      END
