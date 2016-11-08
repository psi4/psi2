      BLOCK DATA
      CHARACTER*3 ISMAT
      COMMON/SDATA/ISMAT(8,8)
      DATA ISMAT(1,1)            / '  A'         /
      DATA ISMAT(2,1),ISMAT(2,2) / ' AP' , 'APP' /
      DATA ISMAT(3,1),ISMAT(3,2) / ' AG' , ' AU' /
      DATA ISMAT(4,1),ISMAT(4,2) / '  A' , '  B' /
      DATA ISMAT(5,1),ISMAT(5,2) / ' AG' , ' AU' /
      DATA ISMAT(5,3),ISMAT(5,4) / ' BG' , ' BU' /
      DATA ISMAT(6,1),ISMAT(6,2) / ' A1' , ' A2' /
      DATA ISMAT(6,3),ISMAT(6,4) / ' B1' , ' B2' /
      DATA ISMAT(7,1),ISMAT(7,2) / '  A' , ' B1' /
      DATA ISMAT(7,3),ISMAT(7,4) / ' B2' , ' B3' /
      DATA ISMAT(8,1),ISMAT(8,2) / ' AG' , 'B1G' /
      DATA ISMAT(8,3),ISMAT(8,4) / 'B2G' , 'B3G' /
      DATA ISMAT(8,5),ISMAT(8,6) / ' AU' , 'B1U' /
      DATA ISMAT(8,7),ISMAT(8,8) / 'B2U' , 'B3U' /
      END
