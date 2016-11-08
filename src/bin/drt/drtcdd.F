      BLOCK DATA DRTCDD
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*3 CODES,WORDS
      CHARACTER*1 MULTRF,VALENC
C
      COMMON /DRTCOD/ NCODES,DELA(9),DELB(9),DELELE(9)
     #,               NTYPES,VIRTUL,OCCUPD,VALOCC,RESCOR,RESVIR,FROZEN
     #,               VALVIR,OPENSH,MULTI,SPESHL,MULTRF,VALENC
     #,               WORDS(6,9),CODES(9)
      COMMON /CODE/  FZC, FZV, COR, VIR, DOC, UOC, ALP, BET, SPE
      COMMON /CASES/  CASEV(9)
C
      DATA NCODES /9/
      DATA CODES /'FZC','FZV','COR','VIR','DOC','UOC','ALP','BET','SPE'/
      DATA DELA  /  0  ,  0  ,  1  ,  0  ,  1  ,  0  ,  0  ,  1  ,  0  /
      DATA DELB  /  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1  , -1  ,  0  /
      DATA DELELE/  0  ,  0  ,  2  ,  0  ,  2  ,  0  ,  1  ,  1  ,  0  /
      DATA CASEV /99999,99999,  4  ,  1  ,  4  ,  1  ,  2  ,  3  ,9999 /
      DATA NTYPES /10/
      DATA VIRTUL /1/, OCCUPD/4/, VALOCC/6/, RESCOR/3/, RESVIR/2/
     #,    FROZEN/10/, OPENSH/8/, MULTI /7/, VALVIR/5/, SPESHL/9/
      DATA MULTRF /'/'/
      DATA VALENC /'%'/
      DATA WORDS  /'FRO','ZEN',' CO','RE ','   ','   '
     #,            'FRO','ZEN',' VI','RTU','AL ','   '
     #,            'RES','TRI','CTE','D C','ORE','   '
     #,            'RES','TRI','CTE','D V','IRT','UAL'
     #,            'DOU','BLY',' OC','CUP','IED','   '
     #,            'UNO','CCU','PIE','D  ','   ','   '
     #,            'ALP','HA ','OCC','UPI','ED ','   '
     #,            'BET','A O','CCU','PIE','D  ','   '
     #,            'SPE','CIA','L O','RBI','TAL','   '/
      DATA FZC/1/, FZV/2/, COR/3/, VIR/4/, DOC/5/, UOC/6/, ALP/7/
     #,    BET/8/, SPE/9/
      END
