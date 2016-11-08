      block data blkatm
      character*2 symb
      character*40 atoms
      common/atoml/symb(103),atoms(103)
c
cets someone else can finish this
c
      data (atoms(i),i=1,4)/'HYDROGEN','HELIUM','LITHIUM','BERYLLIUM'/
      data (atoms(i),i=5,8)/'BORON','CARBON','NITROGEN','OXYGEN'/
      data (atoms(i),i=9,12)/'FLUORINE','NEON','SODIUM','MAGNESIUM'/
      data (atoms(i),i=13,15)/'ALUMINUM','SILICON','PHOSPHORUS'/
      data (atoms(i),i=16,18)/'SULFUR','CHLORINE','ARGON'/
      data (atoms(i),i=19,21)/'POTASSIUM','CALCIUM','SCANDIUM'/
      data (atoms(i),i=22,24)/'TITANIUM','VANADIUM','CHROMIUM'/
      data (atoms(i),i=25,28)/'MANGANESE','IRON','COBALT','NICKEL'/
      data (atoms(i),i=29,32)/'COPPER','ZINC','GALLIUM','GERMANIUM'/
      data (atoms(i),i=33,35)/'ARSENIC','SELENIUM','BROMINE'/
      data (atoms(i),i=36,38)/'KRYPTON','RUBIDIUM','STRONTIUM'/
      data (atoms(i),i=39,41)/'YTTRIUM','ZIRCONIUM','NIOBIUM'/
      data (atoms(i),i=42,44)/'MOLYBDENUM','TECHNETIUM','RUTHENIUM'/
      data (atoms(i),i=45,47)/'RHODIUM','PALLADIUM','SILVER'/
      data (atoms(i),i=48,51)/'CADMIUM','INDIUM','TIN','ANTIMONY'/
      data (atoms(i),i=52,55)/'TELLURIUM','IODINE','XENON','CESIUM'/
      data (atoms(i),i=56,57)/'BARIUM','LANTHANUM'/
      data (symb(i),i=1,9)/'H','HE','LI','BE','B','C','N','O','F'/
      data (symb(i),i=10,17)/'NE','NA','MG','AL','SI','P','S','CL'/
      data (symb(i),i=18,25)/'AR','K','CA','SC','TI','V','CR','MN'/
      data (symb(i),i=26,33)/'FE','CO','NI','CU','ZN','GA','GE','AS'/
      data (symb(i),i=34,41)/'SE','BR','KR','RB','SR','Y','ZR','NB'/
      data (symb(i),i=42,49)/'MO','TC','RU','RH','PD','AG','CD','IN'/
      data (symb(i),i=50,57)/'SN','SB','TE','I','XE','CS','BA','LA'/
c
      end
