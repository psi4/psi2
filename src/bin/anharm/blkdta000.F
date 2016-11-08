      BLOCK DATA
      CHARACTER*1  RAXIS
      CHARACTER*25 RTYPE
      COMMON/VIB206/RAXIS(3),RTYPE(6)
C
      DATA RAXIS(1) / 'A' /
      DATA RAXIS(2) / 'B' /
      DATA RAXIS(3) / 'C' /
      DATA RTYPE(1) / ' DIATOMIC                ' /
      DATA RTYPE(2) / ' LINEAR POLYATOMIC       ' /
      DATA RTYPE(3) / ' ASYMMETRIC TOP          ' /
      DATA RTYPE(4) / ' SYMMETRIC TOP (PROLATE) ' /
      DATA RTYPE(5) / ' SYMMETRIC TOP (OBLATE)  ' /
      DATA RTYPE(6) / ' SPHERICAL TOP           ' /
      end
