      BLOCK DATA
      COMMON/IJVEC/IX(35),IY(35),IZ(35),JX(35),JY(35),JZ(35)
C
C                                         X  Y  Z  X  X  Y  Y  Z  Z  X
C                       X  Y  Z  X  X  Y  X  Y  Z  X  X  Y  Y  Z  Z  Y
C              S X Y Z  X  Y  Z  Y  Z  Z  X  Y  Z  Y  Z  X  Z  X  Y  Z
C
ctph  DATA IX /1,5,1,1, 9, 1, 1, 5, 5, 1,13, 1, 1, 9, 9, 5, 1, 5, 1, 5/
ctph  DATA IY /1,1,5,1, 1, 9, 1, 5, 1, 5, 1,13, 1, 5, 1, 9, 9, 1, 5, 5/
ctph  DATA IZ /1,1,1,5, 1, 1, 9, 1, 5, 5, 1, 1,13, 1, 5, 1, 5, 9, 9, 5/
ctph  DATA JX /0,1,0,0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
ctph  DATA JY /0,0,1,0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
ctph  DATA JZ /0,0,0,1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IX /1,6,1,1,11, 1, 1, 6, 6, 1,16, 1, 1,11,11, 6, 1, 6, 1, 6,
     1         21, 1, 1,16,16, 6, 1, 6, 1,11,11, 1,11, 6, 6/
      DATA IY /1,1,6,1, 1,11, 1, 6, 1, 6, 1,16, 1, 6, 1,11,11, 1, 6, 6,
     1          1,21, 1, 6, 1,16,16, 1, 6,11, 1,11, 6,11, 6/
      DATA IZ /1,1,1,6, 1, 1,11, 1, 6, 6, 1, 1,16, 1, 6, 1, 6,11,11, 6,
     1          1, 1,21, 1, 6, 1, 6,16,16, 1,11,11, 6, 6,11/
      DATA JX /0,1,0,0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     1         4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      DATA JY /0,0,1,0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     1         0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      DATA JZ /0,0,0,1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     1         0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
C
      END
