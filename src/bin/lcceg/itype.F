C
C ----------------------------------------------------------------------
C
      FUNCTION ITYPE(I,J,K,L)
      IF (J-K)150,10,80
10    IF (I-J) 50,20,50
20    IF (K-L) 40,30,40
30    ITYPE=1
      RETURN
40    ITYPE=2
      RETURN
50    IF (K-L) 70,60,70
60    ITYPE=3
      RETURN
70    ITYPE=8
      RETURN
80    IF (I-J) 120,90,120
90    IF (K-L) 110,100,110
100   ITYPE=4
      RETURN
110   ITYPE=6
      RETURN
120   IF (K-L) 140,130,140
130   ITYPE=10
      RETURN
140   ITYPE=12
      RETURN
150   IF (J-L) 220,190,160
160   IF (I-K) 180,170,180
170   ITYPE=7
      RETURN
180   ITYPE=13
      RETURN
190   IF (I-K) 210,200,210
200   ITYPE=5
      RETURN
210   ITYPE=9
      RETURN
220   IF (K-L) 240,230,240
230   ITYPE=11
      RETURN
240   ITYPE=14
      RETURN
      END
