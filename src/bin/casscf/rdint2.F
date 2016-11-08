      SUBROUTINE RDINT2(BUFO,NDIMB,ITAPE,MAXINT,LBLI,BUFI)
      IMPLICIT REAL*8 (A-H,O-Z)                                         SMM25090
      DIMENSION BUFO(NDIMB),LBLI(MAXBF4),BUFI(MAXBF2)                   SMM25100
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          SMM25110
      COMMON/SIGNS/IOFF(32897),JPRNT
    1 FORMAT(2X,4I6,5X,4I6,5X,F15.8)                                    SMM25130
C                                                                       SMM25140
C   READ IN TWO ELECTRON INTEGRALS                                      SMM25150
      CALL RFILE(ITAPE)                                                 SMM25160
      CALL ZERO(BUFO,MAXINT)                                            SMM25170
      INTI=0                                                            SMM25180
      NBLI=0                                                            SMM25190
  200 IBL=0                                                             SMM25200
      NBLI=NBLI+1                                                       SMM25210
      CALL SREAD(ITAPE,LBLI,MAXBF4)                                     SMM25220
  201 IBL=IBL+1                                                         SMM25230
      JA=LBLI(IBL+IBL-1)                                                SMM25240
      IF(JA.EQ.0) GO TO 205                                             SMM25250
      JAA=LBLI(IBL+IBL)                                                 SMM25260
      CALL UNPACK(I,J,K,L,IX,JA,JAA)                                    SMM25270
      IJ=IOFF(I)+J                                                      SMM25300
      KL=IOFF(K)+L                                                      SMM25310
C     IJKL=IJ*(IJ-1)/2+KL                                               SMM25320
      IJKL=IOFF(IJ)+KL                                                  SMM25330
      BUFO(IJKL)=BUFI(MAXBUF+IBL)                                       SMM25340
C*    WRITE(6,1) I,J,K,L,IJ,KL,IX,IJKL,BUFO(IJKL)                       SMM25350
  202 CONTINUE                                                          SMM25360
      IF(IBL.EQ.MAXBUF) GO TO 203                                       SMM25370
      GO TO 201                                                         SMM25380
  203 INTI=INTI+MAXBUF                                                  SMM25390
      GO TO 200                                                         SMM25400
  205 INTI=INTI+IBL                                                     SMM25410
C                                                                       SMM25420
      CALL RCLOSE(ITAPE,3)                                              SMM25430
      RETURN                                                            SMM25440
      END                                                               SMM25450
