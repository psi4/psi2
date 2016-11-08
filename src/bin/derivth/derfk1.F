      SUBROUTINE DERFK1(D1TOT,D2TOT,DFH,DFS,BMAT,SS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer prcntr, prdup, errcod, ip, prup
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2
      COMMON/CALCS/ISTYPE,ICTYPE,IRAMAN,IAFOCK,IDFOCK,IDINT
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/MFSEC/MFILE,NSECT
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/TOLER/DLIM
      DIMENSION D1TOT(3,NATOM),D2TOT(N3N,N3N)
      DIMENSION DFH(NBATRI,N3N),DFS(NBATRI,3),BMAT(NBATRI,N3N,NTREAD)
      DIMENSION SS(NBATRI)
      DIMENSION IBF2AT(256)
      DATA ZERO / 0.0D+00 /
C
C   LOCATION OF MATRICES IN TAPE42
C   (1) E1A
C   (2) E2A
C   (3) SA (N3NXX*NTRIL)
C   (4) HA (N3NXX*NTRIL)
C   (5) TA (N3NXX*NTRIL*NBSET)
C
C
      CALL BFTOAT(IBF2AT)
      ITAP42=42
      CALL RFILE(ITAP42)
C
C12-1-88  Write IBF2AT to file42 first for CI gradients since NGRCPHF
C12-1-88  expects it.
C
      IF(ICTYPE.EQ.2) THEN
c12-1-88CALL SWRIT(ITAP42,IBF2AT,NBFAO*2)
        CALL SWRIT(ITAP42,IBF2AT,NBASIS*2)
      ENDIF
C
      CALL SWRIT(ITAP42,D1TOT,N3N*2)
      CALL SWRIT(ITAP42,D2TOT,N3N*N3N*2)
C
c   
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        do 301 i=1,3
          write(6,*) 'in derfk1:  dfs i=',i
          call print(dfs(1,i),nbatri,nbfao,6)
 301    continue
      endif
C   SAVE SA MATRICES ON TAPE42
cets030691 ci cphf codes need dfs rather than ss
      if(ictype.eq.2) then
        do 302 ixyz=1,3
          call swrit(itap42,dfs(1,ixyz),intowp(nbatri))
 302    continue
      else
        ip = prcntr('IS_ON DEBUG')
        KABC=0
        DO 106 IABC=1,NATOM
          DO 105 IXYZ=1,3
            KABC=KABC+1
            DO 102 I=1,NBATRI
  102       SS(I)=ZERO
            DO 104 I=1,NBFAO
              IA=IOFF(I)
              DO 103 J=1,I
                IF(IBF2AT(I).NE.IABC.AND.IBF2AT(J).NE.IABC) GO TO 103
                IJ=IA+J
                SX=DFS(IJ,IXYZ)
                IF(IBF2AT(I).NE.IABC) SX=-SX
                SS(IJ)=SX
  103         CONTINUE
  104       CONTINUE
c
            if(ip.eq.1) then
              write(6,*) 'in derfk1:  ss iabc=',iabc,' ixyz=',ixyz
              CALL PRINT(SS,NBATRI,NBFAO,6)
            endif
c
            CALL SWRIT(ITAP42,SS,NBATR2)
  105     CONTINUE
  106   CONTINUE
      endif
C
C   SAVE HA MATRICES ON TAPE42
      ip = prcntr('IS_ON DEBUG')
      DO 107 IABC=1,N3N
c
        if(ip.eq.1) then
          write(6,*) 'in derfk1:  dfh iabc=',iabc
          CALL PRINT(DFH(1,IABC),NBATRI,NBASIS,6)
        endif
c
        CALL SWRIT(ITAP42,DFH(1,IABC),NBATR2)
  107 CONTINUE
C
C   SAVE TA MATRICES ON TAPE42
      ip = prcntr('IS_ON DEBUG')
      DO 109 ITYP=1,NTREAD
        DO 108 IABC=1,N3N
c
          if(ip.eq.1) then
            write(6,*) 'in derfk1: bmat matrix,ityp,iabc = ',ityp,iabc
            call print(bmat(1,iabc,ityp),nbatri,nbasis,6)
          endif
c
          CALL SWRIT(ITAP42,BMAT(1,IABC,ITYP),NBATR2)
  108   CONTINUE
  109 CONTINUE
C
      CALL RCLOSE(ITAP42,3)
      RETURN
      END
