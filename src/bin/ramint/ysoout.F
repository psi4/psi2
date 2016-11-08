      SUBROUTINE YSOOUT
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ISPAC/INDIN(120),INDOUT(30)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      integer rshift
      iword=0
      do 20 j=1,nt
         if (mod(j-1,4).eq.0) then
            iword=iword+1
            indin(j)=and(indout(iword),255)
         else
            indin(j)=and(rshift(indout(iword),8*mod(j-1,4)),255)
         end if
   20 continue
c
      RETURN                                                            
C                                                                       
C--YSOIN-------------------------------------------------YSOIN          
C                                                                       
      ENTRY YSOIN                                                       
C                                                                       
      iword=0
      DO 40 J=1,NT                                                      
         if (mod(j-1,4).eq.0) then
            iword=iword+1
            indout(iword)=indin(j)
         else
            indout(iword)=or(indout(iword),lshift(indin(j),
     #                       8*mod(j-1,4)))
         end if
   40 CONTINUE                                                          
c
      RETURN
      END
