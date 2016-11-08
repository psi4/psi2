


C///////////////////////////////////////////////////////////////////////
      SUBROUTINE RCLOSE(ITAPE,JCODE)	
      IMPLICIT INTEGER (A-Z)	

C JCODE = 4     CLOSE AND DELETE FILE	
C JCODE = 3     CLOSE AND SAVE FILE	
	
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'	

	IF (JCODE.EQ.4) THEN
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='DELETE')
                      ELSE
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='KEEP')
             END IF
		
      IF (IERR.NE.0) THEN
             WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE	
             WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE	
             END IF
 	
      RETURN	
      END
