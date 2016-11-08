
C///////////////////////////////////////////////////////////////////////
      subroutine nounfl
C Set cpu control register so no floating-point underflow exceptions

      include '/usr/include/fpdc.h'

	fpstat(fpux) = .FALSE.
	call fpsets (fpstat)

	return
	end
