
#include <stdio.h>

#if defined(AIX)
gfname(status)
#else
gfname_(status)
#endif
int *status;
{
  char buffer[1024];

  if (get_param("FILES:name","%s",buffer) == -1) {
    *status = 0;
    return;
    }

  *status = 1;
  fprintf(stdout,"%s",buffer);
  }
