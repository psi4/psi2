#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_label():  Reads the main file30 label.
**
**   takes no arguments.
**
**   returns: char *
*/

char *file30_rd_label(void)
{
  int i;
  char *label;
  label = (char *) malloc(20*sizeof(char));
  for(i=0; i < 20; i++) label[i] = '\0';
  strncpy(label,info30_.label,20);
  return label;
}
