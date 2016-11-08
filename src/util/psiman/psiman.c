
#include <stdio.h>
#include <ctype.h>

#define MAX_STRING 256

main(argc,argv)
int argc;
char **argv;
{
  int i;
  int print=0;
  char *os,*level;
  int advanced;
  char name[MAX_STRING];
  char manfile[MAX_STRING];
  int c;
  int use_section = 0;

  name[0] = '\0';

  /* Defaults. */
  level = "beg";
  os = "unix";

  for (i=1; i<argc; i++) {
    if (!strncmp(argv[i],"-u",2)) {
        os = "unix";
        }
    else if (!strncmp(argv[i],"-c",2)) {
        os = "cms";
        }
    else if (!strncmp(argv[i],"-a",2)) {
        level = "adv";
        }
    else if (!strncmp(argv[i],"-b",2)) {
        level = "beg";
        }
    else if (!strncmp(argv[i],"-p",2)) {
        print = 1;
        }
    else if (  (!strncmp(argv[i],"-?",2))
             ||(!strncmp(argv[i],"-h",2))) {
      psiman_usage();
      }
    else if (!strncmp(argv[i],"-",1)) {
      printf("psiman: unknown option: %s\n",argv[i]);
      psiman_usage();
      }
    else if ((strlen(argv[i])==1)&&isdigit(*argv[i])) {
      if (use_section) {
        printf("psiman: too many manual sections, give only one\n");
        psiman_usage();
        }
      use_section = atoi(argv[i]);
      }
    else {
      /* strcat(names," ");
         strcat(names,argv[i]);
       */
      if (name[0] != '\0') {
        printf("psiman: too many names, give only one\n");
        psiman_usage();
        }
      strcpy(name,argv[i]);
      }
    }

  if (name[0] == '\0') {
    psiman_usage();
    }

  if (chdir(PSI_MANDIR) == -1) {
    perror("psiman: couldn't chdir into the Psi manual directory:");
    exit(1);
    }

  if (use_section) {
    sprintf(manfile,"man%d/%s.%d",use_section,name,use_section);
    if (fexists(manfile)) {
      return print_man(manfile,print,os,level);
      }
    else {
      printf("psiman: couldn't find a manual page for \"%s\" in section %d\n",
             name,
             use_section);
      }
    }
  else {
    for (use_section=1; use_section<=8; use_section++) {
      sprintf(manfile,"man%d/%s.%d",use_section,name,use_section);
      if (fexists(manfile)) {
        return print_man(manfile,print,os,level);
        }
      }
    printf("psiman: couldn't find a manual page for \"%s\"\n", name);
    }

  return(1);
  }

fexists(file)
char *file;
{
  FILE *tmp;

  if (tmp = fopen(file,"r")) {
    fclose(tmp);
    return(1);
    }
  else {
    return(0);
    }
  }


print_man(manfile,print,os,level)
char *manfile;
int  print;
char *os;
char *level;
{
  char command[MAX_STRING];

  if (print) {
    sprintf(command,"cat macro.%s.%s %s |soelim|eqn|psroff -t|rsh %s lpr -h",
            os,level,manfile,HOST);
    }
  else {
    sprintf(command,"cat macro.%s.%s %s |soelim|neqn|nroff|%s",
            os,level,manfile,MORE);

    }
  return system(command);
  }

psiman_usage()
{
  printf("psiman usage: `psiman options [section] name',\n");
  printf("where options are:\n");
  printf("  -h print help info\n");
  printf("  -? same as -h\n");
  printf("  -u display the unix version of the manual\n");
  printf("  -c display the cms version of the manual\n");
  printf("  -a display the advanced version of the manual\n");
  printf("  -b display the beginners version of the manual\n");
  printf("  -p print the manual page using ptroff\n");
  printf("where section is:\n");
  printf("   1 for the general Psi commands\n");
  printf("   2 for the inputs sections common to several Psi commands\n");
  printf("   3 procedures\n");
  printf("   4 examples for the use of Psi\n");
  printf("   if section is omitted, all sections are searched\n");
  printf("where name is:\n");
  printf("   the name of the manual page\n");
  exit(1);
  }
