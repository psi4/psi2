/*-----------------------------------------------------------------------

  guess.c: Function that reads the guessing parameters from input and then 
           either uses them to form an initial guess or calculates an initial 
           guess from multiplicity and the charge using a diagonalization 
           of the core guess.

  Modified 7/5/99 by C. David Sherrill to make the user input override
           any default behavior: PSI should assume the user is smarter 
           than the program.
------------------------------------------------------------------------*/
#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>
struct symm *s;
int i,j,bool,k,m,l,o;
int errcod;
int nc,no,nh,nn;
int size;
int optri;
char *guess_opt,*occ_fix_str;
int  *cldpi,*openpi;

void guess()
{
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");

/*--- Initialize occupation arrays ----*/
  for(i=0; i < num_ir; i++) {
      scf_info[i].nclosed=0;
      scf_info[i].nopen=0;
    }

/*---Is there user input for occupations?  If so, do what we're told---*/
  if (ip_exist("DOCC",0)) {
      fprintf(outfile, "\nUsing DOCC and SOCC to determine occupations\n\n");
      errcod = ip_count("DOCC",&size,0);
      if(errcod == IPE_OK && size != num_ir) {
          fprintf(outfile,"\n DOCC array is the wrong size\n");
	  fprintf(outfile," is %d, should be %d\n",size,num_ir);
	  exit(size);
	}
      if(errcod != IPE_OK && !iopen) {
	  fprintf(outfile,"\n try adding some electrons buddy!\n");
	  fprintf(outfile," need DOCC\n");
	  exit(1);
	}
	      
      if(iopen) {
	  errcod = ip_count("SOCC",&size,0);
	  if(errcod == IPE_OK && size != num_ir) {
	      fprintf(outfile,"\n SOCC array is the wrong size\n");
	      fprintf(outfile," is %d, should be %d\n",size,num_ir);
	      exit(size);
	    }
		  
	  errcod = ip_count("HOCC",&size,0);
	  if(errcod == IPE_OK && size != num_ir) {
	      fprintf(outfile,"\n HOCC array is the wrong size\n");
	      fprintf(outfile," is %d, should be %d\n",size,num_ir);
	      exit(size);
	    }
	}
      n_open = 0;
	      
      for (i=0; i < num_ir; i++) {
	  errcod = ip_data("DOCC","%d",&scf_info[i].nclosed,1,i);
	  if(iopen) errcod=ip_data("SOCC","%d",&scf_info[i].nopen,1,i);
	  if(iopen) errcod=ip_data("HOCC","%d",&scf_info[i].nhalf,1,i);
        }
    } 

/*---Is there a prior calculation? if so read occupations from file----*/
  else if(inflg ==  1){
      fprintf(outfile,"\nReading Occupations from file30\n\n");
      cldpi = file30_rd_clsdpi();
      openpi = file30_rd_openpi();
	  
      for(o = 0;o < num_ir;o++){
          scf_info[o].nclosed=cldpi[o];
          scf_info[o].nopen=openpi[o];
        }
    }

/*---There is no DOCC,SOCC or a prior calculation, need to make a guess---*/
  else{

      errcod = ip_string("GUESS",&guess_opt,0);
      /* CDS: need to set the default if it's not read, else will seg fault
       * below because there won't be any guess_opt (will pt to NULL)
       * remember to change the malloc size if you change the default string
       */
      if (errcod == IPE_KEY_NOT_FOUND) {
          guess_opt = (char *) malloc(sizeof(char)*5);
          strcpy(guess_opt, "CORE");
        }

      errcod = ip_data("MULTP","%d",&multp,0);
      if(errcod != IPE_OK)
	  multp = 1;
      if(multp > 1 && hsos != 1) {
	  fprintf(outfile,
                         "\n\n\n******* Must Specify Opentype = highspin");
	  fprintf(outfile," for multiplicity: %d\n\n",multp);
	  exit(1);
	}
      if(multp == 1 && hsos == 1) {
	  fprintf(outfile,"\n\n\n****** Please Set Opentype = singlet");
	  fprintf(outfile," for open shell singlets\n\n");
	  exit(1);
	}
	      
      occ_fix = 0;
      errcod = ip_boolean("OCC_FIX",&occ_fix,0);
      
  /* read in the charge */
	  
      errcod = ip_data("CHARGE","%d",&charge,0);
      if(errcod != IPE_OK)
	  charge = 0;
	      
  /* caclulate the number of electrons */
	      	      
      if(!strcmp(guess_opt,"HUCKEL")){
	  if(iter == 0){
	      fprintf(outfile, "\n\nHUCKEL not implemented yet");
	      fprintf(outfile, "\nDefaulting to a core guess\n\n");
	    }
	}
      if(!strcmp(guess_opt,"CORE")){
	  if(iter == 0){
	      fprintf(outfile, "\n\nUsing Core Guess to");
	      fprintf(outfile, " determine occupations\n\n");
	    }
	}
      else{
	  if(iter == 0){
	      fprintf(outfile, "\n\nUsing Core Guess to");
	      fprintf(outfile," determine occupations\n\n");  
	    }
	}
	        
  /* sort orbitals into two arrays, one with the symmetry label and one 
     with the eigenvalues */
	  
      sortev();
	 
  /* Guess the orbital occupations */

      nelec = 0;
      for(i=0; i < natom;i++)
          nelec = nelec + zvals[i];
      nelec = nelec - charge;

  /* Multiplicity is closed shell singlet */
	  
      if(multp == 1 && iopen == 0) {
          if(nelec%2==1) {
              fprintf(outfile,"\nImpossible multiplicity with charge\n");
	      fprintf(outfile,"and # of electrons specified\n");
	      exit(1); 
            }
	  for(m = 0;m < (nelec/2);m++) 
	      scf_info[symm_tot[m]].nclosed++; 
        }

  /* Multiplicity is either open shell singlet or triplet */
	  
      else if(multp == 3 || (multp == 1 && singlet == 1)) {
          if(nelec%2==1) {
              fprintf(outfile,"\nImpossible multiplicity with charge\n");
              fprintf(outfile,"and # of electrons specified\n");
              exit(1); 
            }
          for(m = 0;m < ((nelec/2)-1);m++)
              scf_info[symm_tot[m]].nclosed++;
              for(k = (nelec/2)-1;k < (nelec/2)+1;k++)
              scf_info[symm_tot[k]].nopen++; 
        }

  /* Multiplicity is doublet */
	  
      else if(multp == 2) {
          if(nelec%2==0) {
              fprintf(outfile,"\nImpossible multiplicity with charge\n");
              fprintf(outfile,"and # of electrons specified\n");
              exit(2); 
            }
          for(m = 0;m < (nelec-1)/2;m++)
              scf_info[symm_tot[m]].nclosed++;
          scf_info[symm_tot[((nelec-1)/2)]].nopen++; 
        }

  /* Implement more if you want */
	  
      else {
          fprintf(outfile,"\nUnrecognized or unsupported multiplicity\n");
          exit(1);
        }
  }

/*---We're done getting the occupations, now do some bookkeeping---*/

  /* output occupations to outfile */
  if(iter == 1 && occ_fix == 1){
      fprintf(outfile,"\n\n*****Warning: Occupations will be recalculated");
      fprintf(outfile," at every SCF iterations*****\n\n");
  }

  if((iter == 0) || (print & 2 && occ_fix == 1)){
      fprintf(outfile,"\n");  
      fprintf(outfile,"\n  Symmetry block:  ");
      for(l=0;l<num_ir;l++) {
          fprintf(outfile,"%4s  ",scf_info[l].irrep_label);
        }
      fprintf(outfile,"\n  DOCC:            ");
      for(l=0;l<num_ir;l++) {
          fprintf(outfile,"%3d   ",scf_info[l].nclosed);
        }
      fprintf(outfile,"\n  SOCC:            ");
      for(l=0;l<num_ir;l++) {
          fprintf(outfile,"%3d   ",scf_info[l].nopen);
        }
      fprintf(outfile,"\n\n\n");
    }

/* Setup occupation data for the calculation */

  n_open = 0;

/* Convert occupation data to symmetry info */
      
      for(i=0; i < num_ir; i++) {
	  s=&scf_info[i];
	  
	  if (s->nopen || s->nhalf) n_open++;
	  if (nn = s->num_so) {
	      nc = s->nclosed;
	      no = s->nopen;
	      nh = s->nhalf;
	      
	      if (  (nn < nc + no + nh)
		    ||(nc < 0)
		    ||(no < 0)
		    ||(nh < 0)) {
		  const char* fmt
		      = "cscf: invalid number of electrons in irrep %d\n";
		  fprintf(stderr,fmt,i);
		  fprintf(outfile,fmt,i);
		  exit(3);
		} 
	      
	      if(iopen) {
		  s->fock_eff = (double *) init_array(ioff[nn]);
		  s->fock_open = (double *) init_array(ioff[nn]);
		  s->pmato = (double *) init_array(ioff[nn]);
		  s->dpmato = (double *) init_array(ioff[nn]);
		  s->gmato = (double *) init_array(ioff[nn]);
		}
	      if(twocon) {
		  s->pmat2 = (double *) init_array(ioff[nn]);
		  s->pmato2 = (double *) init_array(ioff[nn]);
		}
	      
	      for (j=0; j < nc ; j++) {
		  s->occ_num[j] = 2.0;
		}
	      for (j=nc; j < nc+no ; j++) {
		  s->occ_num[j] = 1.0;
		}
	      if(nh) s->occ_num[nc+no] = nh*0.5;
	    }
	}
      
      optri = n_open*(n_open+1)/2;
      
      if (iopen) {
	  int mm1=1,mm2=1;
	  if(optri == 0){
	      fprintf(outfile,
                 "\nNot an open shell molecule. Re-check opentype!!!!\n\n");
	      exit(1);
	    }
	 
	  if(iter == 0 || print & 1){
	  fprintf(outfile,"\n  open-shell energy coeffs\n");
	  fprintf(outfile,"  open shell pair    alpha         beta\n");}
	  optri = (n_open*(n_open+1))/2;
      
	  alpha = (double *) init_array(ioff[optri]);
	  beta = (double *) init_array(ioff[optri]);
	  
	  if (twocon) {
	      if(n_open == 2) {
		  alpha[0] = 0.0;
		  alpha[1] = 0.0;
		  alpha[2] = 0.0;
		  beta[0] = 0.0;
		  beta[1] = 1.0;
		  beta[2] = 0.0;
		}
	      else {
		  fprintf(outfile,"this program cannot handle same symmetry\n");
		  fprintf(outfile," tcscf. try SCFX\n");
		  exit(1);
		}
	    }
	  else if(singlet) {
	      if(n_open == 2) {
		  alpha[0] = 0.0;
		  alpha[0] = 0.0;
		  alpha[1] = 0.0;
		  alpha[2] = 0.0;
		  beta[0] = 1.0;
		  beta[1] = -3.0;
		  beta[2] = 1.0;
		}
	      else {
		  fprintf(outfile,"this program cannot handle same symmetry\n");
		  fprintf(outfile," singlets. try SCFX\n");
		  exit(1);
		}
	    }
	  else if(hsos) {
	      for(i=0; i < optri ; i++) {
		  alpha[i]=0.0;
		  beta[i]=1.0;
		}
	    }
	  else {
	      for (i=0; i < optri ; i++) {
		  errcod = ip_data("ALPHA","%lf",&alpha[i],1,i);
		  errcod = ip_data("BETA","%lf",&beta[i],1,i);
		  beta[i] = -beta[i];
		}
	    }
	  for (i=0; i < optri; i++) {
	      if(iter == 0 || print & 1)
		  fprintf(outfile,"        %d  %d       %f     %f\n",mm1,mm2,
			  alpha[i],-beta[i]);
	      mm2++;
	      if (mm2 > mm1) {
		  mm1++;
		  mm2 = 1;
		}
	    }
	}
}       



