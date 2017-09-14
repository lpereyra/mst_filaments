#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

char   message[200];
struct cosmoparam cp;
struct SnapST    snap;
int    size_real;
int    size_int;
int    nfrac;
#ifdef MCRITIC
type_real m_critica;
#endif

void init_variables(int argc, char **argv){
  FILE *pfin;
  char filename[200];

  RED("Initializing variables...\n");

  size_real = sizeof(type_real);
  fprintf(stdout,"type_real: %d\n",size_real);

  size_int  = sizeof(type_int);
  fprintf(stdout,"type_int:  %d\n",size_int);

  sprintf(filename,"%s",argv[1]);
  if(!(pfin=fopen(filename,"r")))
  {
    sprintf(message,"can't open file `%s` \n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&snap.nfiles))
  {
    sprintf(message,"can't read file `%s`\nneed # of snapshots\n",filename);RED(message);
    exit(0);
  }
    
  if(!fscanf(pfin,"%s  \n",snap.root))
  {
    sprintf(message,"can't read file `%s`\nneed snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",snap.name))
  {
    sprintf(message,"can't read file `%s`\nneed snapname\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf \n",&cp.soft))
  {
    sprintf(message,"can't read file `%s`\nneed softening of simulation\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&nfrac))
  {
    sprintf(message,"can't read file `%s`\nneed identification steps\n",filename);RED(message);
    exit(0);
  }

  #ifdef MCRITIC

  #ifdef PRECDOUBLE
    if(!fscanf(pfin,"%lf \n",&m_critica))
  #else
    if(!fscanf(pfin,"%f \n",&m_critica))
  #endif
    {
      sprintf(message,"can't read file `%s`\nneed M_CRITIC\n",filename);RED(message);
      exit(0);
    }
  #endif

  fclose(pfin);
  
  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
  #ifdef MCRITIC
  sprintf(message,"M_CRIT [10^10 Msol / h]  %f\n",m_critica);RED(message);
  #endif

  BLUE("********** Makefile Options ***********\n");
  #ifdef DEBUG
  BLUE("  DEBUG\n");
  #endif
  #ifdef PERIODIC
  BLUE("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  BLUE("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  BLUE("  LONGIDS\n");
  #endif
  #ifdef MPC
  sprintf(message,"  POSFACTOR = 1000.\n");BLUE(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef energies
  BLUE("  Energies\n");
  #endif
  #ifdef LOCK
  BLUE("  USING LOCKS\n");
  #endif

  GREEN("END\n");
}
