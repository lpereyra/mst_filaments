#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "s2plot.h"

#define POSFACTOR 1000
const float lbox = 100.; // in Mpc
int Npart,ngrup,nseg,gmax;
int NNN;
float xmin,ymin,zmin;
float xmax,ymax,zmax;
char  identificacion[200];
float fof[2];

struct io_header{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header;

struct grup_data
{
  int        save;
  int        id;
  float      Pos[3];
  int        NumPart;
  int        *list;
};

struct partstd
{
  int size;
  int *idfil;
  float Pos[3];
};

struct partsimu
{
  float x;
  float y;
  float z;
};

struct segmentstd
{
  int   size;
  int   *list;
  int   flag;
  float len;
  float elong;
  float rms;
  int   size_part;
  int   *list_part;
};

struct grup_data *Gr;
struct partstd *P;
struct partsimu *PSnap;
struct segmentstd *seg;

void read_part();
void read_grup_fof();
void read_segment();
void read_gadget();

int main(int argc, char **argv)
{
  int i,j,k;
  float width;
  char xopt[] = "BCDET";
  char yopt[] = "BCDET";
  char zopt[] = "BCDET";
  char filename[200];
  FILE *pf;
  float *xtmp, *ytmp, *ztmp;
  //int *grupo_mask;
  int migrupo;
  int contador;
  int colors[6] = {2,3,4,5,6,7}; //colors = ('r', 'g', 'b', 'c', 'm', 'y')

  NNN = atoi(argv[1]);
  fof[0] = 1.0                ;
  fof[1] = 200.0              ;

  read_grup_fof();
  read_segment();
  read_part();

  xmin = ymin = zmin = 0.0;
  xmax = ymax = zmax = lbox;

  s2opend("/?",argc,argv);  // Open the display
  ss2spt(1);      					// Generate new projection type
  ss2sbc(0., 0., 0.);				// Set background colour
  s2swin(xmin,xmax,ymin,ymax,zmin,zmax);	// Set the window coordinates
  s2box(xopt,0,0,yopt,0,0,zopt,0,0);   		// Draw coordinate box
  srand48(-15); 	          // semilla random

  read_gadget();

  contador = 0;
  for(i=0;i<nseg;i++)
  {

    xtmp = (float *) malloc(seg[i].size*sizeof(float));
    ytmp = (float *) malloc(seg[i].size*sizeof(float));
    ztmp = (float *) malloc(seg[i].size*sizeof(float));

    for(k=0;k<seg[i].size;k++)
    {
      xtmp[k] = Gr[seg[i].list[k]].Pos[0];
      ytmp[k] = Gr[seg[i].list[k]].Pos[1];
      ztmp[k] = Gr[seg[i].list[k]].Pos[2];
      migrupo = Gr[seg[i].list[k]].save;

      if(k==0) continue;

      if((xtmp[k]-xtmp[k-1])> 0.5*lbox) xtmp[k] -= lbox;
      if((xtmp[k]-xtmp[k-1])<-0.5*lbox) xtmp[k] += lbox;
      if((ytmp[k]-ytmp[k-1])> 0.5*lbox) ytmp[k] -= lbox;
      if((ytmp[k]-ytmp[k-1])<-0.5*lbox) ytmp[k] += lbox;
      if((ztmp[k]-ztmp[k-1])> 0.5*lbox) ztmp[k] -= lbox;
      if((ztmp[k]-ztmp[k-1])<-0.5*lbox) ztmp[k] += lbox;

      //assert(migrupo<(gmax+1));
    }

    
    //if(seg[i].len<5.*POSFACTOR || seg[i].len>15.*POSFACTOR) goto JUMP;
    
    if(seg[i].flag!=2) goto JUMP;

    s2sci(colors[contador%6]);

    s2slw(1.0);

    //if(contador==1)
    for(j=0;j<seg[i].size_part;j++)
      if(drand48()<0.1)
        s2pt1(PSnap[seg[i].list_part[j]].x, \
              PSnap[seg[i].list_part[j]].y, \
              PSnap[seg[i].list_part[j]].z,1);      
    
    //if(contador==1)
      printf("%d\n",seg[i].size);

    s2slw(2.0);

    //if(contador==1)
      s2line(seg[i].size,xtmp,ytmp,ztmp);	// Draw the poly-line

    //if(contador==1)
    //{
      s2slw(1.0);
      s2sci(1);
      for(k=0;k<seg[i].size;k++)
        for(j=0;j<Gr[seg[i].list[k]].NumPart;j++)
       	  if(drand48()<0.1)
            s2pt1(PSnap[Gr[seg[i].list[k]].list[j]].x, \
                  PSnap[Gr[seg[i].list[k]].list[j]].y, \
                  PSnap[Gr[seg[i].list[k]].list[j]].z,1);
    //}

    contador++;

    JUMP:;

    free(xtmp);
    free(ytmp);
    free(ztmp);

  }

  fprintf(stdout,"Quedan %d\n",contador);
  fflush(stdout);

  s2show(1);

  return(EXIT_SUCCESS);

}

void read_part()
{
  char filename[200];
  int  i,j;
  FILE *pfin;

  sprintf(filename,"%.4d_particulas_0.79_0.17.bin",NNN);
  pfin = fopen(filename,"rb"); 
  fread(&j,sizeof(int),1,pfin);

  assert(j==nseg);

  fread(&Npart,sizeof(int),1,pfin);

  fprintf(stdout,"Abre el archivo %s con %d Part\n",filename,Npart);
  fflush(stdout);
  
  for(i=0;i<nseg; i++)
  {
    fread(&seg[i].size_part,sizeof(int),1,pfin);
    seg[i].list_part = (int *) malloc(seg[i].size_part*sizeof(int));

    for(j=0;j<seg[i].size_part;j++)
      fread(&seg[i].list_part[j],sizeof(int),1,pfin);
  }

  fprintf(stdout,"Particulas %d\n",Npart);
  fflush(stdout);

  fclose(pfin);

}

void read_grup_fof()
{
  char  filename[200];
  int   i;
  float ll1;
  FILE  *pfin;
 
  ll1 = 1./cbrt(fof[1]+1);
  
  sprintf(filename,"../fof_grid/%.2f_centros.bin",ll1);
  pfin = fopen(filename,"rb"); 

  fread(&ngrup,sizeof(int),1,pfin);

  fprintf(stdout,"Grupos %d\n",ngrup);
  fflush(stdout);

  Gr = (struct grup_data *) malloc(ngrup*sizeof(struct grup_data));

  for(i=0;i<ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(int),1,pfin);
    fread(&Gr[i].id,sizeof(int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[1],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[2],sizeof(float),1,pfin);
    fread(&Gr[i].NumPart,sizeof(int),1,pfin);

    Gr[i].Pos[0] /= POSFACTOR;
    Gr[i].Pos[1] /= POSFACTOR;
    Gr[i].Pos[2] /= POSFACTOR;
  }

  fclose(pfin);
}

void read_segment()
{
  char  filename[200];
  int   i,k;
  float ll0,ll1;
  FILE  *pf;
 
  ll0 = 1./cbrt(fof[0]+1);
  ll1 = 1./cbrt(fof[1]+1);

  //sprintf(filename,"../fof_grid/mst/segmentos_%.2f_%.2f.bin",ll0,ll1);
  sprintf(filename,"%.4d_segmentos_%.2f_%.2f.bin",NNN,ll0,ll1);

  pf = fopen(filename,"rb"); 

  fread(&nseg,sizeof(int),1,pf);

  fprintf(stdout,"Segmentos %d\n",nseg);
  fflush(stdout);

  seg = (struct segmentstd *) malloc(nseg*sizeof(struct segmentstd));

  for(i=0;i<nseg;i++)
  {
    fread(&seg[i].size,sizeof(int),1,pf);

    seg[i].list = (int *) malloc(seg[i].size*sizeof(int));

    for(k=0;k<seg[i].size;k++)
    {
      fread(&seg[i].list[k],sizeof(int),1,pf);
    }

    seg[i].size_part = 0;    
  }

  fclose(pf);

  sprintf(filename,"%.4d_propiedades_%.2f_%.2f.bin",NNN,ll0,ll1);
  pf = fopen(filename,"rb"); 

  fread(&k,sizeof(int),1,pf);

  assert(k==nseg);

  fprintf(stdout,"Propiedades Segmentos %d\n",nseg);
  fflush(stdout);

  for(i=0;i<nseg;i++)
  {  
    fread(&seg[i].flag,sizeof(int),1,pf);
    fread(&k,sizeof(int),1,pf);
    fread(&seg[i].len,sizeof(float),1,pf);
    fread(&seg[i].elong,sizeof(float),1,pf);
    fread(&seg[i].rms,sizeof(float),1,pf);
    
    assert(k==seg[i].size);
  }

  fclose(pf);

}

void read_gadget(void)
{
  char filename[200];
  int  ifile,ind;
  FILE *pf;
  int d1,d2;
  int  i,j,k,n;
  float r[3];
  XYZ xyz;                                     /* Position */
  COLOUR col;                                  /* Colour */
  float isize;                                 /* Dot size */
  char itrans;                                 /* Transparency type */
  float ialpha;                                /* Alpha channel */
  float ll1;
  float *xtmp, *ytmp, *ztmp;
  int   *index;

  //s2slw(1.0);
  sprintf(filename,"../../snapshot_200");
  pf = fopen(filename,"r");

  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  PSnap = (struct partsimu *) malloc(header.npart[1]*sizeof(struct partsimu));
  index = (int *) malloc((header.npart[1]+1)*sizeof(int)); // le sumo 1

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, j = 0; k < 6; k++)
  {
    for(n = 0; n < header.npart[k]; n++)
    {
      fread(&r[0], sizeof(float), 3, pf);
      if(k == 1)
      { //ONLY KEEP DARK MATTER PARTICLES
        /* 
       	if(drand48()<0.1)
        {
      		xyz.x = r[0];
      		xyz.y = r[1];
		      xyz.z = r[2];
		      col.r = 1.0; 
		      col.g = 1.0;
		      col.b = 1.0;
		      itrans = 't';                             
		      ialpha = 0.5;                             
		      isize  = 1.0;
		      ns2vpa(xyz, col, isize, itrans, ialpha); 
        }
        */     
        PSnap[j].x = r[0];
        PSnap[j].y = r[1];
        PSnap[j].z = r[2];
        j++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  fseek(pf,d1,SEEK_CUR);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, j = 0; k < 6; k++)
  {
    for(n = 0; n < header.npart[k]; n++)
    {
      fread(&i, sizeof(int), 1, pf);
      if(k == 1) /*ONLY KEEP DARK MATTER PARTICLES*/
      {
        index[i] = j;
        j++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fclose(pf);

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);

  ll1 = 1./cbrt(fof[1]+1);
  sprintf(filename,"../fof_grid/%.2f_fof.bin",ll1);
  pf = fopen(filename,"r"); 
  fread(&i,sizeof(int),1,pf);

  fprintf(stdout,"Archivo de Grupos %s\n",filename);
  fflush(stdout); 
  
  assert(i==ngrup);

  for(i=0;i<ngrup;i++)
  {   
    fread(&d1,sizeof(int),1,pf);
    fread(&d2,sizeof(int),1,pf);
    fread(&n,sizeof(int),1,pf);

    assert(n==Gr[i].NumPart);
    Gr[i].list = (int *) malloc(Gr[i].NumPart*sizeof(int));

    for(j=0;j<n;j++)
    {
      fread(&ind,sizeof(int),1,pf);
      ind = index[ind];
      Gr[i].list[j] = ind;
    }

  }
  fclose(pf);

  free(index);

}
