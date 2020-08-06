#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "tipsydefs.h"
#include "gadgetdefs.h"
#include <math.h>

struct gas_particle *gas_particles, *gp;
struct dark_particle *dark_particles, *dp;
struct star_particle *star_particles, *sp;
struct aux_gas_data *aux_gas_particles, *auxgp;
struct aux_star_data *aux_star_particles, *auxsp;
struct dump header;
int *pid;
char basename[200], binname[200], auxname[200], idnumname[200], outname[200];

double unit_Mass = 1;
double unit_Density = 1;
double unit_Velocity = 1;
double unit_Length = 1;
double unit_Time = 1;
double UnitMass_in_g = 1.989e43;
double UnitLength_in_cm = 3.08568e+24;
double UnitVelocity_in_cm_per_s = 100000;
double atime;

//#define M_PI 3.1415926535897932
#define BOXSIZE 500.0
#define OMEGA0 0.258
#define OMEGAL 0.742
#define HUBBLEPARAM 0.72
#define CM_PER_MPC 3.086e24


#ifdef MARK64
void read_header(struct dump *, FILE *);
#endif
void cosmounits();
void get_snap_string(int snapnum, char *snapstr)
{
  char str[2];

  sprintf(str, "%d", snapnum);
  if (snapnum < 0)
    {printf("snapnum < 0! Force Exit!\n"); exit(-1);}
  else if (snapnum < 10)
    strcat(strcpy(snapstr,"00"), str);
  else if (snapnum < 100)
    strcat(strcpy(snapstr,"0"), str);
  else if (snapnum < 1000)
    strcat(strcpy(snapstr,""), str);
  else
    {printf("snapnum > 999! Force Exit!\n"); exit(-1);}
}

void get_filenames(char *snapbase, int snapnum)
{
  char snapstr[4];
  get_snap_string(snapnum, snapstr);
  sprintf(outname, "%s_%s", snapbase, snapstr);
  sprintf(binname, "%s_%s.bin", snapbase, snapstr);
  sprintf(auxname, "%s_%s.aux", snapbase, snapstr);
  sprintf(idnumname, "%s_%s.idnum", snapbase, snapstr);
}

void ReadData()
{
  FILE *fin;
  if(!(fin = fopen(binname, "r")))
    {printf("Can't open bin file! Force exit.\n"); exit(-1);}
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&header, fin);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&header,sizeof(header),1,fin) == 0)
    {printf("BREAK!\n"); exit(-1);}
#endif
  printf("*** Header information ***\n");
  printf("header.time: %f\n", header.time);
  printf("header.nbodies: %d\n", header.nbodies);
  printf("header.ndim: %d\n", header.ndim);
  printf("header.nsph: %d\n", header.nsph);
  printf("header.ndark: %d\n", header.ndark);
  printf("header.nstar: %d\n", header.nstar);
  printf("**************************\n");
  if(header.nsph != 0) {
    gp = (struct gas_particle *)
      malloc(header.nsph*sizeof(*gp));
    auxgp = (struct aux_gas_data *)
      malloc(header.nsph*sizeof(*auxgp));
  }
  if(header.ndark != 0) {
    dp = (struct dark_particle *)
      malloc(header.ndark*sizeof(*dp));
  }
  if(header.nstar != 0) {
    sp = (struct star_particle *)
      malloc(header.nstar*sizeof(*sp));
    auxsp = (struct aux_star_data *)
      malloc(header.nstar*sizeof(*auxsp));
  }
  pid = (int *) 
    malloc((header.nsph + header.ndark + header.nstar)*sizeof(int));

  fread((char *)gp,sizeof(struct gas_particle),
	header.nsph,fin) ;
  fread((char *)dp,sizeof(struct dark_particle),
	header.ndark,fin) ;
  fread((char *)sp,sizeof(struct star_particle),
	header.nstar,fin) ;
  printf("Reading .bin file done.\n");
  fclose(fin);

  if(!(fin = fopen(auxname, "r")))
    {printf("Can't open aux file! Force exit.\n"); exit(-1);}
  fread((char *)auxgp,sizeof(struct aux_gas_data),
	header.nsph,fin);
  fread((char *)auxsp,sizeof(struct aux_star_data),
	header.nstar,fin);
  fclose(fin);
  printf("Reading .aux file done.\n");

  if(!(fin = fopen(idnumname, "r")))
    {printf("Can't open idnum file! Force exit.\n"); exit(-1);}
  fread(pid, sizeof(int),
	(header.nsph + header.ndark + header.nstar),fin);
  fclose(fin);
  printf("Reading .idnum file done.\n");
}

void WriteData(char *filename)
{
  int i, k;
  int blksize, ntot;
  float fp[3];
  float fp1;

  FILE *fout;
  fout = fopen(filename, "w");
#define SKIP  {fwrite(&blksize,sizeof(int),1,fout);}

  cosmounits();
  /* Header: */
  blksize = 256;
  SKIP;
  for(i=0;i<6;i++) io_header.npart[i] = 0;
  io_header.npart[0] = header.nsph;
  io_header.npart[1] = header.ndark;
  io_header.npart[4] = header.nstar;
  for(i=0;i<6;i++) io_header.npartTotal[i] = io_header.npart[i];
  for(i=0;i<6;i++) io_header.mass[i] = 0.0;
  io_header.mass[1] = dp[0].mass * unit_Mass;
  io_header.time = header.time;
  io_header.redshift = 1. / header.time - 1.;
  io_header.flag_sfr = 1;
  io_header.flag_feedback = 1;
  io_header.flag_cooling = 1;
  io_header.num_files = 1;
  io_header.BoxSize = BOXSIZE;
  io_header.Omega0 = OMEGA0;
  io_header.OmegaLambda = OMEGAL;
  io_header.HubbleParam = HUBBLEPARAM;
  fwrite(&io_header, sizeof(io_header), 1, fout);
  SKIP;
  /* POS */
  ntot = header.nsph + header.nstar + header.ndark;
  blksize = ntot * 3 * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++){
    for(k=0;k<3;k++)
      fp[k] = gp[i].pos[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  for(i=0;i<header.ndark;i++){
    for(k=0;k<3;k++)
      fp[k] = dp[i].pos[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  for(i=0;i<header.nstar;i++){
    for(k=0;k<3;k++)
      fp[k] = sp[i].pos[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  SKIP;

  /* VEL */
  blksize = ntot * 3 * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++){
    for(k=0;k<3;k++)
      fp[k] = gp[i].vel[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  for(i=0;i<header.ndark;i++){
    for(k=0;k<3;k++)
      fp[k] = dp[i].vel[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  for(i=0;i<header.nstar;i++){
    for(k=0;k<3;k++)
      fp[k] = sp[i].vel[k] * unit_Length;
    fwrite(fp, sizeof(float), 3, fout);
  }
  SKIP;

  /* ID */
  blksize = ntot * sizeof(int);
  /* printf("ID = %d %d\n",pid[0], pid[header.nsph]); */
  SKIP;
  fwrite(pid, sizeof(int), ntot, fout);
  SKIP;

  /* MASS */
  blksize = (header.nsph + header.nstar) * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++){
    fp1 = gp[i].mass * unit_Mass;
    fwrite(&fp1, sizeof(float), 1, fout);
  }
  for(i=0;i<header.nstar;i++){
    fp1 = sp[i].mass * unit_Mass;
    fwrite(&fp1, sizeof(float), 1, fout);
  }
  SKIP;

  /* U */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++)
    fwrite(&gp[i].temp, sizeof(float), 1, fout);
  SKIP;

  /* RHO */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++){
    fp1 = gp[i].rho * unit_Density;
    fwrite(&fp1, sizeof(float), 1, fout);
  }
  SKIP;

  /* NE */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++)
    fwrite(&auxgp[i].ne, sizeof(float), 1, fout);
  SKIP;

  /* NH*/
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++)
    fwrite(&auxgp[i].nh, sizeof(float), 1, fout);
  SKIP;

  /* HSML */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++){
    fp1 = gp[i].hsmooth * unit_Length;
    fwrite(&fp1, sizeof(float), 1, fout);
  }
  SKIP;

  /* SFR */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++)
    fwrite(&auxgp[i].sfr, sizeof(float), 1, fout);
  SKIP;

  /* DELAYT */
  blksize = header.nsph * sizeof(float);
  SKIP;
  for(i=0;i<header.nsph;i++)
    fwrite(&auxgp[i].delaytime, sizeof(float), 1, fout);
  SKIP;

  /* RHO */
  blksize = header.nstar * sizeof(float);
  SKIP;
  for(i=0;i<header.nstar;i++)
    fwrite(&auxsp[i].age, sizeof(float), 1, fout);
  SKIP;
  fclose(fout);
}

int main(int argc, char **argv)
{
  char snapbase[200];
  int snapnum;  

  strcpy(snapbase, argv[1]);
  snapnum = atoi(argv[2]);

  get_filenames(snapbase, snapnum);

  printf("I am reading %s, master!\n", outname);

  fprintf(stderr,"read time %f\n",header.time);

  ReadData();

  WriteData(outname);

  if(header.nsph != 0) {
    free(gas_particles);
    free(aux_gas_particles);
  }
  if(header.ndark != 0) {
    free(dark_particles);
  }
  if(header.nstar != 0) {
    free(star_particles);
    free(aux_star_particles);
  }
  free(pid);
  return(0);
}

void read_header(struct dump *head, FILE *ftipsy ) {
  fread((char *)&head->time, sizeof(head->time), 1, ftipsy);
  fread((char *)&head->nbodies, sizeof(head->nbodies), 1, ftipsy);
  fread((char *)&head->ndim, sizeof(head->ndim), 1, ftipsy);
  fread((char *)&head->nsph, sizeof(head->nsph), 1, ftipsy);
  fread((char *)&head->ndark, sizeof(head->ndark), 1, ftipsy);
  fread((char *)&head->nstar, sizeof(head->nstar), 1, ftipsy);
  //fread((char *)&head->pad, sizeof(head->pad), 1, ftipsy);
}


void cosmounits(void)
{
    unit_Time=sqrt(8*M_PI/3)*CM_PER_MPC/(100*HUBBLEPARAM*1.e5);
    unit_Density=1.8791E-29*HUBBLEPARAM*HUBBLEPARAM;
    unit_Length=BOXSIZE*CM_PER_MPC*1.e-3;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length/(HUBBLEPARAM*HUBBLEPARAM);
    unit_Velocity=unit_Length/unit_Time;
    
    unit_Mass /= UnitMass_in_g;
    unit_Velocity /= (UnitVelocity_in_cm_per_s/sqrt(atime));
    unit_Length /= UnitLength_in_cm;
    unit_Density /= (UnitMass_in_g/pow(UnitLength_in_cm, 3));
    return;
}
