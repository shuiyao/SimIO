#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64	

#define  BOLTZMANN   1.3806e-16
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  GAMMA         (5.0/3)

#define  H_MASSFRAC    0.76

#define NMETALS 4


struct io_header_1
{
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
};

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass, Rho, Temp, Ne, Nh, Hsml, Sfr, age, metal[NMETALS], DelayTime, Tmax, Phi;
  int Nspawn;
  int    Flag;
};
