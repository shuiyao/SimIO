#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7
#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define  H_MASSFRAC    0.76


double  UnitLength_in_cm=3.085678e21,
        UnitMass_in_g=1.989e43,
               UnitVelocity_in_cm_per_s=1.e5,
               UnitTime_in_s,
               UnitTime_in_Megayears,
               UnitDensity_in_cgs,
               UnitPressure_in_cgs,
               UnitCoolingRate_in_cgs,
               UnitEnergy_in_cgs,
               G,
               Hubble;


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
  int      num_files;	// number of files the snapshot is split into
  double   BoxSize;	// in kpc/h
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 	// H0/100
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


double  Time, Redshift;
double unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity;

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char output_fname[200], input_fname[200], basename[200], hsmlfile[200];
  int  type, files;
  
  if(argc !=2) 
    {
      fprintf(stderr,"usage: snaphead snapshotfile\n");
      exit(-1);
    }

  strcpy(basename, argv[1]);
  hsmlfile[0]= 0;

  sprintf(input_fname, "%s", basename);

  cosmo_setup(input_fname);

  exit(0);
}

int cosmo_setup(char *fname)
{
  FILE *fd;
  char   buf[200];
  int dummy,i;

	sprintf(buf,"%s.0",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  sprintf(buf,"%s",fname);
          if(!(fd=fopen(buf,"r"))) {
	    fprintf(stderr,"can't open file `%s`\n",buf);
	    exit(0);
	  }
	}

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);

      UnitLength_in_cm /= header1.HubbleParam;
      UnitMass_in_g /= header1.HubbleParam;

  Time= header1.time;
  Redshift= header1.time;

      cosmounits();

      for(i=0;i<6;i++) fprintf(stdout,"np[%d]= %d ",i,header1.npartTotal[i]);
      fprintf(stdout,"\nt= %g  z= %g  L= %g  Omega= %g  L= %g  h= %g\n",header1.time,header1.redshift,header1.BoxSize,header1.Omega0,header1.OmegaLambda,header1.HubbleParam);
}



cosmounits()
{
        float L=-1.,h=-1.,totMass=-1.;
	char aline[80],junk[80];
    double Pi=3.14159265358979323846;
    double km=1.E5;
    double Mpc=3.086E24;
    double m_p=1.6726231E-24;     /* proton mass */
    double k_B=1.380622E-16;      /* Boltzman constant */

	totMass = header1.Omega0;
	h = header1.HubbleParam;
	L = header1.BoxSize/h;
	L *= 1.e-3; 	/* convert to comoving Mpc (not h^-1) */

    unit_Time=sqrt(8*Pi/3)*Mpc/(100*h*km);
    unit_Density=1.8791E-29*h*h;
    unit_Length=L*Mpc;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length;
    unit_Velocity=unit_Length/unit_Time;
/*        fprintf(stderr,"COSMO PARAMS:  L=%g Mpc, h=%g, Omega=%g\n",L,h,totMass);
        fprintf(stderr,"UNITS: T=%g rho=%g L=%g M=%g v=%g\n",unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity);*/

        return 0;
}

