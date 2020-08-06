#include "gadgetdefs.h"
#include "tipsydefs_n2.h"
#include "extern.h"

int cosmo_setup(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int dummy;

      if(files>1)
	sprintf(buf,"%s.0",fname);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  fprintf(stderr,"can't open file `%s`\n",buf);
	  exit(0);
	}

      fprintf(stderr,"reading header of `%s'\n",buf);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);

      return(0);

}

int parse_input(int argc, char **argv, char *basename) 
{
	int i;
	char *p;
	int usage(), isalpha();

	i = 1;
	while (i < argc) {
		if (!strcmp(argv[i],"-u")) {
			++i;
			unit_flag = 1;
		}
		else if (!strcmp(argv[i],"-id")) {
			++i;
			id_flag = 1;
		}
		else if (!strcmp(argv[i],"-Nfiles")) {
			++i;
			files = atoi(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-eps")) {
			++i;
			epsilon = atof(argv[i]);
			++i;	  
		}
		else if (!strcmp(argv[i],"-nskip")) {
			++i;
			nskip = atoi(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-Nout")) {
			++i;
			Nout = atoi(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-maxdens")) {
			++i;
			maxdens = atof(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-mindens")) {
			++i;
			mindens = atof(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-maxT")) {
			++i;
			maxT = atof(argv[i]);
			++i;
		}
		else if (!strcmp(argv[i],"-minT")) {
			++i;
			minT = atof(argv[i]);
			++i;
		}
		else if (*argv[i] == '-') {
			p = argv[i];
			++p;
			if (*p == 'd' || *p == 'g' || *p == 's') {
				bDark = 0;
				bGas = 0;
				bStar = 0;
			}
			else {printf("shuiyao:-dgs?"); usage();}
			while (isalpha(*p)) {
				switch (*p) {
					case 'd': bDark = 1; break;
					case 'g': bGas = 1; break;
					case 's': bStar = 1; break;
					default: usage();
				}
				++p;
			}
			++i;
		}
		else {
			strcpy(basename, argv[i]);
			++i;
		}
	}
	if( epsilon < 0. ) {printf("shuiyao:epsilon < 0.\n"); usage();}
	if( files < Nout ) {printf("shuiyao:files < Nout\n"); usage();}
	return 0;
}



/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
/*  fprintf(stderr,"allocating memory...\n");*/

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */
  return(0);

}


int free_memory(void)
{
  Id++;
  P++;
  free(Id);
  free(P);
  return(0);
}


/* this converts from Gadget's units to tipsy "standard" units.  */

int unit_conversion(int type)
{
  int i;
  double MeanWeight;

  for(i=1; i<=NumPart; i++) {
/* Convert units to tipsy standard: boxsize=1, totalmass=Omega */
    if( unit_flag ) {
/*       if (type==4) */
/* 	printf("Before: Mass=%g, %g, %g\n", P[i].Mass, UnitMass_in_g, unit_Mass); */
      P[i].Mass *= UnitMass_in_g/unit_Mass;
/*       if (type==4) */
/* 	printf("After: Mass=%g\n", P[i].Mass); */
      P[i].Pos[0] *= UnitLength_in_cm/unit_Length;
      P[i].Pos[1] *= UnitLength_in_cm/unit_Length;
      P[i].Pos[2] *= UnitLength_in_cm/unit_Length;
      P[i].Pos[0] -= 0.5;
      P[i].Pos[1] -= 0.5;
      P[i].Pos[2] -= 0.5;
      if(Id[i]==666325)
	printf("\n%g %g\n", P[i].Pos[0], P[i].Vel[0]);
#ifdef VELCORR
      double a3 = header.time * header.time * header.time;
      P[i].Vel[0] *= sqrt(a3);
      P[i].Vel[1] *= sqrt(a3);
      P[i].Vel[2] *= sqrt(a3);
#endif      
      P[i].Vel[0] *= UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(header.time);
      P[i].Vel[1] *= UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(header.time);
      P[i].Vel[2] *= UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(header.time);
      // shuiyao: NOTE HEADER.TIME > 0, (INIT_TIME != 0)
      if(Id[i]==666325)
	printf("%g %g %g %g\n", P[i].Vel[0], UnitVelocity_in_cm_per_s, unit_Velocity, sqrt(header.time));
    }
    if(type==0) {
      if( unit_flag ) {
	P[i].Rho *= UnitMass_in_g/(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm*unit_Density);
	P[i].Hsml *= UnitLength_in_cm/unit_Length;
	}
      MeanWeight= 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* P[i].Ne) * PROTONMASS;
      P[i].Temp *= MeanWeight/BOLTZMANN * (GAMMA-1.)
	* UnitEnergy_in_cgs/ UnitMass_in_g; /* temp now in K */
    }
  }
  return(0);
}



cosmounits()
{
        float L=-1.,h=-1.,totMass=-1.;
    double Pi=3.14159265358979323846;
    double km=1.E5;

	totMass = header1.Omega0;
	h = header1.HubbleParam;
	L = header1.BoxSize/h;
	L *= 1.e-3; 	/* convert to comoving Mpc (not h^-1) */

/* set some useful Gadget unit conversions */
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs=UnitMass_in_g/pow(UnitLength_in_cm,3);
    UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
    Hubble = HUBBLE * UnitTime_in_s;

/* conversion to tipsy "standard" units */
    unit_Time=sqrt(8*Pi/3)*CM_PER_MPC/(100*h*km);
    unit_Density=1.8791E-29*h*h;
    unit_Length=L*CM_PER_MPC;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length;
    unit_Velocity=unit_Length/unit_Time;

        fprintf(stderr,"COSMO PARAMS:  L=%g Mpc, h=%g, Omega=%g\n",L,h,totMass);
        fprintf(stderr,"UNITS: T=%g rho=%g L=%g M=%g v=%g\n",unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity);

        return 0;
}

int usage()
{
    fprintf(stderr,"USAGE:\n");
    fprintf(stderr,"snap2tipsy Snap_file\n");
    fprintf(stderr,"    -eps <float>  Softening length in comoving kpc/h.\n");
    fprintf(stderr,"    [-u/-units] Convert to tipsy standard units.\n");
    fprintf(stderr,"    [-dgs] Include gas, dark matter, or stars in computation.\n");
    fprintf(stderr,"    [-id] Output id number file (<Snap_file>.idnum).\n");
    fprintf(stderr,"    [-Nfiles <int>] Number of input files. Default: 1\n");
    fprintf(stderr,"    [-nskip <int>] Output every nskip-th particle. Default: 1\n");
    fprintf(stderr,"    [-Nout <int>] Number of output files (<Nfiles). Default: 1\n");
    fprintf(stderr,"    [-maxdens <float>] Only output gas particles BELOW this density (tipsy units).  Default: 1.e30\n");
    fprintf(stderr,"    [-mindens <float>] Only output gas particles ABOVE this density (tipsy units).  Default: 0\n");
    fprintf(stderr,"    [-maxT <float>] Only output gas particles BELOW this temperature.  Default: 1.e30\n");
    fprintf(stderr,"    [-minT <float>] Only output gas particles ABOVE this temperature.  Default: 0\n");
    fprintf(stderr,"Outputs tipsy binary to <Snap_file>.bin\n");
    exit(-1);
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

void write_header(struct dump *head, FILE *ftipsy ) {
  //head->pad = 0;
  fwrite((char *)&head->time, sizeof(head->time), 1, ftipsy);
  fwrite((char *)&head->nbodies, sizeof(head->nbodies), 1, ftipsy);
  fwrite((char *)&head->ndim, sizeof(head->ndim), 1, ftipsy);
  fwrite((char *)&head->nsph, sizeof(head->nsph), 1, ftipsy);
  fwrite((char *)&head->ndark, sizeof(head->ndark), 1, ftipsy);
  fwrite((char *)&head->nstar, sizeof(head->nstar), 1, ftipsy);
  //fwrite((char *)&head->pad, sizeof(head->pad), 1, ftipsy);
}

