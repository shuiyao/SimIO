/* Adopted from the tipsy-tools */


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "tipsydefs.h"

int flag_gasmetal = 0;
int flag_starmetal = 0;
int flag_age = 0;
int flag_tmax = 0;
int flag_sfr = 0;
int flag_delayt = 0;
int flag_ne = 0;
int flag_nh = 0;
int flag_nspawn = 0;
int flag_nrec = 0;
#ifdef OUTPUT_ALPHA
int flag_alpha = 0;
#endif

#ifdef MARK64
void read_header(struct dump *, FILE *);
#endif

int main(int argc, char **argv)
{
    struct aux_gas_data *aux_gas_particles, *gp, *lastgp;
    struct aux_star_data *aux_star_particles, *sp, *lastsp;
    struct dump header;
    int startstamp=0;
    char basename[200], filename[200], binname[200];
    char sfrfile[200], gasmetalfile[200], tmaxfile[200]; //gas
    char delaytfile[200], nefile[200], nhfile[200], nspawnfile[200];
#ifdef OUTPUT_ALPHA
    char alphafile[200];
#endif
    char nrecfile[200];
    char starmetalfile[200], agefile[200];
    FILE *fd, *outputfile;

    int parse_input(int argc, char **argv, char *basename);

    printf("We have %d parameters here, master.\n", argc);
/*     strcpy(basename, argv[1]); */
    parse_input(argc, argv, basename);

    strcat(strcpy(binname, basename), ".bin");
    strcat(strcpy(filename, basename), ".aux");
    strcat(strcpy(gasmetalfile, basename), ".gasmetal");
    strcat(strcpy(tmaxfile, basename), ".tmax");
    strcat(strcpy(agefile, basename), ".age");
    strcat(strcpy(starmetalfile, basename), ".starmetal");
    strcat(strcpy(sfrfile, basename), ".sfr");
    strcat(strcpy(nefile, basename), ".ne");
    strcat(strcpy(nhfile, basename), ".nh");
    strcat(strcpy(delaytfile, basename), ".delayt");
    strcat(strcpy(nspawnfile, basename), ".nspawn");
    strcat(strcpy(nrecfile, basename), ".nrec");
#ifdef OUTPUT_ALPHA
    strcat(strcpy(alphafile, basename), ".alpha");
#endif

    printf("Basename: %s\n", basename);
    while(startstamp == 0)
      {
	startstamp = 1;

	/*     First, read header from the .bin file */
	printf("I am reading %s, master!\n", binname);
    
	if(!(fd = fopen(binname, "r")))
	  {printf("Can't open input file! Force exit.\n"); exit(-1);}


#ifdef MARK64
	printf("I am reading for 64 machine!\n");
	read_header(&header, fd);
#else
	printf("I am reading for 32 machine!\n");
	if(fread((char *)&header,sizeof(header),1,fd) == 0)
	  {printf("BREAK!\n"); break;}
#endif
/* 	headerinfo(header); */
	printf("*** Header information ***\n");
	printf("header.time: %f\n", header.time);
	printf("header.nbodies: %d\n", header.nbodies);
	printf("header.ndim: %d\n", header.ndim);
	printf("header.nsph: %d\n", header.nsph);
	printf("header.ndark: %d\n", header.ndark);
	printf("header.nstar: %d\n", header.nstar);
	printf("**************************\n");
	fclose(fd);

	/* 	Secondly, read .aux file */
	printf("I am reading %s, master!\n", filename);
	if(!(fd = fopen(filename, "r")))
	  {printf("Can't open input file! Force exit.\n"); exit(-1);}

	if(header.nsph != 0) {
	    aux_gas_particles = (struct aux_gas_data *)
				malloc(header.nsph*sizeof(*aux_gas_particles));
	    if(aux_gas_particles == NULL) {
		printf("<sorry, no memory for gas particles, master>\n") ;
		return -1;
	    }
	}
	if(header.nstar != 0) {
	    aux_star_particles = (struct aux_star_data *)
				malloc(header.nstar*sizeof(*aux_star_particles));
	    if(aux_star_particles == NULL) {
		printf("<sorry, no memory for star particles, master>\n") ;
		return -1;
	    }
        }
    
	fread((char *)aux_gas_particles,sizeof(struct aux_gas_data),
			 header.nsph,fd) ;
	fread((char *)aux_star_particles,sizeof(struct aux_star_data),
			 header.nstar,fd) ;

	lastgp = aux_gas_particles + header.nsph ;
	lastsp = aux_star_particles + header.nstar ;
	printf("Reading done.\n");

	printf("Writing age, sfr, metals, tmax ...");
	if(flag_age) {
	  outputfile = fopen(agefile, "w");
	  for(sp=aux_star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g\n",sp->age);
	  }
	  fclose(outputfile);
	}
	if(flag_sfr) {
	  outputfile = fopen(sfrfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->sfr);
	  }
	  fclose(outputfile);
	}
	if(flag_tmax) {
	  outputfile = fopen(tmaxfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->tmax);
	  }
	  fclose(outputfile);
	} 
	if(flag_delayt) {
	  outputfile = fopen(delaytfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->delaytime);
	  }
	  fclose(outputfile);
	} 
	if(flag_ne) {
	  outputfile = fopen(nefile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->ne);
	  }
	  fclose(outputfile);
	} 
	if(flag_nh) {
	  outputfile = fopen(nhfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->nh);
	  }
	  fclose(outputfile);
	} 
	if(flag_gasmetal) {
	  outputfile = fopen(gasmetalfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g %g %g %g\n",gp->metal[0],gp->metal[1],gp->metal[2],gp->metal[3]);
	  }
	  fclose(outputfile);
	}
	if(flag_starmetal) {
	  outputfile = fopen(starmetalfile, "w");
	  for(sp=aux_star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g %g %g %g\n",sp->metal[0],sp->metal[1],sp->metal[2],sp->metal[3]);
	  }
	  fclose(outputfile);
	}
	if(flag_nspawn) {
	  outputfile = fopen(nspawnfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%d\n",gp->nspawn);
	  }
	  fclose(outputfile);
	}
	if(flag_nrec) {
	  outputfile = fopen(nrecfile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%d\n",gp->nrec);
	  }
	  for(sp=aux_star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%d\n",sp->nrec);
	  }
	  fclose(outputfile);
	}
#ifdef OUTPUT_ALPHA
	if(flag_alpha) {
	  outputfile = fopen(alphafile, "w");
	  for(gp=aux_gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->alpha);
	  }
	  fclose(outputfile);
	}
#endif
	printf(" done.\n");

	fprintf(stderr,"read time %f\n",header.time) ;
	if(header.nsph != 0) {
	  free(gas_particles);
	}
	if(header.ndark != 0) {
	  free(dark_particles);
	}
	if(header.nstar != 0) {
	  free(star_particles);
	}
	fclose(fd);
    }
    return(0);
}

int headerinfo(struct dump *header)
{
  printf("*** Header information ***\n");
  printf("header.time: %f\n", header->time);
  printf("header.nbodies: %d\n", header->nbodies);
  printf("header.ndim: %d\n", header->ndim);
  printf("header.nsph: %d\n", header->nsph);
  printf("header.ndark: %d\n", header->ndark);
  printf("header.nstar: %d\n", header->nstar);
  printf("**************************\n");
}

int parse_input(int argc, char **argv, char *basename)
{
  int i=1;
  while(i < argc) {
    printf("%s\n", argv[i]);
    if (!strcmp(argv[i], "-age"))
      {flag_age=1; ++i;}
    else if (!strcmp(argv[i], "-tmax"))
      {flag_tmax=1; ++i;}
    else if (!strcmp(argv[i], "-gasmetal"))
      {flag_gasmetal=1; ++i;}
    else if (!strcmp(argv[i], "-starmetal"))
      {flag_starmetal=1; ++i;}
    else if (!strcmp(argv[i], "-sfr"))
      {flag_sfr=1; ++i;}
    else if (!strcmp(argv[i], "-delayt"))
      {flag_delayt=1; ++i;}
    else if (!strcmp(argv[i], "-ne"))
      {flag_ne=1; ++i;}
    else if (!strcmp(argv[i], "-nh"))
      {flag_nh=1; ++i;}
    else if (!strcmp(argv[i], "-nspawn"))
      {flag_nspawn=1; ++i;}
    else if (!strcmp(argv[i], "-nrec"))
      {flag_nrec=1; ++i;}
#ifdef OUTPUT_ALPHA
    else if (!strcmp(argv[i], "-alpha"))
      {flag_alpha=1; ++i;}
#endif
    else {
      strcpy(basename, argv[i]); ++i;
    }
  }
  return 0;
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
