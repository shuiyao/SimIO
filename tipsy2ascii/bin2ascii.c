/* Adopted from the tipsy-tools */


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "tipsydefs.h"

int flag_eps = 0;
int flag_temp = 0;
int flag_rho = 0;
int flag_tform = 0;
int flag_metal = 0;
int flag_basic = 0;
int flag_phi = 0;
int flag_pos = 0;
int flag_vel = 0;
int flag_mass = 0;
int flag_hsmooth = 0;

#ifdef MARK64
void read_header(struct dump *, FILE *);
#endif

int main(int argc, char **argv)
{
    struct gas_particle *gas_particles, *gp, *lastgp;
    struct dark_particle *dark_particles, *dp, *lastdp;
    struct star_particle *star_particles, *sp, *lastsp;
    struct dump header;
    int startstamp=0;
    char basename[200], filename[200], basicoutput[200];
    char metalfile[200], tempfile[200], phifile[200];
    char rhofile[200], tformfile[200], epsfile[200], hsmoothfile[200];
    char velfile[200], massfile[200], posfile[200];
    FILE *fd, *outputfile, *asciifile;

    int parse_input(int argc, char **argv, char *basename);

    printf("We have %d parameters here, master.\n", argc);
/*     strcpy(basename, argv[1]); */
    parse_input(argc, argv, basename);

    strcat(strcpy(filename, basename), ".bin");
    strcat(strcpy(basicoutput, basename), ".ascii");
    strcat(strcpy(posfile, basename), ".pos");
    strcat(strcpy(velfile, basename), ".vel");
    strcat(strcpy(massfile, basename), ".mass");
    strcat(strcpy(metalfile, basename), ".metal");
    strcat(strcpy(tempfile, basename), ".temp");
    strcat(strcpy(rhofile, basename), ".rho");
    strcat(strcpy(tformfile, basename), ".tform");
    strcat(strcpy(epsfile, basename), ".eps");
    strcat(strcpy(hsmoothfile, basename), ".hsmooth");
    strcat(strcpy(phifile, basename), ".phi");

    printf("Basename: %s\n", basename);
    printf("I am reading %s, master!\n", filename);

    if(!(fd = fopen(filename, "r")))
      {printf("Can't open input file! Force exit.\n"); exit(-1);}

    while(startstamp == 0)
      {
	startstamp = 1;
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

	if(header.nsph != 0) {
	    gas_particles = (struct gas_particle *)
				malloc(header.nsph*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		printf("<sorry, no memory for gas particles, master>\n") ;
		return -1;
	    }
	}
	if(header.ndark != 0) {
	    dark_particles = (struct dark_particle *)
				malloc(header.ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, master>\n") ;
		return -1;
	    }
	}
	if(header.nstar != 0) {
	    star_particles = (struct star_particle *)
				malloc(header.nstar*sizeof(*star_particles));
	    if(star_particles == NULL) {
		printf("<sorry, no memory for star particles, master>\n") ;
		return -1;
	    }
        }
    
	fread((char *)gas_particles,sizeof(struct gas_particle),
			 header.nsph,fd) ;
	fread((char *)dark_particles,sizeof(struct dark_particle),
			 header.ndark,fd) ;
	fread((char *)star_particles,sizeof(struct star_particle),
			 header.nstar,fd) ;

	lastgp = gas_particles + header.nsph ;
	lastdp = dark_particles + header.ndark ;
	lastsp = star_particles + header.nstar ;
	printf("Reading done.\n");

#ifdef CHECK_SNAPSHOT	
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    if(gp->pos[0] == 0) printf("BAD GAS POS!\n");
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    if(dp->pos[0] == 0) printf("BAD DARK POS!\n");
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    if(sp->pos[0] == 0) printf("BAD STAR POS!\n");
	  }
#endif

	if(flag_basic) {
	  asciifile = fopen(basicoutput, "w");
	  fprintf(asciifile, "%d %d %d\n" ,header.nbodies, header.nsph,
		  header.nstar) ;
	  fprintf(asciifile,"%d\n",header.ndim) ;
	  fprintf(asciifile,"%g\n",header.time) ;
	  printf("Writing mass ...");
	  for(gp=gas_particles; gp < lastgp ; gp++){
	    fprintf(asciifile,"%g\n",gp->mass);
	  }
	  for(dp=dark_particles; dp < lastdp;  dp++){
	    fprintf(asciifile,"%g\n",dp->mass);
	  }
	  for(sp=star_particles; sp < lastsp; sp++){
	    fprintf(asciifile,"%g\n",sp->mass);
	  }
	  printf(" done.\n");
	  printf("Writing position ...");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(asciifile,"%g\n",gp->pos[0]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(asciifile,"%g\n",dp->pos[0]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(asciifile,"%g\n",sp->pos[0]);
	  }
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(asciifile,"%g\n",gp->pos[1]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(asciifile,"%g\n",dp->pos[1]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(asciifile,"%g\n",sp->pos[1]);
	  }
	  if (header.ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
	      fprintf(asciifile,"%g\n",gp->pos[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
	      fprintf(asciifile,"%g\n",dp->pos[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
	      fprintf(asciifile,"%g\n",sp->pos[2]);
	    }
	  }
	  printf(" done.\n");



	  printf("Writing velocity ...");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(asciifile,"%g\n",gp->vel[0]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(asciifile,"%g\n",dp->vel[0]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(asciifile,"%g\n",sp->vel[0]);
	  }
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(asciifile,"%g\n",gp->vel[1]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(asciifile,"%g\n",dp->vel[1]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(asciifile,"%g\n",sp->vel[1]);
	  }
	  if (header.ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
	      fprintf(asciifile,"%g\n",gp->vel[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
	      fprintf(asciifile,"%g\n",dp->vel[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
	      fprintf(asciifile,"%g\n",sp->vel[2]);
	    }
	  }
	  printf(" done.\n");
	  fclose(asciifile);
	} /* basic file done */

	if(flag_mass) {
	  outputfile = fopen(massfile, "w");
	  printf("Writing mass ...");
	  for(gp=gas_particles; gp < lastgp ; gp++){
	    fprintf(outputfile,"%g\n",gp->mass);
	  }
	  for(dp=dark_particles; dp < lastdp;  dp++){
	    fprintf(outputfile,"%g\n",dp->mass);
	  }
	  for(sp=star_particles; sp < lastsp; sp++){
	    fprintf(outputfile,"%g\n",sp->mass);
	  }
	  printf(" done.\n");
	  fclose(outputfile);
	}

	if(flag_pos) {
	  outputfile = fopen(posfile, "w");
	  printf("Writing posocity ...");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g %g %g\n",gp->pos[0], gp->pos[1], gp->pos[2]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(outputfile,"%g %g %g\n",dp->pos[0], dp->pos[1], dp->pos[2]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g %g %g\n",sp->pos[0], sp->pos[1], sp->pos[2]);	  
	  }
	  printf(" done.\n");
	  fclose(outputfile);
	}

	if(flag_vel) {
	  outputfile = fopen(velfile, "w");
	  printf("Writing velocity ...");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g %g %g\n",gp->vel[0], gp->vel[1], gp->vel[2]);
	  }
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(outputfile,"%g %g %g\n",dp->vel[0], dp->vel[1], dp->vel[2]);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g %g %g\n",sp->vel[0], sp->vel[1], sp->vel[2]);
	  }
	  printf(" done.\n");
	  fclose(outputfile);
	}

	printf("Writing eps, rho, temp, hsmooth, metals, tform ...");
	if(flag_eps) {
	  outputfile = fopen(epsfile, "w");
	  for(dp=dark_particles; dp < lastdp ; dp++) {
	    fprintf(outputfile,"%g\n",dp->eps);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g\n",sp->eps);
	  }
	  fclose(outputfile);
	}
	if(flag_rho) {
	  outputfile = fopen(rhofile, "w");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->rho);
	  }
	  fclose(outputfile);
	}
	if(flag_temp) {
	  outputfile = fopen(tempfile, "w");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->temp);
	  }
	  fclose(outputfile);
	} 
	if(flag_hsmooth) {
	  outputfile = fopen(hsmoothfile, "w");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->hsmooth);
	  }
	  fclose(outputfile);
	}
	if(flag_metal) {
	  outputfile = fopen(metalfile, "w");
	  for(gp=gas_particles; gp < lastgp ; gp++) {
	    fprintf(outputfile,"%g\n",gp->metals);
	  }
	  for(sp=star_particles; sp < lastsp ; sp++) {
	    fprintf(outputfile,"%g\n",sp->metals);
	  }
	  fclose(outputfile);
	}
	if(flag_tform) {
	  outputfile = fopen(tformfile, "w");
	  for(sp=star_particles ; sp < lastsp ; sp++) {
	    if(sp->tform != 0)
	      fprintf(outputfile,"%g\n", 1.);
	    else
	      fprintf(outputfile,"%g\n", -(sp->phi));
	  }
	  fclose(outputfile);
        }
	printf(" done.\n");
	if(flag_phi) {
	  printf("Writing phi ...");
	  outputfile = fopen(phifile, "w");
	  for(gp=gas_particles; gp < lastgp ; gp++){
	    fprintf(outputfile,"%g\n",gp->phi);
	  }
	  /* for(dp=dark_particles; dp < lastdp;  dp++){ */
	  /*   fprintf(outputfile,"%g\n",dp->phi); */
	  /* } */
	  /* for(sp=star_particles; sp < lastsp; sp++){ */
	  /*   fprintf(outputfile,"%g\n",sp->phi); */
	  /* } */
	  fclose(outputfile);
	}
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
    if (!strcmp(argv[i], "-eps"))
      {flag_eps=1; ++i;}
    else if (!strcmp(argv[i], "-temp"))
      {flag_temp=1; ++i;}
    else if (!strcmp(argv[i], "-rho"))
      {flag_rho=1; ++i;}
    else if (!strcmp(argv[i], "-metal"))
      {flag_metal=1; ++i;}
    else if (!strcmp(argv[i], "-tform"))
      {flag_tform=1; ++i;}
    else if (!strcmp(argv[i], "-hsmooth"))
      {flag_hsmooth=1; ++i;}
    else if (!strcmp(argv[i], "-basic"))
      {flag_basic=1; ++i;}
    else if (!strcmp(argv[i], "-phi"))
      {flag_phi=1; ++i;}
    else if (!strcmp(argv[i], "-vel"))
      {flag_vel=1; ++i;}
    else if (!strcmp(argv[i], "-pos"))
      {flag_pos=1; ++i;}
    else if (!strcmp(argv[i], "-mass"))
      {flag_mass=1; ++i;}
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
