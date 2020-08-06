/* bin2ascii  filename  [nskip]
     * filename = name of tipsy binary file
     * nskip > 0: output every nskip particles, SPH only
       nskip < 0: output every nskip particles, ALL particles
     > file in ascii format

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tipsydefs_n2.h"

void read_header(struct dump *, FILE *);

int
main(argc, argv)
     int argc;
     char **argv;
{
  int i;
  int nskip=100, int_hold;
  struct gas_particle gp[1];
  struct dark_particle dp[1];
  struct star_particle sp[1];
  struct dump header;
  FILE *fp ;
  
  if(argc != 2 && argc != 3){
    fprintf(stderr, "Usage: bin2ascii filename [nskip]\n");
    exit(-1);
  }
  if( argc == 3 ) nskip = atoi(argv[2]);
  
  if( (fp=fopen(argv[1],"r")) == NULL ) {
    fprintf(stderr,"Cannot open file %s\n",argv[1]) ;
    exit(-1);
  }
  
  //fread((char *)&header,sizeof(header),1,fp) ;
  read_header(&header, fp);
  fprintf(stderr,"reading tipsy data %d+%d+%d=%d\n",header.nsph,header.ndark,header.nstar,header.nbodies) ;
  
  if( nskip > 0 ) {
    for(i=0; i < header.nsph ; i++) {
      fread(gp,sizeof(struct gas_particle),1,fp) ;
      int_hold = (int)(gp->metals);
      if(int_hold>0) gp->metals = gp->metals - int_hold;
      if( i%nskip==0) fprintf(stdout,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",i,gp->mass,gp->pos[0],gp->pos[1],gp->pos[2],gp->vel[0],gp->vel[1],gp->vel[2],gp->hsmooth,gp->phi,gp->temp,gp->hsmooth,gp->rho,gp->metals);
    }
  }
  else if( nskip < 0 ) {
    for(i=0; i < header.nsph ; i++) {
      fread(gp,sizeof(struct gas_particle),1,fp) ;
      //if( i%(-nskip)==0 ) fprintf(stdout,"%d %g %g %g %g %g %g %g %g %g\n",i,dp->mass,dp->pos[0],dp->pos[1],dp->pos[2],dp->vel[0],dp->vel[1],dp->vel[2],dp->eps,dp->phi);
    }
    for(i=0; i < header.ndark ; i++) {
      fread(dp,sizeof(struct dark_particle),1,fp) ;
      //if( i%(-nskip)==0 ) fprintf(stdout,"%d %g %g %g %g %g %g %g %g %g\n",i+header.nsph,dp->mass,dp->pos[0],dp->pos[1],dp->pos[2],dp->vel[0],dp->vel[1],dp->vel[2],dp->eps,dp->phi);
    }
    for(i=0; i < header.nstar ; i++) {
      fread(sp,sizeof(struct star_particle),1,fp) ;
      int_hold = (int)(sp->metals);
      if(int_hold>0) sp->metals = sp->metals - int_hold;
      if( i%(-nskip)==0 ) fprintf(stdout,"%d %g %g %g %g %g %g %g %g %g %g\n",i+header.nsph+header.ndark,sp->mass,sp->pos[0],sp->pos[1],sp->pos[2],sp->vel[0],sp->vel[1],sp->vel[2],sp->eps,sp->phi,sp->metals);
    }
  }
  else fprintf(stderr,"nskip cannot be 0\n");
  fclose(fp);
  
  return 0;
}

void read_header(struct dump *head, FILE *ftipsy ) {
  fread((char *)&head->time, sizeof(head->time), 1, ftipsy);
  fread((char *)&head->nbodies, sizeof(head->nbodies), 1, ftipsy);
  fread((char *)&head->ndim, sizeof(head->ndim), 1, ftipsy);
  fread((char *)&head->nsph, sizeof(head->nsph), 1, ftipsy);
  fread((char *)&head->ndark, sizeof(head->ndark), 1, ftipsy);
  fread((char *)&head->nstar, sizeof(head->nstar), 1, ftipsy);
  fread((char *)&head->pad, sizeof(head->pad), 1, ftipsy);
}
