#include "defs.h"

double  UnitLength_in_cm=3.085678e21, UnitMass_in_g=1.989e43, UnitVelocity_in_cm_per_s=1.e5;

double  UnitTime_in_s, UnitDensity_in_cgs, UnitEnergy_in_cgs, G, Hubble;

float epsilon=-1.;
int nskip=1, files=1, Nout=1;
int unit_flag=1,id_flag=1;
int bGas=1, bDark=1, bStar=1;
float mindens=0.,minT=0.;
float maxdens=1.e30,maxT=1.e30;


struct io_header_1 header1;

int     NumPart;

struct particle_data *P;

int   *Id;

double unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity;

struct dump header ;

int ngasout;

struct gas_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real rho;
  Real temp;
  Real hsmooth;
  Real metals[4] ;
  Real sfr ;
  Real phi ;
  Real delaytime ; 
  Real tmax; 
} gasp;

struct dark_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real eps;
  Real phi ;
} darkp;

struct star_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real metals[4] ;
  Real tform ;
  Real eps;
  Real phi ;
} starp;

int load_snapshot(char *, int);

int main(int argc, char **argv)
{
	char output_fname[200], input_fname[200], basename[200], buf[200], id_fname[200], aux_fname[200];
	int  type, ifile, i, junk;
	FILE *fd,*outfile,*idfile,*auxfile;

 	int parse_input(),cosmounits(),output_tipsy_gas(),load_snapshot(),free_memory(),output_tipsy_dark(),output_tipsy_star(),usage(),unit_conversion();

	parse_input(argc,argv,basename); 

  	sprintf(input_fname, "%s", basename);
  	header.ndim =  3;
	for( i=0; i<Nout; i++ ) {
  		if( Nout==1 ) {
			sprintf(output_fname, "%s.bin", basename);
  			if( id_flag ) sprintf(id_fname, "%s.idnum", basename);
			sprintf(aux_fname, "%s.aux", basename);
  		}
		else {
			sprintf(output_fname, "%s.%d.bin", basename,i);
  			if( id_flag ) sprintf(id_fname, "%s.%d.idnum", basename,i);
			sprintf(aux_fname, "%s.%d.aux", basename,i);
		}
  		if(!(outfile = fopen(output_fname,"w"))) {
      			printf("can't open file `%s'\n", output_fname);
      			exit(-1);
  		}
  		if(!(idfile = fopen(id_fname,"w"))) {
      			printf("can't open file `%s'\n", id_fname);
      			exit(-1);
  		}
  		if(!(auxfile = fopen(aux_fname,"w"))) {
      			printf("can't open file `%s'\n", aux_fname);
      			exit(-1);
  		}

/* Load tipsy header */
  		header.nsph = header.ndark = header.nstar = 0;
  		for( ifile=i*files/Nout; ifile<(i+1)*files/Nout; ifile++ ) {
			if( ifile >= files ) continue;
  			if(files>1) sprintf(buf,"%s.%d",input_fname,ifile);
  			else sprintf(buf,"%s",input_fname);
  			if(!(fd=fopen(buf,"r"))) {
    				fprintf(stderr,"can't open file `%s`\n",buf);
    				exit(-1);
  			}
  			fread(&junk, sizeof(int), 1, fd);
  			fread(&header1, sizeof(header1), 1, fd);
  			fclose(fd);
  			header.time = header1.time;
  			if( bGas ) header.nsph += header1.npart[0]/nskip;
  			if( bDark ) header.ndark += header1.npart[1]/nskip;
  			if( bStar ) header.nstar += header1.npart[4]/nskip;
		}
/* output tipsy header */
  		header.nbodies = header.nsph+header.ndark+header.nstar;
  		fwrite(&header, sizeof(header), 1, outfile);
 
/* Set some cosmology units stuff to physical units */ 
  		UnitLength_in_cm /= header1.HubbleParam;
  		UnitMass_in_g /= header1.HubbleParam;
  		cosmounits();
  		epsilon *= 1.4*UnitLength_in_cm/unit_Length;
  		fprintf(stderr,"outfile %d of %d: time= %g  nbodies= %d  nsph= %d  ndark= %d  nstar= %d  eps= %g\n",i+1,Nout,header.time,header.nbodies,header.nsph,header.ndark,header.nstar,epsilon);

  		if( id_flag ) fprintf(stderr,"Outputting IDs to %s\n",output_fname);

  		fprintf(stderr,"Outputting auxiliary info to %s\n",aux_fname);

  		if( bGas ) {
  			fprintf(stderr,"Outputting gas particles in file number ");
  			ngasout = 0;
  			for( ifile=i*files/Nout; ifile<(i+1)*files/Nout; ifile++ ) {
				if( ifile >= files ) continue;
    				if(files>1) sprintf(buf,"%s.%d",input_fname,ifile);
    				else sprintf(buf,"%s",input_fname);
    				load_snapshot(buf, type=0);
				unit_conversion(type=0);
    				output_tipsy_gas(outfile,idfile,auxfile);
    				free_memory(); 
    				fprintf(stderr,"%d ",ifile);
  			}
  			fprintf(stderr,"\n");
  		}
  		else fprintf(stderr,"No gas particles output\n");

  		if( bDark ) {
  			fprintf(stderr,"Outputting dark particles in file number ");
  			for( ifile=i*files/Nout; ifile<(i+1)*files/Nout; ifile++ ) {
				if( ifile >= files ) continue;
    				if(files>1) sprintf(buf,"%s.%d",input_fname,ifile);
    				else sprintf(buf,"%s",input_fname);
    				load_snapshot(buf, type=1);
				unit_conversion(type=1);
    				output_tipsy_dark(outfile,idfile,auxfile);
    				free_memory();
    				fprintf(stderr,"%d ",ifile);
  			}
  			fprintf(stderr,"\n");
  		}
  		else fprintf(stderr,"No dark particles output\n");

  		if( bStar ) {
  			fprintf(stderr,"Outputting star particles in file number ");
  			for( ifile=i*files/Nout; ifile<(i+1)*files/Nout; ifile++ ) {
				if( ifile >= files ) continue;
    				if(files>1) sprintf(buf,"%s.%d",input_fname,ifile);
    				else sprintf(buf,"%s",input_fname);
    				load_snapshot(buf, type=4);
				unit_conversion(type=4);
    				output_tipsy_star(outfile,idfile,auxfile);
    				free_memory();
    				fprintf(stderr,"%d ",ifile);
  			}
  			fprintf(stderr,"\n");
  		}
		else fprintf(stderr,"No star particles output\n");

  /* fix header if gas particles not output due to density/temp criterion */
  		if( header.nsph != ngasout ) {
			fprintf(stderr,"Fixing header: %d != %d\n",header.nsph,ngasout);
			header.nsph = ngasout;
  			header.nbodies = header.nsph+header.ndark+header.nstar;
			rewind(outfile);
  			fwrite(&header, sizeof(header), 1, outfile);
		}
		fclose(outfile);
		fclose(idfile);
		fclose(auxfile);
	}

  	exit(0);
}



int output_tipsy_gas(FILE *outfile, FILE *idfile)
{
  int   i;

  for(i=1; i<=NumPart; i+=nskip) {
	  gasp.mass =   P[i].Mass;
	  /*      printf("gas mass %f ...\n",P[i].Mass); */
	  gasp.pos[0] = P[i].Pos[0];
	  gasp.pos[1] = P[i].Pos[1];
	  gasp.pos[2] = P[i].Pos[2];
	  gasp.vel[0] = P[i].Vel[0];
	  gasp.vel[1] = P[i].Vel[1];
	  gasp.vel[2] = P[i].Vel[2];
	  gasp.temp =   P[i].Temp;
	  gasp.hsmooth = 0.5*P[i].Hsml;	/* X2 difference in hsm definitions */
	  gasp.rho = P[i].Rho;
	  gasp.metals[0] = P[i].metal[0];
	  gasp.metals[1] = P[i].metal[1];
	  gasp.metals[2] = P[i].metal[2];
	  gasp.metals[3] = P[i].metal[3];
	  gasp.phi = P[i].Sfr;	/* store SFR (units of Mo/yr) in phi */
	  //	  gasp.delaytime = P[i].DelayTime;
	  //	  gasp.tmax = P[i].Tmax;
          gasp.delaytime = 0;
          gasp.tmax = 0;
	  if( gasp.rho < mindens || gasp.rho > maxdens || gasp.temp < minT || gasp.temp > maxT ) continue;
	  
	  //printf("%8d % 5.3e % 5.3e % 5.3e % 5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n",Id[i],gasp.pos[0],gasp.pos[1],gasp.pos[2],gasp.rho,gasp.temp,gasp.hsmooth,gasp.metals[0],gasp.metals[1],gasp.metals[2],gasp.metals[3],gasp.tmax,gasp.delaytime);
	  fwrite(&gasp, sizeof(struct gas_particle), 1, outfile);
	  if( id_flag ) fwrite(&Id[i], sizeof(int), 1, idfile) ;
	  ngasout++ ;
    }

  return 0;
}

int output_tipsy_dark(FILE *outfile, FILE *idfile)
{
  int   i;

  for(i=1; i<=NumPart; i+=nskip) {
	  darkp.mass =   P[i].Mass;
	  darkp.pos[0] = P[i].Pos[0];
	  darkp.pos[1] = P[i].Pos[1];
	  darkp.pos[2] = P[i].Pos[2];
	  darkp.vel[0] = P[i].Vel[0];
	  darkp.vel[1] = P[i].Vel[1];
	  darkp.vel[2] = P[i].Vel[2];
	  darkp.eps = epsilon;
	  darkp.phi = 0.;
	  
	  fwrite(&darkp, sizeof(struct dark_particle), 1, outfile) ;
	  if( id_flag ) fwrite(&Id[i], sizeof(int), 1, idfile) ;
    }

  return 0;
}

int output_tipsy_star(FILE *outfile, FILE *idfile)
{
  int   i;

  for(i=1; i<=NumPart; i+=nskip) {
	  starp.mass =   P[i].Mass;
	  /*	      printf("star mass %f ...\n",P[i].Mass); */
	  starp.pos[0] = P[i].Pos[0];
	  starp.pos[1] = P[i].Pos[1];
	  starp.pos[2] = P[i].Pos[2];
	  starp.vel[0] = P[i].Vel[0];
	  starp.vel[1] = P[i].Vel[1];
	  starp.vel[2] = P[i].Vel[2];
	  starp.metals[0] = P[i].metal[0];
	  starp.metals[1] = P[i].metal[1];
	  starp.metals[2] = P[i].metal[2];
	  starp.metals[3] = P[i].metal[3];
	  starp.tform = P[i].age;	/* expansion factor at formation */
	  starp.eps = epsilon;
	  starp.phi = 0.0;

	  fwrite(&starp, sizeof(struct star_particle), 1, outfile) ;
	  if( id_flag ) fwrite(&Id[i], sizeof(int), 1, idfile) ;
    }

  return 0;
}

