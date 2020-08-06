/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

#include <stdio.h>
#include "gadgetdefs.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int flag_rho = 0;
int flag_temp = 0;
int flag_id = 0;
int flag_phew = 0;

int flag_read_dark = 0;
int flag_read_star = 0;

/* Unit Conversion for .tipsy */
/*   - Read in .hdf file, write it into .tipsy file. */
/*   - Usage: hdf52tipsy snap_* (no .hdf5) */

/* WARNING: */
/* Tipsy metallicity array contains only C, O, Si, Fe. Make sure the corresponding field is loaded from the HDF5 output. Default: [C, O, Si, Fe] = [2, 4, 7, 10] (NMETALS = 11) */

static long NumPart;

double unit_Time;
double unit_Length;
double unit_Density;
double unit_Mass;
double unit_Velocity;
double unit_Temp;

struct particle_data *P;
struct gadget_dump h;

char basename[200];

int main(int argc, char **argv){  
  int write_hdf5_fields_to_ascii(char *snap);
  int load_hdf5(char *snap);
  int parse_input(int argc, char **argv, char *basename);

  /* Load HDF5 snapshot */
  /* Inherited from Romeel's specexsnap code */
  /* Should return 1. a gadget header h; 2. Particle data P */
  printf("We have %d parameters here, master.\n", argc);
  strcpy(basename, argv[1]);
  parse_input(argc, argv, basename);
  
  load_hdf5(basename);

  write_hdf5_fields_to_ascii(basename);
  /* write_ascii(argv[1]); */
  return 1;
}

int write_hdf5_fields_to_ascii(char *snap){
  /* Pass through flags, write for each flag */
  int i;
  char filename[200];
  FILE *outputfile;
  
  if(flag_rho) {
    sprintf(filename, "%s.rho", snap);
    outputfile = fopen(filename, "w");
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].Rho);
    }
    fclose(outputfile);
  }
  if(flag_temp) {
    sprintf(filename, "%s.temp", snap);        
    outputfile = fopen(filename, "w");
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].Temp);
    }
    fclose(outputfile);
  }
  if(flag_id) {
    sprintf(filename, "%s.id", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%d\n", P[i].ID);
    }
    fclose(outputfile);
  } 
  if(flag_phew) {
    sprintf(filename, "%s.key", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%d\n", P[i].Key);
    }
    fclose(outputfile);
  } 
  if(flag_phew) {
    sprintf(filename, "%s.mcloud", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].Mcloud);
    }
    fclose(outputfile);
  } 
  if(flag_phew) {
    sprintf(filename, "%s.windmass", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].WindMass);
    }
    fclose(outputfile);
  } 
  if(flag_phew) {
    sprintf(filename, "%s.lastsftime", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].LastSFTime);
    }
    fclose(outputfile);
  } 
  if(flag_phew) {
    sprintf(filename, "%s.vinit", snap);        
    outputfile = fopen(filename, "w");
    fprintf(outputfile, "%d\n", h.npart[0]);
    for(i=0; i < h.npart[0]; i++) {
      fprintf(outputfile,"%g\n", P[i].Vinit);
    }
    fclose(outputfile);
  }
  return 1; 
} // int write_hdf5_fields_to_ascii(char *snap)

int load_hdf5(char *snap){

  char infile[256];
  FILE *fp;
  int multipart = 0;

  int allocate_memory();
  
  sprintf(infile,"%s.hdf5",snap);
  if(!(fp=fopen(infile,"r"))){
    sprintf(infile,"%s.0.hdf5",snap);
    multipart = 1;
    if(!(fp=fopen(infile,"r"))){
      fprintf(stderr,"Error opening file '%s' \n",infile);
      exit(0);
    }
  }
  fclose(fp);

  long i, cnt;
  int j, k;
  long noffset;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute, hdf5_grp, hdf5_dataset;
  /* int ngas  = 1; // Why? */
  long ngas  = 0; 
  long ndark  = 0; 
  long nstar  = 0; 
  float *posvel, *single, *metals;
  int *intsingle;

  hdf5_file = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.redshift);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.Omega0);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, h.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Metals");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_metals);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_StellarAge");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_stellarage);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_cooling);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Feedback");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_feedback);
  H5Aclose(hdf5_attribute);

  /* Missing several fields including flag_sfr, flag_feedback, nparttotal, etc. */

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  NumPart = h.npartTotal[0] + h.npartTotal[1] + h.npartTotal[4];

  h.flag_metals = NMETALS;

  gheader = h;

  fprintf(stderr,"NumFiles = %d\n", h.num_files);  
  fprintf(stderr,"NumPart = %ld\n", NumPart);
  fprintf(stderr,"Time = %g; Redshift = %g\n", h.time, h.redshift);  
  
  allocate_memory();

  for(k=0; k<h.num_files; k++){
    if(multipart)
      sprintf(infile,"%s.%d.hdf5",snap,k);

    hdf5_file = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(multipart){
      hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npart);
      H5Aclose(hdf5_attribute);
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
      H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, h.mass);
      H5Aclose(hdf5_attribute);
      H5Gclose(hdf5_headergrp);
    }

    //GAS
    fprintf(stdout, "Reading GAS for file #%d\n", k);    
    if(!(posvel = (float *)malloc(sizeof(float)*h.npart[0]*3))){
      fprintf(stderr, "Failed to allocate memory for posvel\n");
      exit(-1);
    }
    if(!(single = (float *)malloc(sizeof(float)*h.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for single\n");
      exit(-1);
    }
    if(!(metals = (float *)malloc(sizeof(float)*h.npart[0]*h.flag_metals))){
      fprintf(stderr, "Failed to allocate memory for single\n");
      exit(-1);
    }
    if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for intsingle\n");
      exit(-1);
    }
    
    hdf5_grp = H5Gopen1(hdf5_file, "/PartType0");
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Pos[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Vel[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].ID = intsingle[cnt];
      cnt += 1;
    }

    if(h.mass[0] == 0 && h.npart[0] > 0){
      hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
	P[i].Mass = single[cnt];
	cnt += 1;
      }
    }
    else{
      for(i=ngas; i<h.npart[0]+ngas; i++)
	P[i].Mass = h.mass[0];
    }
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "InternalEnergy");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Temp = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Density");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Rho = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Ne = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "NeutralHydrogenAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Nh = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "SmoothingLength");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Hsml = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "StarFormationRate");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++) {
      P[i].Sfr = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Tmax = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "DelayTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].DelayTime = single[cnt];
      cnt += 1;
    }

    if(H5Lexists(hdf5_grp,"FractionH2",H5P_DEFAULT)==0){
      for(i=ngas; i<h.npart[0]+ngas; i++){
	P[i].fH2 = 0.0;
      }
    }
    else{
      hdf5_dataset = H5Dopen1(hdf5_grp, "FractionH2");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
	P[i].fH2 = single[cnt];
	cnt += 1;
      }
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<h.flag_metals; j++)
	P[i].metal[j] = metals[cnt*h.flag_metals + j];
      cnt += 1;
    }

#ifdef PHEW_EXTRA_OUTPUT
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWKey");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Key = intsingle[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWMcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Mcloud = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWWindMass");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].WindMass = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWLastSFTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].LastSFTime = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWVinit");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Vinit = single[cnt];
      cnt += 1;
    }
#endif    

    H5Gclose(hdf5_grp);
    free(metals);
    free(single);
    free(posvel);
    ngas += h.npart[0];

    // DARK
    if(flag_read_dark){
      fprintf(stdout, "Reading DARK from file #%d\n", k);    
      noffset = gheader.npartTotal[0];
      if(!(posvel = (float *)malloc(sizeof(float)*h.npart[1]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*h.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType1");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(h.mass[1] == 0 && h.npart[1] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+ndark; i<noffset+h.npart[1]+ndark; i++)
	  P[i].Mass = h.mass[1];

      H5Gclose(hdf5_grp);
      free(single);
      free(posvel);
      ndark += h.npart[1];
    } // flag_read_dark

    if(flag_read_star){
      fprintf(stdout, "Reading STAR from file #%d\n", k);
      // STAR
      noffset = gheader.npartTotal[0] + gheader.npartTotal[1];
      if(!(posvel = (float *)malloc(sizeof(float)*h.npart[4]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*h.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(metals = (float *)malloc(sizeof(float)*h.npart[4]*h.flag_metals))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType4");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(h.mass[4] == 0 && h.npart[4] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+nstar; i<noffset+h.npart[4]+nstar; i++)
	  P[i].Mass = h.mass[4];

      hdf5_dataset = H5Dopen1(hdf5_grp, "StellarFormationTime");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].Sfr = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].Tmax = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<h.flag_metals; j++)
	  P[i].metal[j] = metals[cnt*h.flag_metals + j];
	cnt += 1;
      }
    
      H5Gclose(hdf5_grp);
      free(posvel);
      free(single);
      free(metals);
      free(intsingle);
      nstar += h.npart[4];
    } // flag_read_star

    H5Fclose(hdf5_file);
    fprintf(stderr, "File: %d ngas = %d(%5.3f)\n", k, h.npart[0],
	    (float)(h.npart[0])/(float)(NumPart));
    fprintf(stderr, "File: %d ndark = %d(%5.3f)\n", k, h.npart[1],
	    (float)(h.npart[1])/(float)(NumPart));
    fprintf(stderr, "File: %d nstar = %d(%5.3f)\n", k, h.npart[4],
	    (float)(h.npart[4])/(float)(NumPart));
  }
  return 0;
}

int allocate_memory(void)
{
  fprintf(stdout, "Allocating %6.3f GB Memory for %ld particles.\n",
	  NumPart * sizeof(struct particle_data) / (1024. * 1024. * 1024.),
	  NumPart);

  if(!(P=malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

  fprintf(stdout, "Memory allocated.\n");
  /* P--;    */
  /* start with offset 1 */
  return 0;
}

int parse_input(int argc, char **argv, char *basename)
{
  int i=1;
  while(i < argc) {
    printf("%s\n", argv[i]);
    if (!strcmp(argv[i], "-rho"))
      {flag_rho=1; ++i;}
    else if (!strcmp(argv[i], "-temp"))
      {flag_temp=1; ++i;}
    else if (!strcmp(argv[i], "-id"))
      {flag_id=1; ++i;}
    else if (!strcmp(argv[i], "-phew"))
      {flag_phew=1; ++i;}
    else {
      strcpy(basename, argv[i]); ++i;
    }
  }
  return 0;
}
