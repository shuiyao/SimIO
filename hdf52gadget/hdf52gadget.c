#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "gadgetdefs.h"
#include <math.h>
#include <hdf5.h>

static long NumPart;
struct particle_data *P;

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

int main(int argc, char **argv)
{
  int load_hdf5(char *snap);
  int WriteData(char *snap);  

  load_hdf5(argv[1]);

  WriteData(argv[1]);

  return(0);
}

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

  struct gadget_dump h;

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

#ifdef PHEW    
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWMcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Mcloud = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWRcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Rcloud = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWLastSFTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].LastSFTime = single[cnt];
      cnt += 1;
    }
#endif    

    H5Gclose(hdf5_grp);
    free(metals);
    free(single);
    free(posvel);
    ngas += h.npart[0];

    // DARK
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

int WriteData(char *snap)
{
  long i;
  long blksize;
  long ntot, noffset;

  char filename[256];

  FILE *fout;

  sprintf(filename,"%s.gad", snap);  
  fout = fopen(filename, "w");
#define SKIP  {fwrite(&blksize,sizeof(int),1,fout);}

  /* Header: */
  blksize = 256;
  SKIP;
  /* for(i=0;i<6;i++) io_header.npart[i] = 0; */
  /* io_header.npart[0] = header.nsph; */
  /* io_header.npart[1] = header.ndark; */
  /* io_header.npart[4] = header.nstar; */
  /* for(i=0;i<6;i++) io_header.npartTotal[i] = io_header.npart[i]; */
  /* for(i=0;i<6;i++) io_header.mass[i] = 0.0; */
  /* io_header.mass[1] = dp[0].mass * unit_Mass; */
  /* io_header.time = header.time; */
  /* io_header.redshift = 1. / header.time - 1.; */
  /* io_header.flag_sfr = 1; */
  /* io_header.flag_feedback = 1; */
  /* io_header.flag_cooling = 1; */
  /* io_header.num_files = 1; */
  /* io_header.BoxSize = BOXSIZE; */
  /* io_header.Omega0 = OMEGA0; */
  /* io_header.OmegaLambda = OMEGAL; */
  /* io_header.HubbleParam = HUBBLEPARAM; */
  gheader.mass[1] = P[0].Mass;
  fwrite(&gheader, sizeof(gheader), 1, fout);

  ntot = NumPart;
  
  SKIP;
  /* POS */
  blksize = ntot * 3 * sizeof(float);
  SKIP;
  for(i=0;i<ntot;i++){
    fwrite(P[i].Pos, sizeof(float), 3, fout);
  }
  /* for(i=0;i<header.nstar;i++){ */
  /*   for(k=0;k<3;k++) */
  /*     fp[k] = sp[i].pos[k] * unit_Length; */
  /*   fwrite(fp, sizeof(float), 3, fout); */
  /* } */
  SKIP;

  /* VEL */
  blksize = ntot * 3 * sizeof(float);
  SKIP;
  for(i=0;i<ntot;i++){
    fwrite(P[i].Vel, sizeof(float), 3, fout);
  }
  SKIP;

  /* ID */
  blksize = ntot * sizeof(int);
  SKIP;
  for(i=0;i<ntot;i++){
    fwrite(&P[i].ID, sizeof(int), 1, fout);
  }
  SKIP;

  /* MASS */
  blksize = (gheader.npart[0] + gheader.npart[4]) * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++){
    fwrite(&P[i].Mass, sizeof(float), 1, fout);
  }
  for(i=0;i<gheader.npart[4];i++){
    fwrite(&P[i].Mass, sizeof(float), 1, fout);
  }
  SKIP;

  /* U */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++)
    fwrite(&P[i].Temp, sizeof(float), 1, fout);
  SKIP;

  /* RHO */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++){
    fwrite(&P[i].Rho, sizeof(float), 1, fout);
  }
  SKIP;

  /* NE */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++)
    fwrite(&P[i].Ne, sizeof(float), 1, fout);
  SKIP;

  /* NH*/
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++)
    fwrite(&P[i].Nh, sizeof(float), 1, fout);
  SKIP;

  /* HSML */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++){
    fwrite(&P[i].Hsml, sizeof(float), 1, fout);
  }
  SKIP;

  /* SFR */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++)
    fwrite(&P[i].Sfr, sizeof(float), 1, fout);
  SKIP;

  /* DELAYT */
  blksize = gheader.npart[0] * sizeof(float);
  SKIP;
  for(i=0;i<gheader.npart[0];i++)
    fwrite(&P[i].DelayTime, sizeof(float), 1, fout);
  SKIP;

  /* RHO */
  blksize = gheader.npart[4] * sizeof(float);
  SKIP;
  noffset = gheader.npart[0] + gheader.npart[1];
  for(i=noffset;i<gheader.npart[4]+noffset;i++)
    fwrite(&P[i].Sfr, sizeof(float), 1, fout);
  SKIP;
  fclose(fout);
  return 1;
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
