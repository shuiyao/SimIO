/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

#include <stdio.h>
#include "gadgetdefs.h"
#include "tipsydefs.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
struct dump header;

void cosmounits(void)
{
  double HubbleParam, BoxSize;
  HubbleParam = gheader.HubbleParam;
  BoxSize = gheader.BoxSize;
  unit_Time=sqrt(8*M_PI/3)*CM_PER_MPC/(100*HubbleParam*1.e5);
  unit_Density=1.8791E-29*HubbleParam*HubbleParam;
  unit_Length=BoxSize*CM_PER_MPC*1.e-3;
  unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length/(HubbleParam*HubbleParam);
  unit_Velocity=unit_Length/unit_Time;
  unit_Temp = pow(UNIT_V, 2);
  fprintf(stdout, "unit_Length = %g\n", unit_Length);
  fprintf(stdout, "unit_Mass = %g\n", unit_Mass);
  fprintf(stdout, "unit_Velocity = %g\n", unit_Velocity);
  // = pow(UNIT_L, 2) / pow((UNIT_L/Unit_V), 2); 
  return;
}

int main(int argc, char **argv){  
  int load_hdf5(char *snap);
  int write_tipsy(char *snap);

  /* Load HDF5 snapshot */
  /* Inherited from Romeel's specexsnap code */
  /* Should return 1. a gadget header h; 2. Particle data P */
  load_hdf5(argv[1]);

  cosmounits(); // when gheader is at place

  /* Write into tipsy files */
  write_tipsy(argv[1]);
  return 1;
}

/* Write the tipsy file. */
/* Currently only supports GAS particles. */
int write_tipsy(char *snap){
  
  int k;
  long i, pidx, noffset;
  char binfile[256], auxfile[256], idnumfile[256];

  struct gas_particle *gp;
  struct aux_gas_data *auxgp;
  struct star_particle *sp;
  struct aux_star_data *auxsp;
  struct dark_particle *dp;
#ifdef PHEW
  struct aw_gas_data *awgp;
#endif  
  int *pids;
  FILE *fout;
  double MeanWeight;
  /* double a3inv; */
  /* double tmp, XH, XHe; */

  sprintf(binfile,"%s.bin",snap);
  sprintf(auxfile,"%s.aux",snap);
  sprintf(idnumfile,"%s.idnum",snap);
#ifdef PHEW
  char awfile[256];
  sprintf(awfile,"%s.aw",snap);
#endif  

  theader.time = gheader.time;
  theader.ndim = 3;
  theader.nsph = gheader.npartTotal[0];
  theader.ndark = gheader.npartTotal[1];
  theader.nstar = gheader.npartTotal[4];
  theader.nbodies = theader.nsph + theader.ndark + theader.nstar;

  if(!(pids=malloc(theader.nbodies*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory for pids.\n");
      exit(0);
    }

  // Allocate memory for gp, auxgp
  if(!(gp=malloc(gheader.npartTotal[0]*sizeof(struct gas_particle))))
    {
      fprintf(stderr,"failed to allocate memory for gp.\n");
      exit(0);
    }
  if(!(auxgp=malloc(gheader.npartTotal[0]*sizeof(struct aux_gas_data))))
    {
      fprintf(stderr,"failed to allocate memory for auxgp.\n");
      exit(0);
    }
  if(!(dp=malloc(gheader.npartTotal[1]*sizeof(struct dark_particle))))
    {
      fprintf(stderr,"failed to allocate memory for dp.\n");
      exit(0);
    }
  if(!(sp=malloc(gheader.npartTotal[4]*sizeof(struct star_particle))))
    {
      fprintf(stderr,"failed to allocate memory for sp.\n");
      exit(0);
    }
  if(!(auxsp=malloc(gheader.npartTotal[4]*sizeof(struct aux_star_data))))
    {
      fprintf(stderr,"failed to allocate memory for auxsp.\n");
      exit(0);
    }
#ifdef PHEW  
  if(!(awgp=malloc(gheader.npartTotal[0]*sizeof(struct aw_gas_data))))
    {
      fprintf(stderr,"failed to allocate memory for awgp.\n");
      exit(0);
    }
#endif

  /* if(gheader.time != 0){ // Cosmos Run */
  /*   a3inv = 1. / (gheader.time * gheader.time * gheader.time); */
  /* } */
  /* else{ // Non-cosmos Run */
  /*   a3inv = 1.; */
  /*   gheader.time = 1.; */
  /* } */

  fprintf(stderr, "Writing Tipsy: Ngas = %d\n", gheader.npartTotal[0]);
  for(i=0;i<gheader.npartTotal[0];i++){
    pids[i] = P[i].ID;
    gp[i].mass = P[i].Mass * UNIT_M / unit_Mass;
    for(k=0;k<3;k++){
      gp[i].pos[k] = P[i].Pos[k] * UNIT_L / unit_Length - 0.5;
      if(gp[i].pos[k] >= 0.5) gp[i].pos[k] = 0.499999;
      if(gp[i].pos[k] <= -0.5) gp[i].pos[k] = -0.499999;				
      gp[i].vel[k] = P[i].Vel[k] * UNIT_V / unit_Velocity / sqrt(gheader.time);
      gp[i].vel[k] *= sqrt(gheader.time * gheader.time * gheader.time);
      // V.hdf5 = V.gad * sqrt(a3inv); see io.c
      }
    gp[i].phi = 0.0;
    gp[i].rho = P[i].Rho * UNIT_M / (pow(UNIT_L, 3) * unit_Density);
    gp[i].hsmooth = P[i].Hsml * UNIT_L / unit_Length;

    gp[i].metals = P[i].metal[0]; 

    /* gp[i].metals += P[i].metal[2]; // C */
    /* gp[i].metals += P[i].metal[4]; // O */
    /* gp[i].metals += P[i].metal[7]; // Si */
    /* gp[i].metals += P[i].metal[10]; // Fe */
    /* gp[i].metals *= 1.28571; */
    // elements: He, C, Mg, O, Fe, Si, H, N, Ne, S, Ca

    /* XH = 1.0 - (P[i].metal[0] + P[i].metal[1]); */
    /* XHe = P[i].metal[1]; */
    /* tmp = 0.0; */
    /* tmp += P[i].metal[2] / 12.011; // C */
    /* tmp += P[i].metal[4] / 15.999; // O */
    /* tmp += P[i].metal[7] / 28.085; // Si */
    /* tmp += P[i].metal[10] / 55.845; // Fe */
    /* tmp *= 1.28571; */
    /* tmp += XH / 1.008; */
    /* tmp += XHe / 4.003; */
    /* MeanWeight = ATOMIC_UNIT / (tmp + P[i].Ne * XH); */
    /* MeanWeight= 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC*P[i].Ne)*PROTONMASS; */
    /* gp[i].temp = P[i].Temp / GAMMA_MINUS1 * pow(P[i].Rho * a3inv, GAMMA_MINUS1); */
    /* gp[i].temp = P[i].Temp * MeanWeight / BOLTZMANN * GAMMA_MINUS1 * unit_Temp; */
    /* gp[i].temp *= (gheader.time * gheader.time); // The unit of U contains this. */

    /* From GIZMO:cooling/cooling.c */
    /* static double mhboltz = PROTONMASS / BOLTZMANN; */          
    /* static double yhelium = (1 - HYDROGEN_MASSFRAC) / (4 * HYDROGEN_MASSFRAC);       */
    /* T = log10(GAMMA_MINUS1 * u * mhboltz * (1 + 4 * yhelium) / (1 + ne + yhelium));   */

    MeanWeight = (1 + 4 * XHE) / (1 + P[i].Ne + XHE);
    gp[i].temp = P[i].Temp * unit_Temp;
    gp[i].temp *= GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * MeanWeight;

    auxgp[i].metal[0] = P[i].metal[2];
    auxgp[i].metal[1] = P[i].metal[4];
    auxgp[i].metal[2] = P[i].metal[7];
    auxgp[i].metal[3] = P[i].metal[10];
    auxgp[i].sfr = P[i].Sfr;
    auxgp[i].tmax = P[i].Tmax;
    auxgp[i].delaytime = P[i].DelayTime;
    auxgp[i].ne = P[i].Ne;
    auxgp[i].nh = P[i].Nh;
    auxgp[i].nspawn = 0;
    auxgp[i].nrec = 0;

#ifdef PHEW
    awgp[i].rho = gp[i].rho;
    awgp[i].temp = gp[i].temp;
    awgp[i].mcloud = P[i].Mcloud; // fractional
    awgp[i].rcloud = P[i].Rcloud * UNIT_L / unit_Length;
    awgp[i].lastsftime = P[i].LastSFTime;    
    awgp[i].wind_flag = (P[i].Mcloud > 0) ? 1 : 0;
#endif    
  }
  // DARK
  fprintf(stderr, "Writing Tipsy: Ndark = %d\n", gheader.npartTotal[1]);
  noffset = gheader.npartTotal[0];
  for(i=0;i<gheader.npartTotal[1];i++){
    pidx = i + noffset;
    pids[pidx] = P[pidx].ID;
    dp[i].mass = P[pidx].Mass * UNIT_M / unit_Mass;
    for(k=0;k<3;k++){
      dp[i].pos[k] = P[pidx].Pos[k] * UNIT_L / unit_Length - 0.5;
      if(dp[i].pos[k] >= 0.5) dp[i].pos[k] = 0.499999;
      if(dp[i].pos[k] <= -0.5) dp[i].pos[k] = -0.499999;				
      dp[i].vel[k] = P[pidx].Vel[k] * UNIT_V / unit_Velocity / sqrt(gheader.time);
      dp[i].vel[k] *= sqrt(gheader.time * gheader.time * gheader.time);
      // V.hdf5 = V.gad * sqrt(a3inv); see io.c
      }
    dp[i].phi = 0.0;
    dp[i].eps = 0.0;
  }
  // STAR
  fprintf(stderr, "Writing Tipsy: Nstar = %d\n", gheader.npartTotal[4]);
  noffset = gheader.npartTotal[0] + gheader.npartTotal[1];
  for(i=0;i<gheader.npartTotal[4];i++){
    pidx = i + noffset;
    pids[pidx] = P[pidx].ID;    
    sp[i].mass = P[pidx].Mass * UNIT_M / unit_Mass;
    for(k=0;k<3;k++){
      sp[i].pos[k] = P[pidx].Pos[k] * UNIT_L / unit_Length - 0.5;
      if(sp[i].pos[k] >= 0.5) sp[i].pos[k] = 0.499999;
      if(sp[i].pos[k] <= -0.5) sp[i].pos[k] = -0.499999;				
      sp[i].vel[k] = P[pidx].Vel[k] * UNIT_V / unit_Velocity / sqrt(gheader.time);
      sp[i].vel[k] *= sqrt(gheader.time * gheader.time * gheader.time);
      // V.hdf5 = V.gad * sqrt(a3inv); see io.c
      }
    sp[i].phi = 0.0;
    sp[i].eps = 0.0;
    sp[i].tform = 0.0;
    sp[i].metals = 0.0;
    sp[i].metals += P[pidx].metal[2]; // C
    sp[i].metals += P[pidx].metal[4]; // O
    sp[i].metals += P[pidx].metal[7]; // Si
    sp[i].metals += P[pidx].metal[10]; // Fe
    sp[i].metals *= 1.28571;
    // elements: He, C, Mg, O, Fe, Si, H, N, Ne, S, Ca

    auxsp[i].metal[0] = P[pidx].metal[2];
    auxsp[i].metal[1] = P[pidx].metal[4];
    auxsp[i].metal[2] = P[pidx].metal[7];
    auxsp[i].metal[3] = P[pidx].metal[10];
    auxsp[i].tmax = P[pidx].Tmax;
    auxsp[i].age = P[pidx].Sfr;
    auxsp[i].nspawn = 0;
    auxsp[i].nrec = 0;
  }

  fout = fopen(binfile, "w");
  fwrite(&theader, sizeof(theader), 1, fout);
  fwrite(gp, sizeof(struct gas_particle), gheader.npartTotal[0], fout);
  fwrite(dp, sizeof(struct dark_particle), gheader.npartTotal[1], fout);
  fwrite(sp, sizeof(struct star_particle), gheader.npartTotal[4], fout);
  fclose(fout);

  fout = fopen(auxfile, "w");
  fwrite(auxgp, sizeof(struct aux_gas_data), gheader.npartTotal[0], fout);
  fwrite(auxsp, sizeof(struct aux_star_data), gheader.npartTotal[4], fout);  
  fclose(fout);

  fout = fopen(idnumfile, "w");
  fwrite(pids, sizeof(int), theader.nbodies, fout);
  fclose(fout);

#ifdef PHEW  
  fout = fopen(awfile, "w");
  fwrite(awgp, sizeof(struct aw_gas_data), gheader.npartTotal[0], fout);
  fclose(fout);
  free(awgp);
#endif  

  free(auxsp);
  free(sp);
  free(dp);  
  free(auxgp);
  free(gp);
  free(pids);
  return 0;
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
