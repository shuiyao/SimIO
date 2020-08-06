#include <stdio.h>
#include "gadgetdefs.h"
#include "tipsydefs.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* SH2019 */

/* Unit Conversion for .tipsy */
/*   - Read in .tipsy file, write it into .hdf5 file. */
/*   - Usage: tipsy2hdf5 snap_* (no suffix) */

/* WARNING: */
/* Tipsy metallicity array contains only C, O, Si, Fe. Make sure the corresponding field is loaded from the HDF5 output. Default: [C, O, Si, Fe] = [2, 4, 7, 10] (NMETALS = 11) */

double unit_Time;
double unit_Length;
double unit_Density;
double unit_Mass;
double unit_Velocity;
double unit_Temp;

struct dump header;
struct gadget_dump gheader ;
struct gas_particle *gp;
struct dark_particle *dp;
struct star_particle *sp;
struct aux_gas_data *aux_gp;
struct aux_star_data *aux_sp;
struct aw_gas_data *aw_gp;
struct block_attributes *blkattr;

#define HYDROGEN_MASSFRAC 0.76
#define PROTONMASS 1.6726e-24

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
  // = pow(UNIT_L, 2) / pow((UNIT_L/Unit_V), 2); 
  return;
}

int main(int argc, char **argv){  
  int read_tipsy_binary(char *snap);/* OUTPUT: struct (*_particle) *gp, *dp, *sp */
  int read_tipsy_aux(char *snap);/* OUTPUT: struct aux_gas_data *aux_gp */
  int create_gadget_header(int nfiles);
  int write_hdf5(char *snap, int nfiles);
  void init_block_info();

  /* Initialize the HDF5 blocks */
  init_block_info();
  
  /* Load data from TIPSY binaries */ 
  read_tipsy_binary(argv[1]);
  read_tipsy_aux(argv[1]);

  /* Create a fake gadget header */
  if(argc > 2) create_gadget_header(atoi(argv[2]));
  else create_gadget_header(1);

  cosmounits(); // when gheader is at place

  /* Write into HDF5 files */
  if(argc > 2)
    write_hdf5(argv[1], atoi(argv[2])); /* argv[2] Files */
  else
    write_hdf5(argv[1], 1); /* One File */

  fprintf(stdout, "SUCCESS, YAY!\n");
  return 1;
}

int create_gadget_header(int nfiles){
  int k;
  for(k = 0; k < 6; k ++) gheader.npart[k] = 0;
  gheader.npart[0] = header.nsph;
  gheader.npart[1] = header.ndark;
  gheader.npart[4] = header.nstar;
  for(k = 0; k < 6; k ++){
    gheader.npartTotal[k] = gheader.npart[k];
    gheader.npartTotalHighWord[k] = 0.0;
    gheader.mass[k] = 0;
  }
  gheader.time = header.time;
  gheader.redshift = 1. / header.time - 1.;

  gheader.flag_sfr = 1;
  gheader.flag_feedback = 1;
  gheader.flag_cooling = 1;
  gheader.flag_stellarage = 1;
  gheader.flag_metals = 1;  
  gheader.flag_tmax = 1;
  gheader.flag_delaytime = 1;
  gheader.flag_potential = 1;
  gheader.flag_fH2 = 0;

  /* Here are values for the cluster 300 run. */
  /* gheader.BoxSize = 1.e6; */
  /* gheader.Omega0 = 0.307115; */
  /* gheader.OmegaLambda = 0.692885; */
  /* gheader.HubbleParam = 0.6777; */
  /* Values for the p6n36 series. */
  gheader.BoxSize = BOXSIZE;
  gheader.Omega0 = OMEGA_0;
  gheader.OmegaLambda = OMEGA_L;
  gheader.HubbleParam = HUBBLE_PARAM;

  gheader.num_files = nfiles;
  return 1;
}

int read_tipsy_binary(char *snap)
{
  char filename[256];
  FILE *fp;

  sprintf(filename, "%s.bin", snap);
  fp = fopen(filename, "r");
  
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&header, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&header,sizeof(header),1,fp) == 0)
    {printf("BREAK!\n"); exit(-1);}
#endif
  //  headerinfo(&header);
  printf("*** Header Information ***\n");
  printf("header.time: %f\n", header.time);
  printf("header.nbodies: %d\n", header.nbodies);
  printf("header.ndim: %d\n", header.ndim);
  printf("header.nsph: %d\n", header.nsph);
  printf("header.ndark: %d\n", header.ndark);
  printf("header.nstar: %d\n", header.nstar);
  printf("**************************\n");
  if(header.nsph != 0) {
    gp = (struct gas_particle *)
      malloc(header.nsph*sizeof(*gp));
    if(gp == NULL) {
      printf("<sorry, no memory for gas particles, master>\n") ;
      return -1;
    }
  }
  if(header.ndark != 0) {
    dp = (struct dark_particle *)
      malloc(header.ndark*sizeof(*dp));
    if(dp == NULL) {
      printf("<sorry, no memory for dark particles, master>\n") ;
      return -1;
    }
  }
  if(header.nstar != 0) {
    sp = (struct star_particle *)
      malloc(header.nstar*sizeof(*sp));
    if(sp == NULL) {
      printf("<sorry, no memory for star particles, master>\n") ;
      return -1;
    }
  }
  fread((char *)gp,sizeof(struct gas_particle),	header.nsph,fp) ;
  fread((char *)dp,sizeof(struct dark_particle), header.ndark,fp) ;
  fread((char *)sp,sizeof(struct star_particle), header.nstar,fp) ;
  fclose(fp);
  printf("---> Reading .bin done.\n");
  fprintf(stdout,"read time %f\n",header.time);
  return 1;
}

int read_tipsy_aux(char *snap)
{
/* Read in auxiliary particle data file */
  char filename[256];
  FILE *fp;

  sprintf(filename, "%s.aux", snap);
  fp = fopen(filename, "r");
  if(header.nsph > 0){
    aux_gp = (struct aux_gas_data *) malloc(header.nsph*sizeof(*aux_gp));
    fread(aux_gp,sizeof(struct aux_gas_data),header.nsph,fp);
  }
  if(header.nstar > 0){
    aux_sp = (struct aux_star_data *) malloc(header.nstar*sizeof(*aux_sp));
    fread(aux_sp,sizeof(*aux_sp),header.nstar,fp);
  }
  fclose(fp);
  fprintf(stdout, "---> Reading .aux done.\n");  
  return 1;
}

int read_tipsy_aw(char *snap)
{
/* Read in auxiliary particle data file */
  char filename[256];
  FILE *fp;

  sprintf(filename, "%s.aw", snap);
  fp = fopen(filename, "r");
  
  aw_gp = (struct aw_gas_data *) malloc(header.nsph*sizeof(*aw_gp));
  fread(aw_gp,sizeof(struct aw_gas_data),header.nsph,fp);
  fclose(fp);
  fprintf(stdout, "---> Reading .aw done.\n");  
  return 1;
}

int write_hdf5(char *snap, int nfiles){
  void write_hdf5_header_attributes(hid_t handle);
  
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0;  /* identifiers */
  hid_t hdf5_dataset = 0, hdf5_datatype = 0, hdf5_dataspace = 0;
  hsize_t dims[2];
  herr_t hdf5_status;
  int rank = 0;

  int i, k, type, blocknr;
  char hdf5_blkname[50];
  float *f1block, *f3block, *fZblock;
  int *iblock;
  short int *siblock;
  /* *float, *float[3], *int, *short int */

#ifdef INTERNAL_ENERGY
  float mu;
#endif  

  char filename[256], buf[50];

  sprintf(filename, "%s.hdf5", snap);

  /* Create a new file using default properties. */
  hdf5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  /* H5F_ACC_TRUNC means OVERWRITE, use H5F_ACC_EXCEL to check if existed */

  /* Create the dataset. */
  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
  /* H5Gcreate: The last parameter is size_hint, 0 means DEFAULT */   
  write_hdf5_header_attributes(hdf5_headergrp);
  /* Close the groups. */
  hdf5_status = H5Gclose(hdf5_headergrp);

  for(type = 0; type < 6; type++){
    fprintf(stdout, "gheader.npart[%d] = %d\n", type, gheader.npart[type]);
    if(gheader.npart[type] > 0){
      sprintf(buf, "/PartType%d", type);
      hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);

      for(blocknr = 0; blocknr < NBLOCK; blocknr ++){
	if(!((1 << type) & blkattr[blocknr].typeflag)) continue;
	/* This type do not have that block! e.g. type=1(DARK) do not have blocknr=7 (Density) */

	dims[0] = gheader.npart[type];
	dims[1] = blkattr[blocknr].dim;
	if(dims[1] == 1) rank = 1;
	else rank = 2;
	
	switch (blkattr[blocknr].dtype) {
#ifdef LONGIDS	  
	case 1:
	  hdf5_datatype = H5Tcopy(H5T_NATIVE_ULONG);
	  if((iblock = malloc(gheader.npart[type] * sizeof(long))) == NULL){
	    fprintf(stderr, "Failed to allocate memory for iblock. Exit.\n"); exit(-1);}
	  break;
#else	  
	case 1:
	  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
	  if((iblock = malloc(gheader.npart[type] * sizeof(int))) == NULL){
	    fprintf(stderr, "Failed to allocate memory for iblock. Exit.\n"); exit(-1);}
	  break;
#endif	  
	case 2:
	  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
	  if(dims[1] == 1){
	    if((f1block = malloc(gheader.npart[type] * sizeof(float))) == NULL){
	      fprintf(stderr, "Failed to allocate memory for f1block. Exit.\n"); exit(-1);}
	  }
	  else if(dims[1] == 3){
	    if((f3block = malloc(gheader.npart[type] * sizeof(float) * 3)) == NULL){
	      fprintf(stderr, "Failed to allocate memory for f3block. Exit.\n"); exit(-1);}
	  }
	  else if(dims[1] == NMETALS){
	    if((fZblock = malloc(gheader.npart[type] * sizeof(float) * NMETALS)) == NULL){
	      fprintf(stderr, "Failed to allocate memory for fZblock. Exit.\n"); exit(-1);}
	  }
	  break;
	case 3:
	  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
	  if((siblock = malloc(gheader.npart[type] * sizeof(short int))) == NULL){
	    fprintf(stderr, "Failed to allocate memory for siblock. Exit.\n"); exit(-1);}
	  break;
	/* case 4: */
	/*   hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE); */
	/*   if((dblock = malloc(gheader.npart[type] * sizeof(double))) == NULL) */
	/*     fprintf(stderr, "Failed to allocate memory for iblock. Exit.\n"); exit(-1); */
	/*   break;	   */

	} // switch

	/* Here: Fill the Buffer! */
	if(type == 0){
	  for (i = 0; i < gheader.npart[type]; i ++){
	    switch (blocknr){
	    case 0:
	      for(k=0; k<3; k++) {
		f3block[i*3+k] = (gp[i].pos[k] + 0.5) * unit_Length / UNIT_L;
		if(f3block[i*3+k] > gheader.BoxSize)
		  f3block[i*3+k] -= gheader.BoxSize;
		if(f3block[i*3+k] < 0)
		  f3block[i*3+k] += gheader.BoxSize;
	      }
	      break;
	    case 1:
	      for(k=0; k<3; k++) f3block[i*3+k] = gp[i].vel[k] * unit_Velocity / UNIT_V / gheader.time; break;
	      // V.hdf5 = V.gad * sqrt(a3inv); V.tipsy = V.gad / sqrt(a)
	    case 2: break; // ParticleID
	    case 3: f1block[i] = gp[i].mass * unit_Mass / UNIT_M; break;
	    case 4: f1block[i] = gp[i].phi * (unit_Velocity / UNIT_V) * (unit_Velocity / UNIT_V); break;
	    case 5:
#ifdef INTERNAL_ENERGY
	      mu = 4.0 / (3.*0.76 + 1. + 4.*0.76*aux_gp[i].ne) * 1.6726e-24;
	      f1block[i] = gp[i].temp * 2.070987e-16 / mu * 1.e-10; // 1.e-10 = unit_M / unit_E = 1.989e43 / 1.989e53, 2.07e-16 = k/(gamma-1)
	      break;
#else	      
	      f1block[i] = gp[i].temp; break; // Temperature Already !
#endif	      
	    case 6: f1block[i] = gp[i].rho * unit_Density * pow(UNIT_L, 3) / UNIT_M; break;
	    case 7: f1block[i] = aux_gp[i].ne; break;
	    case 8: f1block[i] = aux_gp[i].nh; break;
	    case 9: f1block[i] = gp[i].hsmooth * unit_Length / UNIT_L; break;
	    case 10: f1block[i] = aux_gp[i].sfr; break;
	    case 11: f1block[i] = aux_gp[i].delaytime; break;
	    case 12: f1block[i] = aux_gp[i].tmax; break;
	    case 13:
	      for(k=0; k<NMETALS; k++) fZblock[i*NMETALS+k] = 0;
	      fZblock[i*NMETALS+1] = aux_gp[i].metal[0]; // C
	      fZblock[i*NMETALS+3] = aux_gp[i].metal[1]; // O
	      fZblock[i*NMETALS+5] = aux_gp[i].metal[2]; // Si
	      fZblock[i*NMETALS+4] = aux_gp[i].metal[3]; // Fe	      
	      break;
	    default:
	      fprintf(stderr, "Type[%d] has no block *%s*. Quit.\n", type, blkattr[blocknr].name); exit(-1);
	    }
	  }
	} // GAS
	if(type == 1){
	  for (i = 0; i < gheader.npart[type]; i ++){
	    switch (blocknr){
	    case 0:
	      for(k=0; k<3; k++){
		f3block[i*3+k] = (dp[i].pos[k] + 0.5) * unit_Length / UNIT_L;
		if(f3block[i*3+k] > gheader.BoxSize)
		  f3block[i*3+k] -= gheader.BoxSize;
		if(f3block[i*3+k] < 0)
		  f3block[i*3+k] += gheader.BoxSize;
	      }
	      break;
	    case 1:
	      for(k=0; k<3; k++) f3block[i*3+k] = dp[i].vel[k] * unit_Velocity / UNIT_V / gheader.time; break;
	      // V.hdf5 = V.gad * sqrt(a3inv); V.tipsy = V.gad / sqrt(a)
	    case 2: break; // ParticleID
	    case 3: f1block[i] = dp[i].mass * unit_Mass / UNIT_M; break;
	    case 4: f1block[i] = dp[i].phi * (unit_Velocity / UNIT_V) * (unit_Velocity / UNIT_V); break;
	    default:
	      fprintf(stderr, "Type[%d] has no block *%s*. Quit.\n", type, blkattr[blocknr].name); exit(-1);
	    }
	  }
	} // DARK
	if(type == 4){
	  for (i = 0; i < gheader.npart[type]; i ++){
	    switch (blocknr){
	    case 0:
	      for(k=0; k<3; k++) {
		f3block[i*3+k] = (sp[i].pos[k] + 0.5) * unit_Length / UNIT_L;
		if(f3block[i*3+k] > gheader.BoxSize)
		  f3block[i*3+k] -= gheader.BoxSize;
		if(f3block[i*3+k] < 0)
		  f3block[i*3+k] += gheader.BoxSize;
	      }
	      break;
	    case 1:
	      for(k=0; k<3; k++) f3block[i*3+k] = sp[i].vel[k] * unit_Velocity / UNIT_V / gheader.time; break;
	      // V.hdf5 = V.gad * sqrt(a3inv); V.tipsy = V.gad / sqrt(a)
	    case 2: break; // ParticleID
	    case 3: f1block[i] = sp[i].mass * unit_Mass / UNIT_M; break;
	    case 4: f1block[i] = sp[i].phi * (unit_Velocity / UNIT_V) * (unit_Velocity / UNIT_V); break;
	    case 12: f1block[i] = aux_sp[i].tmax; break;
	    case 13: for(k=0; k<NMETALS; k++) fZblock[i*NMETALS+k] = aux_sp[i].metal[k]; break;
	    case 14: f1block[i] = aux_sp[i].age; break;
	    default:
	      fprintf(stderr, "Type[%d] has no block *%s*. Quit.\n", type, blkattr[blocknr].name); exit(-1);
	    }
	  }
	} // STAR
	/* ---------------- */

	fprintf(stdout, "Writing: Type[%d], Block[%d]: %s.\n", type, blocknr, blkattr[blocknr].name);
	hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
	/* Create the first data space, and open it for access */
	/* H5Screate_simple: Create data space. dims[rank] tells the dimension of the data */
	strcpy(hdf5_blkname, blkattr[blocknr].name);
	hdf5_dataset = H5Dcreate(hdf5_grp[type], hdf5_blkname, hdf5_datatype, hdf5_dataspace, H5P_DEFAULT);
	/* Create the dataset named hdf5_blkname, from hdf5_grp, occupying the hdf5_dataspace */
	/* H5Dcreate1 (old version): Last two parameters: space_id and dcpl_id */
	/* H5Dcreate: hdf5_datatype could be H5T_NATIVE_{UINT, DOUBLE, FLOAT64, UINT64, ... } */

	switch (blkattr[blocknr].dtype){
#ifdef LONGIDS	  
	case 1:
	  hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, iblock);
	  /* memspace -> f3block; dataspace -> hdf5_dataspace (in file) */
	  break;
#else	  
	case 1:
	  hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, iblock);	  
	  break;
#endif	  
	case 2:
	  if(dims[1] == 1){
	    hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, f1block);
	    free(f1block);
	  }
	  else if(dims[1] == 3){
	    hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, f3block);
	    free(f3block);
	  }
	  else if(dims[1] == NMETALS){
	    hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, fZblock);
	    free(fZblock);
	  }
	  break;
	case 3:
	  hdf5_status = H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, siblock);
	  free(siblock);
	  break;
	} // switch

	H5Dclose(hdf5_dataset);
	H5Sclose(hdf5_dataspace);
	H5Tclose(hdf5_datatype);
      } // blocknr loop

      hdf5_status = H5Gclose(hdf5_grp[type]);

    } // gheader.npart[type] > 0
  } // type loop
   
  hdf5_status = H5Fclose(hdf5_file);

  return 1;
}

void init_block_info()
{
  int blocknr;
  if((blkattr = (struct block_attributes *) malloc(NBLOCK * sizeof(*blkattr))) == NULL){
    fprintf(stderr, "Failed to allocate memory for blkattr. Quit.\n");
    exit(-1);
  }
  char blkname[NBLOCK][50] =
    {"Coordinates", "Velocities", "ParticleIDs", "Masses", "Potential",
     "InternalEnergy", "Density", "ElectronAbundance", "NeutralHydrogenAbundance", "SmoothingLength",
     "StarFormationRate", "DelayTime", "TemperatureMax", "Metallicity", "StellarFormationTime"
    };
  int dimensions[NBLOCK] =
    {3, 3, 1, 1, 1,
     1, 1, 1, 1, 1,
     1, 1, 1, NMETALS, 1
    };
  int dtypes[NBLOCK] =
    {2, 2, 1, 2, 2,
     2, 2, 2, 2, 2,
     2, 2, 2, 2, 2
    };
  int flags[NBLOCK] =
    {31, 31, 31, 31, 31,
     1, 1, 1, 1, 1,
     1, 1, 17, 17, 16
    };

  for(blocknr = 0; blocknr < NBLOCK; blocknr ++){
    blkattr[blocknr].dim = dimensions[blocknr];
    blkattr[blocknr].dtype = dtypes[blocknr];
    blkattr[blocknr].typeflag = flags[blocknr];
    strcpy(blkattr[blocknr].name, blkname[blocknr]);
  }
}

void write_hdf5_header_attributes(hid_t handle)  
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, gheader.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, gheader.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, gheader.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, gheader.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_doubleprecision);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Potential", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_potential);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_fH2", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_fH2);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_TMax", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_tmax);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_DelayTime", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_delaytime);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
