/* Last Update on 2014-03-19 */
/* Add a factor of 1.4 in the .eps field */

/* Update on 2014-02-20 */
/* In this version, data are dumped to one single core (ThisTask = 0) that writes. */
/* In the old version, zero values are generated in random snapshots on UMASS, no clue yet. */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef OUTPUT_TIPSY

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"
#include "tipsydefs_n2.h"



/*! \file io_tipsy.c
 *  \brief Output of a snapshot file to disk in tipsy binary/aux/idnum format.
 */

static int n_type[6];
static long long ntot_type_all[6];
static MyFloat unit_Time, unit_Density, unit_Length, unit_Mass, unit_Velocity;
static FILE *binfile, *auxfile, *idfile;
#ifdef OUTPUT_TIPSY_AW
static FILE *awfile;
#endif
static int task;
static int n_for_this_task;

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savetipsy(int num)
{
  size_t bytes;
  int i,n;
  void tipsy_write_gas(),tipsy_write_dark(),tipsy_write_star(),tipsy_write_aux_gas(),tipsy_write_aux_star(),tipsy_write_idnum();
#ifdef OUTPUT_TIPSY_AW
  void tipsy_write_aw_gas();
#endif
  double t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

#if defined(SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();
  /* ensures that new tree will be constructed */
  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);
#endif

  if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(2);
    }

      /* determine global and local particle numbers */
      for(n = 0; n < 6; n++)
        n_type[n] = 0;

      for(n = 0; n < NumPart; n++)
        n_type[P[n].Type]++;

      sumup_large_ints(6, n_type, ntot_type_all);

      if(DumpFlag) {
        if(ThisTask == 0){
	  printf("\nwriting tipsy snapshot files... \n");
	  /* Open files */

	  char buf[500];
	  sprintf(buf, "%s%s_%03d.bin", All.OutputDir, All.SnapshotFileBase, num);
	  binfile = fopen(buf, "w");
	  sprintf(buf, "%s%s_%03d.aux", All.OutputDir, All.SnapshotFileBase, num);
	  auxfile = fopen(buf, "w");
	  sprintf(buf, "%s%s_%03d.idnum", All.OutputDir, All.SnapshotFileBase, num);
	  idfile = fopen(buf, "w");
#ifdef OUTPUT_TIPSY_AW
	  sprintf(buf, "%s%s_%03d.aw", All.OutputDir, All.SnapshotFileBase, num);
	  awfile = fopen(buf, "w");
#endif
	}
	
/* Output tipsy header */
      if(ThisTask == 0) {
        tipsy_header.time = All.Time; 
        tipsy_header.ndim = 3;
	for( i=0, tipsy_header.nbodies=0; i<6; i++ ) tipsy_header.nbodies += ntot_type_all[i];
	tipsy_header.nsph = ntot_type_all[0];
	tipsy_header.ndark = ntot_type_all[1];
	tipsy_header.nstar = ntot_type_all[4];
//	fprintf(stdout,"savetipsy: %g %d %d %d %d\n",tipsy_header.time,tipsy_header.nbodies,tipsy_header.nsph,tipsy_header.ndark,tipsy_header.nstar);
	my_fwrite(&tipsy_header, sizeof(tipsy_header), 1, binfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);

/* Output tipsy binary */
      t0 = second();
      cosmounits();

      tipsy_write_gas(0);
      tipsy_write_dark(1);
      tipsy_write_star(4);

      MPI_Barrier(MPI_COMM_WORLD);
      if(ThisTask == 0)
	fclose(binfile);
      t1 = second();
      if(ThisTask == 0) printf("outputting tipsy binary file took = %g sec\n", timediff(t0, t1));

/* Output auxiliary file */
      t0 = second();
      tipsy_write_aux_gas(0);
      tipsy_write_aux_star(4);
      
      MPI_Barrier(MPI_COMM_WORLD);
      if(ThisTask == 0)
	fclose(auxfile);
      t1 = second();
      if(ThisTask == 0) printf("outputting auxiliary file took = %g sec\n", timediff(t0, t1));

#ifdef OUTPUT_TIPSY_AW
/* Output analytic wind file */
      t0 = second();
      tipsy_write_aw_gas(0);
      
      MPI_Barrier(MPI_COMM_WORLD);
      if(ThisTask == 0)
	fclose(awfile);
      t1 = second();
      if(ThisTask == 0) printf("outputting analytic wind file took = %g sec\n", timediff(t0, t1));
#endif

/* Output idnum file */
      t0 = second();
      tipsy_write_idnum();

      MPI_Barrier(MPI_COMM_WORLD);
      if(ThisTask == 0)
	fclose(idfile);
      t1 = second();
      if(ThisTask == 0) printf("outputting idnum file took = %g sec\n", timediff(t0, t1));
  }
      myfree(CommBuffer);
      return;
}

/* Outputs tipsy gas particles, with Gadget particle type itype (usually =0) */
void tipsy_write_gas(int itype)
{
  struct gas_particle gasp;
  MyOutputFloat Temp,metal_tot;
  MyFloat MeanWeight,a3inv=1;
  int i,j;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;

  struct gas_particle *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(gasp);

/* Output gas particles */
  if(All.ComovingIntegrationOn) a3inv = 1. / (All.Time*All.Time*All.Time);

  for(task = 0; task < NTask; task++){
    // IMPORT/EXPORT N_TYPE
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT, 
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT, 
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);
    //By Now, On ALL cores, n_for_this_task = task.n_type[itype]

    while(n_for_this_task > 0){
      write_buffer = (struct gas_particle *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

    if(task == ThisTask){
      // FILL BUFFER
      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	gasp.mass = P[i].Mass*All.UnitMass_in_g/unit_Mass;
        gasp.pos[0] = P[i].Pos[0]*All.UnitLength_in_cm/unit_Length-0.5;
        gasp.pos[1] = P[i].Pos[1]*All.UnitLength_in_cm/unit_Length-0.5;
        gasp.pos[2] = P[i].Pos[2]*All.UnitLength_in_cm/unit_Length-0.5;
        gasp.vel[0] = P[i].Vel[0]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        gasp.vel[1] = P[i].Vel[1]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        gasp.vel[2] = P[i].Vel[2]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        gasp.phi = P[i].p.Potential*All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s/(unit_Velocity*unit_Velocity);

	MeanWeight= 4.0/(3*HYDROGEN_MASSFRAC+1+4*HYDROGEN_MASSFRAC*SphP[i].Ne)*PROTONMASS;

#ifndef DENSITY_INDEPENDENT_SPH
	Temp = DMAX(All.MinEgySpec, SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));	// internal energy so far
#else
	Temp = DMAX(All.MinEgySpec, SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].EgyWtDensity * a3inv, GAMMA_MINUS1));
#endif
#if defined(ANALYTIC_WINDS) && !defined(WIND_ENTROPY)
	if(SphP[i].wind_flag == 1)
	  Temp = DMAX(All.MinEgySpec, SphP[i].Entropy);
#endif	
	Temp *= MeanWeight/BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs/ All.UnitMass_in_g;		// now convert to T
        gasp.temp = Temp;
        gasp.hsmooth = SphP[i].Hsml*All.UnitLength_in_cm/unit_Length;
	/* gasp.hsmooth = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density) * All.cf_afac3; SoundSpeed */
        gasp.rho = SphP[i].d.Density*All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*unit_Density);
#ifdef ANALYTIC_WINDS
	if(SphP[i].wind_flag == 1)
	  /* gasp.rho = SphP[i].ambient_density*All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*unit_Density); */
	  gasp.rho = get_cloud_density(i)*All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*unit_Density);
#endif	
#ifdef METALS
	for( j=0,metal_tot=0.; j<NMETALS; j++ ) metal_tot += P[i].Metallicity[j];
        gasp.metals = metal_tot;
#else
	gasp.metals = 0.;
#endif
	*write_buffer++ = gasp;

	pc --;
	if( pc == 0 ) break;
      }
      offset = i+1; // tested, i++ is not done
	// Send data ONLY to write_task ONLY WHEN ThisTask != 0
	if(ThisTask != 0)
	  MPI_Ssend(CommBuffer, sizeof(gasp)*blocklen, MPI_BYTE, 
		    0, TAG_PDATA, MPI_COMM_WORLD);
	// Ssend: blocking synchronous send
    } // task != ThisTask
    else{
      // Receive data ONLY when task != (ThisTask == 0)
      if(ThisTask == 0)
	MPI_Recv(CommBuffer, sizeof(gasp)*blocklen, MPI_BYTE,
		 task, TAG_PDATA, MPI_COMM_WORLD, &status);
    }
    if(ThisTask == 0)
      my_fwrite(CommBuffer, sizeof(gasp), blocklen, binfile);
    // NOTE: n_for_this_task must be imported from "task" to writeTask(0)
    n_for_this_task -= blocklen;
    } // while loop
  } // task loop
  return;
}

/* Outputs tipsy dark particles, with Gadget particle type itype (ex: 1,2,3) */
void tipsy_write_dark(int itype)
{
  struct dark_particle darkp;
  int i;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;

  struct dark_particle *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(darkp);

/* Output halo particles */
  for(task = 0; task < NTask; task++){
    // IMPORT/EXPORT N_TYPE
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT, 
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT, 
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

    while(n_for_this_task > 0){
      write_buffer = (struct dark_particle *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

    if(task == ThisTask){
      // FILL BUFFER

      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	darkp.mass = P[i].Mass*All.UnitMass_in_g/unit_Mass;
        darkp.pos[0] = P[i].Pos[0]*All.UnitLength_in_cm/unit_Length-0.5;
        darkp.pos[1] = P[i].Pos[1]*All.UnitLength_in_cm/unit_Length-0.5;
        darkp.pos[2] = P[i].Pos[2]*All.UnitLength_in_cm/unit_Length-0.5;
        darkp.vel[0] = P[i].Vel[0]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        darkp.vel[1] = P[i].Vel[1]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        darkp.vel[2] = P[i].Vel[2]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
        darkp.phi = P[i].p.Potential*All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s/(unit_Velocity*unit_Velocity);
        darkp.eps = 1.4 * All.SofteningHalo*All.UnitLength_in_cm/unit_Length;
	*write_buffer++ = darkp;
	pc --;
	if( pc == 0 ) break;
      }
      offset = i+1;
	// Send data
	if(ThisTask != 0)
	  MPI_Ssend(CommBuffer, sizeof(darkp)*blocklen, MPI_BYTE, 
		    0, TAG_PDATA, MPI_COMM_WORLD);
    }
    else{
      if(ThisTask == 0)
	MPI_Recv(CommBuffer, sizeof(darkp)*blocklen, MPI_BYTE,
		 task, TAG_PDATA, MPI_COMM_WORLD, &status);
    }
    if(ThisTask == 0)
      my_fwrite(CommBuffer, sizeof(darkp), blocklen, binfile);
    n_for_this_task -= blocklen;
    }
  }
  return;
}

/* Outputs tipsy star particles, with Gadget particle type itype (usually =4) */
void tipsy_write_star(int itype)
{
  struct star_particle starp;
  MyOutputFloat metal_tot;
  double get_time_from_step();
  int i,j;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;

  struct star_particle *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(starp);

  /* Output star particles */

  for(task = 0; task < NTask; task++){
    // IMPORT/EXPORT N_TYPE
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT, 
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT, 
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

    while(n_for_this_task > 0){
      write_buffer = (struct star_particle *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

    if(task == ThisTask){
      // FILL BUFFER
      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	starp.mass = P[i].Mass*All.UnitMass_in_g/unit_Mass;
	starp.pos[0] = P[i].Pos[0]*All.UnitLength_in_cm/unit_Length-0.5;
	starp.pos[1] = P[i].Pos[1]*All.UnitLength_in_cm/unit_Length-0.5;
	starp.pos[2] = P[i].Pos[2]*All.UnitLength_in_cm/unit_Length-0.5;
	starp.vel[0] = P[i].Vel[0]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
	starp.vel[1] = P[i].Vel[1]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
	starp.vel[2] = P[i].Vel[2]*All.UnitVelocity_in_cm_per_s/unit_Velocity/sqrt(All.Time);
	starp.phi = P[i].p.Potential*All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s/(unit_Velocity*unit_Velocity);
	starp.eps = 1.4 * All.SofteningHalo*All.UnitLength_in_cm/unit_Length;
#ifdef STELLARAGE
	starp.tform = get_time_from_step(All.Ti_Current);
#endif
#ifdef METALS
	for( j=0,metal_tot=0.; j<NMETALS; j++ ) metal_tot += P[i].Metallicity[j];
	starp.metals = metal_tot;
#else
	starp.metals = 0.;
#endif
	*write_buffer++ = starp; // Tricky 1
	pc --;
	if( pc == 0 ) break;
      }
	if(ThisTask != 0)
	  MPI_Ssend(CommBuffer, sizeof(starp)*blocklen, MPI_BYTE, 
		    0, TAG_PDATA, MPI_COMM_WORLD);
    } 
    else{
      if(ThisTask == 0)
	MPI_Recv(CommBuffer, sizeof(starp)*blocklen, MPI_BYTE,
		 task, TAG_PDATA, MPI_COMM_WORLD, &status);
    }
    if(ThisTask == 0)
      my_fwrite(CommBuffer, sizeof(starp), blocklen, binfile);
    n_for_this_task -= blocklen;
    }
    // NOTE: n_for_this_task must be imported from "task" to writeTask(0)
  } // task loop
  return;
}

void tipsy_write_aux_gas(int itype)
{
  struct aux_gas_data auxgp;
  MyFloat ne,nh0,nHeII;
  MyFloat a3inv=1;
  int i,j;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;

  if(All.ComovingIntegrationOn) a3inv = 1 / (All.Time*All.Time*All.Time);

  struct aux_gas_data *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(auxgp);

  for(task = 0; task < NTask; task++){
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT,
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT,
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

    while(n_for_this_task > 0){
      write_buffer = (struct aux_gas_data *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

      if(task == ThisTask){
      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	for( j=0; j<NMETALS; j++ ) auxgp.metal[j] = P[i].Metallicity[j];
	/* auxgp.metal[0] = SphP[i].h_eff; */
	/* auxgp.metal[1] = SphP[i].MaxSignalVel * All.cf_afac3; */
	/* auxgp.metal[2] = SphP[i].NV_DivVel; */
	/* auxgp.metal[3] = SphP[i].NV_dt_DivVel; */
#ifdef SFR
	auxgp.sfr = get_starformation_rate(i);
#else
	auxgp.sfr = 0;
#endif
#ifdef OUTPUTTMAX
	auxgp.tmax = P[i].Tmax;
#else
	auxgp.tmax = 0.0;
#endif
#ifdef WINDS
	auxgp.delaytime = SphP[i].DelayTime;
#else
	auxgp.delaytime = 0.;
#endif
	ne = SphP[i].Ne;
	AbundanceRatios(DMAX(All.MinEgySpec, SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv , GAMMA_MINUS1)), SphP[i].d.Density * a3inv, &ne, &nh0, &nHeII);
	auxgp.ne = ne;
	auxgp.nh = nh0;
#ifdef SFR
	auxgp.nspawn = P[i].NumStarsGen;
#ifdef ANALYTIC_WINDS_DEBUG
	auxgp.nspawn = P[i].TimeBin;
#endif	
#else
	auxgp.nspawn = 0;
#endif
/*       SHUIYAO: 13-04-09 */
#ifdef WINDS
	auxgp.nrec = P[i].Nrec;
#ifdef ANALYTIC_WINDS_DEBUG
	auxgp.nrec = SphP[i].numngb_wind / SphP[i].numngb_evaluated;
#endif	
#else
	auxgp.nrec = 0;
#endif
      /* SHUIYAO: 13-07-16 */
#ifdef OUTPUT_ALPHA
	auxgp.alpha = SphP[i].alpha;
#endif
	*write_buffer++ = auxgp;
	pc --;
	if(pc == 0) break;
      }
      offset = i + 1;
      if(ThisTask != 0)
	MPI_Ssend(CommBuffer, sizeof(auxgp)*blocklen, MPI_BYTE,
		  0, TAG_PDATA, MPI_COMM_WORLD);
      }// task == ThisTask
      else{
	if(ThisTask == 0)
	  MPI_Recv(CommBuffer, sizeof(auxgp)*blocklen, MPI_BYTE,
		   task, TAG_PDATA, MPI_COMM_WORLD, &status);
      }
      if(ThisTask == 0)
	my_fwrite(CommBuffer, sizeof(auxgp), blocklen, auxfile);
      n_for_this_task -= blocklen;
  } // while loop
}// task loop
  return;
}

void tipsy_write_aux_star(int itype)
{
  struct aux_star_data auxsp;
  int i,j;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;

  struct aux_star_data *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(auxsp);
  
  for(task = 0; task < NTask; task++){
    // IMPORT/EXPORT N_TYPE
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT, 
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT, 
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

    while(n_for_this_task > 0){
      write_buffer = (struct aux_star_data *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

    if(task == ThisTask){
      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	for( j=0; j<NMETALS; j++ ) auxsp.metal[j] = P[i].Metallicity[j];
#ifdef STELLARAGE
	auxsp.age = P[i].StellarAge;	// formation expansion factor (All.Time)
#else
	auxsp.age = 0;
#endif
#ifdef OUTPUTTMAX
	auxsp.tmax = P[i].Tmax;
#else
	auxsp.tmax = 0;
#endif
#ifdef SFR
	auxsp.nspawn = P[i].NumStarsGen;
#else
	auxsp.nspawn = 0;
#endif
/*       SHUIYAO: 13-04-09 */
#ifdef WINDS
	auxsp.nrec = P[i].Nrec;
#else
	auxsp.nrec = 0;
#endif
	*write_buffer++ = auxsp;
	pc --;
	if( pc == 0 ) break;
      }
      offset = i+1;
      	if(ThisTask != 0)
	  MPI_Ssend(CommBuffer, sizeof(auxsp)*blocklen, MPI_BYTE, 
		    0, TAG_PDATA, MPI_COMM_WORLD);
    } 
    else{
      if(ThisTask == 0)
	MPI_Recv(CommBuffer, sizeof(auxsp)*blocklen, MPI_BYTE,
		 task, TAG_PDATA, MPI_COMM_WORLD, &status);
    }
    if(ThisTask == 0)
      my_fwrite(CommBuffer, sizeof(auxsp), blocklen, auxfile);
    n_for_this_task -= blocklen;
    }
  } // task loop
  return;
}

void tipsy_write_idnum()
{
  int i;
  MPI_Status status;
  int pc, blockmaxlen, blocklen;
  int itype;

  MyIDType *write_buffer;
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(MyIDType);  

  // itype here
  for(itype = 0; itype < 6; itype++){
    if(itype == 2) continue;
    if(itype == 3) continue;
    if(itype == 5) continue;
    for(task = 0; task < NTask; task++){
    // IMPORT/EXPORT N_TYPE
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++) // Starting from CPU 0 to Ntask-1
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT, 
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT, 
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);
    // Every CPU receive n_type[itype] from CPU i 

    while(n_for_this_task > 0){
      write_buffer = (MyIDType *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;
  
      if( task == ThisTask){
	for( i=0; i<NumPart; i++ ) {
	  if( P[i].Type != itype ) continue;
	  *write_buffer++ = P[i].ID;
	  pc --;
	  if(pc == 0) break;
      }
	// Send data
	if(ThisTask != 0)
	  MPI_Ssend(CommBuffer, sizeof(MyIDType)*blocklen, MPI_BYTE, 
		    0, TAG_PDATA, MPI_COMM_WORLD);
    }
    else{
      if(ThisTask == 0)
	MPI_Recv(CommBuffer, sizeof(MyIDType)*blocklen, MPI_BYTE,
		 task, TAG_PDATA, MPI_COMM_WORLD, &status);
    }
    if(ThisTask == 0)
      my_fwrite(CommBuffer, sizeof(MyIDType), blocklen, idfile);
    n_for_this_task -= blocklen;
	
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

      return;
}
	
/* conversion to tipsy "standard" units: L=1,Omega=1 */
void cosmounits(void)
{
    unit_Time=sqrt(8*M_PI/3)*CM_PER_MPC/(100*All.HubbleParam*1.e5);
    unit_Density=1.8791E-29*All.HubbleParam*All.HubbleParam;
    unit_Length=All.BoxSize*CM_PER_MPC*1.e-3;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length/(All.HubbleParam*All.HubbleParam);
    unit_Velocity=unit_Length/unit_Time;
    return;
}

/* ************ MPI Parallel I/O Stuff *********************** */

/* Opens snapshot tipsy output file with given suffix for snapshot num */
MPI_File My_MPI_open(char *suffix, int num)
{
    MPI_File fh;
    char mpi_err_str[MPI_MAX_ERROR_STRING];
    int  mpi_err_strlen;
    int  mpi_err;
    char buf[500];

    sprintf(buf, "%s%s_%03d.%s", All.OutputDir, All.SnapshotFileBase, num, suffix);
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, buf,
            MPI_MODE_WRONLY | MPI_MODE_CREATE,
            MPI_INFO_NULL, &fh)) != MPI_SUCCESS) {
        MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
        printf("Proc %d: ", ThisTask);
        printf("MPI_File_open failed (%s)\n", mpi_err_str);
        endrun(1013);
    }
    return fh;
}

#endif // OUTPUT_TIPSY

#ifdef OUTPUTCOOLINGRATE
void write_coolingrate(int num)
{
  int i;
  char fname[300], mkcmd[200];
  FILE *fcool;
  double dtcool, dtcond, dmcool, u, dmdtcool;
  sprintf(mkcmd, "mkdir -p %sCOOLING%03d", All.OutputDir, num);
  system(mkcmd);
  sprintf(fname, "%sCOOLING%03d/%s.%d", All.OutputDir, num, "cooling", ThisTask);
  if(!(fcool = fopen(fname, "w")))
	  {
	    printf("can't open file `%s`\n", fname);
	    endrun(1183);
	  }
  for(i=0;i<N_gas;i++){
    u = SphP[i].Entropy / GAMMA_MINUS1 * pow(get_Density_for_Energy(i), GAMMA_MINUS1);
    dtcool = u / SphP[i].dUcool * All.UnitTime_in_s / 3.1536e16; // Gyr
    dtcond = u / SphP[i].dUcond * All.UnitTime_in_s / 3.1536e16; // Gyr
    // dtcond = inf if dUcond == 0.0
    dmdtcool = P[i].Mass / All.HubbleParam / dtcool; // 10^10Msolar/Gyr
    fprintf(fcool, "%d %g %g  %g  %g %g %g\n",
	    P[i].ID, SphP[i].Mgal, dmdtcool,
	    dtcool,
	    SphP[i].dUcool, SphP[i].dUcond, SphP[i].dUvisc
	    );
    }
    fclose(fcool);
}
#endif // OUTPUTCOOLINGRATE

#ifdef OUTPUT_TIPSY_AW
void tipsy_write_aw_gas(int itype)
{
  struct aw_gas_data awgp;
  MyFloat a3inv=1;
  int i,j;
  MPI_Status status;
  int pc, blockmaxlen, blocklen, offset=0;
  double Uint, Temp, MeanWeight, u_ambient, dtime;

  cosmounits();
  if(All.ComovingIntegrationOn) a3inv = 1 / (All.Time*All.Time*All.Time);

  struct aw_gas_data *write_buffer; // aw
  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(awgp);

  for(task = 0; task < NTask; task++){
    if(task == ThisTask){
      n_for_this_task = n_type[itype];
      for(i=0; i<NTask; i++)
	if(i != ThisTask)
	  MPI_Send(&n_for_this_task, 1, MPI_INT,
		   i, TAG_NFORTHISTASK, MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n_for_this_task, 1, MPI_INT,
	       task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

    while(n_for_this_task > 0){
      write_buffer = (struct aw_gas_data *) CommBuffer;
      pc = n_for_this_task;
      if(pc > blockmaxlen)
	pc = blockmaxlen;
      blocklen = pc;

      if(task == ThisTask){
      for( i=offset; i<NumPart; i++ ) {
	if( P[i].Type != itype ) continue;
	/* Set Value for awgp here */
	/* -------------------------------- */
	awgp.wind_flag = SphP[i].wind_flag;
#ifndef DISABLE_WIND_MASSLOSS
	if(awgp.wind_flag > 0) awgp.wind_mass = P[i].Mass;
	else awgp.wind_mass = SphP[i].WindMass;
#else
	awgp.wind_mass = -1.0;
#endif
	MeanWeight= 4.0/(3*HYDROGEN_MASSFRAC+1+4*HYDROGEN_MASSFRAC*SphP[i].Ne)*PROTONMASS;
	// Internal Energy and Cooling Time
	Uint = SphP[i].Entropy / GAMMA_MINUS1 * pow(get_Density_for_Energy(i), GAMMA_MINUS1);
	if(SphP[i].wind_flag == 1) Uint = SphP[i].Entropy;
	awgp.dtcool = Uint / SphP[i].dUcool * All.UnitTime_in_s / 3.1536e16; // Gyr
	// dU/dt from wind impact; or dU/dt from wind
	awgp.dudt = SphP[i].dUcond * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs; // physical
	dtime = ( (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval /
		   All.cf_hubble_a * All.UnitTime_in_s) / All.HubbleParam;
	awgp.dudt /= dtime;
	// Ambient Temperature and Density
	Temp = DMAX(All.MinEgySpec, Uint);
	Temp *= MeanWeight/BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs/ All.UnitMass_in_g;
	if(SphP[i].wind_flag == 1)
	  u_ambient = SphP[i].u_ambient / SphP[i].kernel_sum_mass;
	  Temp = GAMMA_MINUS1 / BOLTZMANN * (u_ambient) * PROTONMASS * 0.62;
	awgp.temp = Temp;
        awgp.rho = SphP[i].d.Density*All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*unit_Density);
	if(SphP[i].wind_flag == 1)
	  awgp.rho = SphP[i].ambient_density*All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm*All.UnitLength_in_cm*unit_Density);
	// Cloud Radius
	if(SphP[i].wind_flag > 0)
	  awgp.rcloud = SphP[i].r_cloud / CM_PER_MPC * 1.e3;
	else awgp.rcloud = SphP[i].Hsml * All.atime / All.HubbleParam;
	/* -------------------------------- */
	*write_buffer++ = awgp;
	pc --;
	if(pc == 0) break;
      }
      offset = i + 1;
      if(ThisTask != 0)
	MPI_Ssend(CommBuffer, sizeof(awgp)*blocklen, MPI_BYTE,
		  0, TAG_PDATA, MPI_COMM_WORLD);
      }// task == ThisTask
      else{
	if(ThisTask == 0)
	  MPI_Recv(CommBuffer, sizeof(awgp)*blocklen, MPI_BYTE,
		   task, TAG_PDATA, MPI_COMM_WORLD, &status);
      }
      if(ThisTask == 0)
	my_fwrite(CommBuffer, sizeof(awgp), blocklen, awfile);
      n_for_this_task -= blocklen;
  } // while loop
}// task loop
  return;
}
#endif
