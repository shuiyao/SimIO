#include "gadgetdefs.h"
#include "extern.h"

/* this routine loads particle data from Gadget's default
 * binary file format. A snapshot may be distributed into multiple files.
 */
int load_snapshot(char *buf, int type)
{

  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph,pc_star;
  FILE   *fd;

  int allocate_memory();

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

      fd = fopen(buf,"r");
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      NumPart= header1.npart[type];
      pc = 1;

      for(k=0, ntot_withmasses=0; k<5; k++) if(header1.mass[k]==0) ntot_withmasses+= header1.npart[k];

      allocate_memory();
      SKIP; // position
      for(k=0,pc_new=pc;k<6;k++) {
	  if(k==type) {
	      for(n=0;n<header1.npart[k];n++) {
		  fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
		  pc_new++;
		}
	    }
	  else
	    fseeko(fd, sizeof(float)*3*header1.npart[k], SEEK_CUR);
	}
      SKIP;

      SKIP; //velocity
      for(k=0,pc_new=pc;k<6;k++)
	{
	  if(k==type)
	    {
	      for(n=0;n<header1.npart[k];n++)
		{
		  fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
		  pc_new++;
		}
	    }
	  else
	    fseeko(fd, sizeof(float)*3*header1.npart[k], SEEK_CUR);
	}
      SKIP;
    

      SKIP; // id
      for(k=0,pc_new=pc;k<6;k++)
	{
	  if(k==type)
	    {
	      for(n=0;n<header1.npart[k];n++)
		{
		  fread(&Id[pc_new], sizeof(int), 1, fd);
		  pc_new++;
		}
	    }
	  else
	    fseeko(fd, sizeof(int)*header1.npart[k], SEEK_CUR);
	}
      SKIP;

      
      /* doing masses stuff -------------------------------------*/

      if(ntot_withmasses>0) SKIP;

	for(k=0, pc_new=pc; k<6; k++) {
	  if (k==type) {
		for(n=0;n<header1.npart[k];n++) {
		  if(header1.mass[k]==0) fread(&P[pc_new].Mass, sizeof(float), 1, fd);
		  else P[pc_new].Mass= header1.mass[k];
                  //if(pc_new%1000==0 || pc_new==1 || pc_new==header1.npart[k]) printf("ALL MASS %10d % 5.3e\n",pc_new,P[pc_new].Mass);
#ifdef OUTPUT_DEBUG
		  if (n==1)
		    printf("Mass=%g\n", P[pc_new].Mass);
#endif

	          pc_new++;
	        }
	  }
	  else if(header1.mass[k]==0) 
		fseeko(fd,header1.npart[k]*sizeof(float),SEEK_CUR);
	  }

      if(ntot_withmasses>0) SKIP;

      /*----------------------------------------------------------*/

      /* Potential Block */
#ifndef POT_BLOCK_BACKWARD
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          if(k==type)
            {
              for(n=0;n<header1.npart[k];n++)
                {
                  fread(&P[pc_new].Phi, sizeof(float), 1, fd);
                  //if(pc_new%1000==0 || pc_new==1 || pc_new==header1.npart[k]) printf("ALL POT %10d % 5.3e\n",pc_new,P[pc_new].Phi);
                  pc_new++;
                }
            }
          else
            fseeko(fd, sizeof(float)*header1.npart[k], SEEK_CUR);
        }
      SKIP;
#endif
      

      if(header1.npart[0]>0 && type==0) // Gas only
	{
	  SKIP; // Internal Energy
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Temp, sizeof(float), 1, fd);
#ifdef OUTPUT_DEBUG
	      if (n==1)
		printf("Temp=%g\n", P[pc_sph].Temp);
#endif
	      pc_sph++;
	    }
	  SKIP;

	  SKIP; // Density
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP; // Ne
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
#ifdef OUTPUT_DEBUG
		  if (n==1)
		    printf("Ne=%g\n", P[pc_sph].Ne);
#endif
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	

      /* now read neutral hydrogen array */

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Nh, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Nh= 1.0;
		pc_sph++;
	      }
	


      /* now read hsm array */

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Hsml, sizeof(float), 1, fd);
#ifdef OUTPUT_DEBUG
		  if (n==1)
		    printf("Hsml=%g\n", P[pc_sph].Hsml);
#endif
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Hsml= 1.0;
		pc_sph++;
	      }
	
      /* now read star formation rate array */

	  if(header1.flag_sfr)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Sfr, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Sfr= 1.0;
		pc_sph++;
	      }

	  // SHUIYAO: DEPENDING ON WHETHER -DOUTPUTDELAYTIME IS TURNED ON!!!
#ifdef OUTPUTDELAYTIME
	  printf("Reading the delayed time.\n");
	  SKIP; // Delayed Time
          for(n=0, pc_sph=pc; n<header1.npart[0];n++) 
	    {
	      fread(&P[pc_sph].DelayTime, sizeof(float), 1, fd);
	      //if(pc_sph%1000==0  || pc_sph==1 || pc_sph==header1.npart[4]) printf("SPH DELAYTIME %10d %5.3e\n",pc_sph,P[pc_sph].DelayTime);

	      pc_sph++;
	    }
          SKIP;
#else
	  printf("No delaytime output.\n");
          for(n=0, pc_sph=pc; n<header1.npart[0];n++) 
	    {
	      P[pc_sph].DelayTime = 0.;
	      pc_sph++;
	    }
#endif

#ifdef AGE_BLOCK_FORWARD
          /* skip over age block, if it exists */
	  printf("Age Block is moved forward.\n");
          if(header1.npart[4]>0) {
            SKIP;
            fseeko(fd,header1.npart[4]*sizeof(float),SEEK_CUR);
            SKIP;
	  }
#endif

	  if(header1.flag_sfr) {
            SKIP;
            for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
              fread(&P[pc_sph].metal[0], sizeof(float), 1, fd);
              fread(&P[pc_sph].metal[1], sizeof(float), 1, fd);
              fread(&P[pc_sph].metal[2], sizeof(float), 1, fd);
              fread(&P[pc_sph].metal[3], sizeof(float), 1, fd);
	      //if(pc_sph%1000==0 || pc_sph==1 || pc_sph==header1.npart[0]) printf("SPH METAL %10d %5.3e %5.3e %5.3e %5.3e\n",pc_sph,P[pc_sph].metal[0],P[pc_sph].metal[1],P[pc_sph].metal[2],P[pc_sph].metal[3]);
              pc_sph++;
            }
	    fseeko(fd,header1.npart[4]*4*sizeof(float),SEEK_CUR);
            SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
              P[pc_sph].metal[0]= 0.0;
              P[pc_sph].metal[1]= 0.0;
              P[pc_sph].metal[2]= 0.0;
              P[pc_sph].metal[3]= 0.0;
              pc_sph++;
            }

          SKIP; // Tmax
          for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
            fread(&P[pc_sph].Tmax, sizeof(float), 1, fd);
#ifdef OUTPUT_DEBUG
	    if (n==1)
	      printf("Tmax=%g\n", P[pc_sph].Tmax);
#endif
            pc_sph++;
          }
	  fseeko(fd,header1.npart[4]*sizeof(float),SEEK_CUR);
          SKIP;

	  SKIP; // Nspawn
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++){
            fread(&P[pc_sph].Nspawn, sizeof(int), 1, fd);
#ifdef OUTPUT_DEBUG
	    if (n==1)
	      printf("Nspawn=%d\n", P[pc_sph].Nspawn);
#endif
	    pc_sph++;
	  }
	  fseeko(fd,header1.npart[4]*sizeof(int),SEEK_CUR);
          SKIP;

#ifndef AGE_BLOCK_FORWARD
          /* skip over age block, if it exists */
          if(header1.npart[4]>0) {
            SKIP;
            fseeko(fd,header1.npart[4]*sizeof(float),SEEK_CUR);
            SKIP;
	  }
#endif
	}

      if(header1.npart[4]>0 && type==4) //ONLY FOR STAR PARTICLES
	{
	  /* skip over u,rho,Ne,Nh,Hsml,Sfr,DelayTime blocks */
	  for(i=0; i<7; i++ ) {
	    SKIP;
	    fseeko(fd,header1.npart[0]*sizeof(float),SEEK_CUR);
	    SKIP;
	  }

#ifdef AGE_BLOCK_FORWARD
	  SKIP;
	  for(n=0, pc_star=pc; n<header1.npart[4];n++)
	    {
	      fread(&P[pc_star].age, sizeof(float), 1, fd);
	      pc_star++;
	    }
	  SKIP;
#endif

	  /* now read age and metallicity for stars into arrays */
	  SKIP;
	  fseeko(fd,header1.npart[0]*4*sizeof(float),SEEK_CUR);
	  for(n=0, pc_star=pc; n<header1.npart[4];n++)
	    {
              fread(&P[pc_star].metal[0], sizeof(float), 1, fd);
              fread(&P[pc_star].metal[1], sizeof(float), 1, fd);
              fread(&P[pc_star].metal[2], sizeof(float), 1, fd);
              fread(&P[pc_star].metal[3], sizeof(float), 1, fd);
	      //if(pc_star%1000==0 || pc_star==1 || pc_star==header1.npart[0]) printf("STAR METAL %10d %5.3e %5.3e %5.3e %5.3e\n",pc_star,P[pc_star].metal[0],P[pc_star].metal[1],P[pc_star].metal[2],P[pc_star].metal[3]);
              pc_star++;
	    }
	  SKIP;

          SKIP;
	  fseeko(fd,header1.npart[0]*sizeof(float),SEEK_CUR);
	  for(n=0, pc_star=pc; n<header1.npart[4];n++){
            fread(&P[pc_star].Tmax, sizeof(float), 1, fd);
	    //if(pc_star%1000==0  || pc_star==1 || pc_star==header1.npart[4]) printf("STAR TMAX %10d %5.3f\n",pc_star,P[pc_star].Tmax);
	    pc_star++;
	  }
	  SKIP;

	  SKIP;
	  fseeko(fd,header1.npart[0]*sizeof(int),SEEK_CUR);
	  for(n=0, pc_star=pc; n<header1.npart[4];n++){
	    fread(&P[pc_star].Nspawn, sizeof(int), 1, fd);
	    //if(pc_star%1000==0 || pc_star==1 || pc_star==header1.npart[4]) printf("STAR NSPAWN %10d %2d\n",pc_star,P[pc_star].Nspawn);
	    pc_star++;
          }
	  SKIP;

#ifndef AGE_BLOCK_FORWARD
	  SKIP;
	  for(n=0, pc_star=pc; n<header1.npart[4];n++)
	    {
	      fread(&P[pc_star].age, sizeof(float), 1, fd);
	      //if(pc_star%1000==0  || pc_star==1 || pc_star==header1.npart[4]) printf("STAR AGE %10d %5.3e\n",pc_star,P[pc_star].age);
	      pc_star++;
	    }
	  SKIP;
#endif
	  		  
	}

#ifdef POT_BLOCK_BACKWARD
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          if(k==type)
            {
              for(n=0;n<header1.npart[k];n++)
                {
                  fread(&P[pc_new].Phi, sizeof(float), 1, fd);
                  //if(pc_new%1000==0 || pc_new==1 || pc_new==header1.npart[k]) printf("ALL POT %10d % 5.3e\n",pc_new,P[pc_new].Phi);
                  pc_new++;
                }
            }
          else
            fseeko(fd, sizeof(float)*header1.npart[k], SEEK_CUR);
        }
      SKIP;
#endif

      fclose(fd);
      return(0);
}

