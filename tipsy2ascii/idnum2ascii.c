#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void main(int argc, char **argv)
{
  int idnum;
  unsigned char *buf;
  char basename[200], inputname[200], outputname[200];
  FILE *inputfile, *outputfile;
  int CORRUPTED_FILE=0;
  strcpy(basename, argv[1]);
  strcat(strcpy(inputname, basename), ".idnum");
  strcat(strcpy(outputname, basename), ".id");

  if(!(inputfile = fopen(inputname, "r")))
    {printf("Can't open input file! Force exit.\n"); exit(-1);}
  outputfile = fopen(outputname, "w");
  
  printf("Reading begins.\n");
  printf("sizeof(int) = %i\n", sizeof(int));
    while(!feof(inputfile))
      {	fread(&idnum, sizeof(int), 1, inputfile);
	if(!feof(inputfile))
	  fprintf(outputfile, "%d\n", idnum);
	if(CORRUPTED_FILE == 0)
	  {if(idnum <= 0) CORRUPTED_FILE = 1;}
      }
    // ID is always one line more!!
  fclose(inputfile);
  fclose(outputfile);
  if (CORRUPTED_FILE == 1) printf("CORRUPTED ID FILE!!!\n");
  else printf("Good File.\n");
  printf("done.\n");
}
