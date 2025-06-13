/****************************************************************************/
/* pdb2cont.c                                                               */
/* link with readpdb3.c                                                     */
/****************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "../sys.h"
# include "../defs.h"
# include "../global.h"
/****************************************************************************/
int disreg[N];
/****************************************************************************/
int main(int argc,char *argv[]){
  int i,j,k;
  FILE *fp;

  if (argc != 2){
    printf("Input: Give file with two columns specifying disordered regions\n");
    printf("       Column 1: bead index (from 0 to N-1) \n)");
    printf("       Column 2: 1 = bead disordered, 0 = not disordered. \n");
    printf("Output: Apply output file to BONDEDPARAM in defs.h to\n");
    printf("        load generated parameter values\n");
    exit(-1);
  }

  init(1);
  printinfo();

  if ( (fp = fopen(argv[1],"r")) != NULL) {
    while (2 == fscanf(fp,"%d %d",&j,&k) && feof(fp) == 0) disreg[j] = k;
    fclose(fp);
  }
  
  for (i = 0; i < N-1; ++i) {
    if (disreg[i] == 1) {
      bn[i] = 3.8;
      thn[i] = thn_dis * deg2rad;
      phn[i][0] = phn_dis[0] * deg2rad;
      phn[i][1] = phn_dis[1] * deg2rad;
      phn[i][2] = phn_dis[2] * deg2rad;
      
      kbond[i] = kbon;
      kbend[i] = kth;
      ktor[i][0] = kph_dis[0]; ktor[i][1] = kph_dis[1]; ktor[i][2] = kph_dis[2];
    }
  }
  
  printf("Writing parameters to file: bonded_out\n");
  write_bonded_param(bn,thn,phn,kbond,kbend,ktor,"./","bonded_out");
  
  return 0;
}
/****************************************************************************/
