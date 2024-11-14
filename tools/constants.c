/*******************************************************************/
/* calculates constants                                            */
/*******************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
/*******************************************************************/
int seq[5000];             /* current amino acid sequence           */
int iBeg[5000],iEnd[5000];  /* first, last aa in chain               */
/*******************************************************************/
int main(int argc,char *argv[])
{
  int j,k;
  int c,NCH=1,NCR=0;
  FILE *fp1;
  
  if (argc == 1) {
    printf("Usage: Give input filename, #chains (NCH) and #crowders (NCR).\n");
    printf("Output: File sys.h with constants N, NCH, NCR, etc\n");
    exit(-1);
  }

  printf("argc %i\n",argc);
  if (argc >= 3) NCH = atoi(argv[2]);
  if (argc == 4) NCR = atoi(argv[3]);
 
  /* sequence */

  fp1 = fopen(argv[1],"r");
  k=0;
  for (j=0; j<NCH; j++) {
    iBeg[j] = k;
    while ((c=getc(fp1)) != '\n') {
      if (c >= 'A' && c <= 'z') {
	seq[k] = c;
	++k;
      }
    }
    iEnd[j] = k-1;    
  }
  fclose(fp1);

  printf("/************* geometry ************/\n");
  printf("# define N %i        /* # monomers */\n",k);
  printf("# define NCH %i      /* # chains   */\n",NCH);
  printf("# define NCR %i      /* # crowders   */\n",NCR);
  printf("/***********************************/\n");  

  fp1=fopen("sys.h","w");
  fprintf(fp1,"/************* geometry ************/\n");
  fprintf(fp1,"# define N %i        /* # monomers */\n",k);
  fprintf(fp1,"# define NCH %i      /* # chains   */\n",NCH);
  fprintf(fp1,"# define NCR %i      /* # crowders   */\n",NCR);
  fprintf(fp1,"/***********************************/\n");  
  fclose(fp1);

  return 0;
}
/**************************************************************/ 


