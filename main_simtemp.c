# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "sys.h"
# include "defs.h"
# include "global.h"
/****************************************************************************/
/* Langevin dynamics simulation of coarse-grained polymer chain             */
/****************************************************************************/
int main (int argc,char *argv[])
{
  int i,j;
  double o[NOBS],so[NTMP][NOBS];
  double nn1=0,nn2=0,rmsd1=0.0,rmsd2=0.0;
  double qcut_a = 58;
  double qcut_b = 76;
    
  printf("Exec: ");
  for (i = 0; i < argc; i++)
    printf("%s ",argv[i]);
  printf("\n\n");

  printf("<main_simtemp> qcut_a %lf qcut_b %lf \n",qcut_a,qcut_b);

  for (i = 0; i < NOBS; i++) {
    for (j = 0; j < NTMP; j++) {
      o[i] = so[j][i] = 0;
    }
  }

  init(1);
  printinfo();

  for (imd = imd0; imd < MDSTEP; imd++) {
    
    if ((imd+1) % IFLIP == 0)
      tflip(Epot+Ekin);

    mdstep();
    
    if ((imd+1) % ISAMP == 0) {
      rmsd1 = rmsd_calc(xnat,ynat,znat,x,y,z,0,N-1);
      rmsd2 = rmsd_calc(xnat2,ynat2,znat2,x,y,z,0,N-1);
      
      o[1]=Ekin; o[2]=Epot; o[3]=Ebon; o[4]=Eben; o[5]=Erep; o[6]=Etor;
      o[7]=Econ1; o[8]=Econ2; o[9]=Ecorr;  o[10]=Ecc; o[11]=Ecb;

      o[12] = rmsd1; 
      o[13] = rmsd2;
      o[14] = (nn1 = no_cont()) / max(npair,1);  
      o[15] = (nn2 = no_cont2()) / max(npair2,1);
      o[16] = (nn1 > qcut_a ? 1 : 0);
      o[17] = (nn2 > qcut_b ? 1 : 0);

      if ((imd+1) > NTHERM) {
        so[ind][0]++; for (i=1 ; i<NOBS; i++) so[ind][i] += o[i];

        histo_bond(0);
        histo_bend(0);
        histo_tors(0,5);
        histoe(0,Epot);
	histo_cont1(0,ind,nn1);
	histo_cont2(0,ind,nn2);
      }
    }
     
    if ((imd+1) % IRT == 0) {
      runtime(imd,o);
    }

    if ((imd+1) % ICONF == 0) {
      write_conf(CONF,OUTDIR,"a");
    }
    
    if ((imd+1) % ICHECK == 0){
      averages(so);
      update_g(so,g);
      histo_bond(1);
      histo_bend(1);
      histo_tors(1,0);
      histoe(1,0);
      histo_cont1(1,0,0);
      histo_cont2(1,0,0);

      dumppdb(PDB,o,NOBS);
      write_checkpnt();
    }
  }

  printf("\nRun over\n\n");
  printf("carterr %i therr %i\n\n",carterr,therr);
  dumppdb("_stop.pdb",o,0);

  printf("\nWriting averages and histograms\n");

  averages(so);
  update_g(so,g);
  histo_bond(1);
  histo_bend(1);
  histo_tors(2,0);
  histoe(1,0);
  histo_cont1(2,0,0);
  histo_cont2(2,0,0);

  return 0;
}
/****************************************************************************/
