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
  double nn1=0,nn2=0,rmsd1=0,rmsd2=0,rg1=0,rg2=0;
  // double rmsd3=0,rmsd4=0;

  double qcut_a = 95;  // fold A cutoff
  double qcut_b = 180; // fold B cutoff
  
  printf("Exec: ");
  for (i = 0; i < argc; i++)
    printf("%s ",argv[i]);
  printf("\n\n");

  if (NTMP > 1) {
    printf("NTMP > 1\nExiting...\n");
    exit(-1);
  }

  for (i = 0; i < NOBS; i++) {
    for (j = 0; j < NTMP; j++) {
      o[i] = so[j][i] = 0;
    }
  }

  init(1);
  printinfo();
  
  for (imd = imd0; imd < MDSTEP; imd++) {
    mdstep();

    if ((imd+1) % ISAMP == 0) {

      /* Energies (o[] index = column number in output files) */

      o[3]=Ekin; o[4]=Epot; o[5]=Ebon; o[6]=Eben; o[7]=Erep; o[8]=Etor;
      o[9]=Econ1; o[10]=Econ2; o[11]=Ecorr;  o[12]=Ecc; o[13]=Ecb; o[14]=Ehp; o[15]=Eel;

      /* Custom */

      /* two chains */

      rg1 = sqrt( gyr2(iBeg[0],iEnd[0]) );
      rg2 = sqrt( gyr2(iBeg[1],iEnd[1]) );
      rmsd1 = rmsd_calc(xnat,ynat,znat,x,y,z,9,68);
      rmsd2 = rmsd_calc(xnat2,ynat2,znat2,x,y,z,7,53);
      // rmsd3 = rmsd_calc(xnat,ynat,znat,x,y,z,102,161);
      // rmsd4 = rmsd_calc(xnat2,ynat2,znat2,x,y,z,100,146);

      o[16] = rg1;
      o[17] = rg2;
      o[18] = rmsd1;
      o[19] = rmsd2;
      // o[18] = rmsd3;
      // o[19] = rmsd4;
      o[20] = ncont_map1(0);
      o[21] = ncont_map2(0);
      o[22] = ncont_map1(1);
      o[23] = ncont_map2(1);
      o[24] = ncont_map2_inter(0,1);
      o[25] = (ncont_map1(0) >= qcut_a ? 1 : 0);
      o[26] = (ncont_map1(1) >= qcut_a ? 1 : 0);
      o[27] = (ncont_map2(-1) >= qcut_b ? 1 : 0);
      //      o[25] = sqrt(dist2(id1[0],id2[0]));
      //      o[26] = sqrt(dist2(id1[1],id2[1]));


      if ((imd+1) > NTHERM) {
        so[ind][0]++; for (i=0 ; i<NOBS; i++) so[ind][i] += o[i];

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
  printf("<hp> non-native hp pairs calculated: qpair %i\n",qpair);
  printf("<el> non-native el pairs calculated: zpair %i\n",zpair);
  printf("carterr %i therr %i\n\n",carterr,therr);
  dumppdb("_stop.pdb",o,NOBS);

  printf("\nWriting averages and histograms\n");

  averages(so);
  histo_bond(1);
  histo_bend(1);
  histo_tors(2,0);
  histoe(1,0);
  histo_cont1(2,0,0);
  histo_cont2(2,0,0);
  
  return 0;
}
/****************************************************************************/
