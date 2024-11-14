# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <sys/types.h>
# include <sys/stat.h>
# include "../sys.h"
# include "../defs.h"
# include "../global.h"
/****************************************************************************/
/* Program analys_conf.c                                                    */
/****************************************************************************/
# define NRUN 48
/****************************************************************************/
int main (int argc,char *argv[])
{
  int i,j,k,run;
  int a1,a2,nruns;
  double o[NOBS],so[NTMP][NOBS];
  double po2[NTMP][NOBS],po[NRUN][NTMP][NOBS];
  double nn1=0,nn2=0,rmsd1,rmsd2;
  long int imd,imd0;
  char fname[100];
  double qcut_a = 58;
  double qcut_b = 76;
  FILE *fp;
  
  if (argc != 4){
    printf("Give: conf filenames in C string style, start (a1), end (a2) indices.\n");
    exit(-1);
  }

  printf("Exec: ");
  for (i = 0; i < argc; i++)
    printf("%s ",argv[i]);
  printf("\n\n");

  FANALYS = 1;
  
  a1 = atoi(argv[2]);
  a2 = atoi(argv[3]);
  nruns = a2 - a1 + 1;
  
  for (i = 0; i < NOBS; i++) {
    for (j = 0; j < NTMP; j++) {
      so[j][i] = po2[j][i] = 0;
      for (k = 0; k < NRUN; k++) po[k][j][i] = 0;
    }
  }

  init(1);
  printinfo();
  Ekin = 0;
  
  printf("<analys_conf> Setting output directory to %s\n",ANADIR);
  mkdir(ANADIR, 0777);
  strcpy(OUTDIR,ANADIR);

  printf("<analys_conf> qcut_a %lf qcut_b %lf \n",qcut_a,qcut_b);

  imd0 = 1;

  for (run = 0; run < nruns; run++) {
    sprintf(fname,argv[1],run + a1);

    imd = 0;
    while (0 == read_conf(0,fname,"")) {
      imd0 += ICONF;
      imd += ICONF;
      
      cart2dof(0);
      Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep=exvol(0))+(Etor=torsion(0))+
	(Econ=cont(0))+(Ehp=hp(0))+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0)); 

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

      //      printf("imd %li %lf\n",imd,Epot);

      if ((imd+1) > NTHERM) {
	so[ind][0]++; for (i = 0 ; i < NOBS; i++) so[ind][i] += o[i];
	po[run][ind][0]++; for (i = 1; i < NOBS; i++) po[run][ind][i] += o[i];
	
	histo_bond(0);
	histo_bend(0);
	histo_tors(0,5);
	histoe(0,Epot);
	histo_cont1(0,ind,nn1);
	histo_cont2(0,ind,nn2);
      }
      
      //      runtime(imd-1,o);
    }

  }
  
  printf("\nRun over\n\n");
  printf("carterr %i therr %i\n\n",carterr,therr);
  printf("\nWriting averages and histograms\n");

  averages(so);
  update_g(so,g);

  for (i = 0; i < NTMP;i++) {
    for (j = 1; j < NOBS; j++) {
      so[i][j] /= so[i][0];
      for (k = 0; k < nruns; k++) {
	po[k][i][j] /= po[k][i][0];
      }
    }
  }

  for (i = 0; i < NTMP; i++) {
    for (j = 1; j < NOBS; j++) {
      for (k = 0; k < nruns; k++) {
	po2[i][j] += (po[k][i][j] - so[i][j]) * (po[k][i][j] - so[i][j]);
      }
    }
  }

  char str[100];
  
  strcpy(str,OUTDIR);
  strcat(str,AVERAGES);

  fp = fopen(str,"w");
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"%i %lf ",i,1./beta[i]);
    for (j = 1; j < NOBS; j++)
      fprintf(fp,"%.5f ",so[i][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  for (k = 1; k < NOBS; k++) {
    sprintf(str,"%s%s_%d",OUTDIR,AVERAGES,k+2);

    fp = fopen(str,"w");    
    for (i = 0; i < NTMP; i++) {
      fprintf(fp,"%lf %i %lf %lf\n",1./beta[i],i,
	      so[i][k],sqrt( po2[i][k]/ nruns / (nruns-1) ) );
    }
    fclose(fp);
  }


  histo_bond(1);
  histo_bend(1);
  histo_tors(2,0);
  histoe(1,0);
  histo_cont1(2,0,0);
  histo_cont2(2,0,0);

  return 0;
}
/****************************************************************************/

