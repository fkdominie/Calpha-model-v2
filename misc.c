# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <sys/types.h>
# include <sys/stat.h>
# include "sys.h"
# include "defs.h"
# include "global.h"
/************* miscellaneous ************************************************/
double pi,pi2,pid2;                /* pi,2*pi,pi/2                          */
long seed=-11;                     /* random number seed                    */
long orig_seed;                    /* original seed                         */
char OUTDIR[100];                  /* output directory                      */
int FANALYS = 0;                   /* flag for analysis (post processing)   */ 
/****************************************************************************/
void read_cont_param(char fn[],int ip1[],int ip2[],int npair,double kcont[]);
void read_bonded_param(char fn[],double v1[],double k1[],double v2[],double k2[],
		       double v3[],double k3[], double k4[]);
/****************************************************************************/
/***** INPUT/OUTPUT *********************************************************/
/****************************************************************************/
void printinfo(void) {
  printf("\n");
  printf("Simulation parameters: \n");
  printf("  N %i NCH %i NCR %i \n",N,NCH,NCR);
  printf("  MDSTEP %li NTHERM %i\n",(long)MDSTEP,NTHERM);
  printf("  IRT %i ICONF %i ISTART %i\n",IRT,ICONF,ISTART);
  printf("  NTMP %i TMIN %lf TMAX %lf \n",NTMP, TMIN, TMAX);
  printf("  BOX %i \n",BOX);
  printf("Force field:\n");
  printf("  FF_BOND %i FF_BEND %i FF_TORS %i\n",FF_BOND,FF_BEND,FF_TORS);
  printf("  FF_CONT %i FF_SEQ %i \n",FF_CONT,FF_SEQ);
  printf("  kbon %f kth %f kph1 %f kph3 %f \n",
	 kbon,kth,kph1,kph3);
  printf("  krep %f kcon %f khp %f \n",krep,kcon,khp);
  printf("  eps %f sigsa %f \n",eps,sigsa);
  if (NCR > 0){
    printf("Interaction parameters crowders:\n");
    printf("  epsilonrep %lf rcrowd %lf srefcr %lf\n",epsilonrep,rcrowd,srefcr);
  }
  if (NCH > 0) {
    printf("MD parameters chains:\n");
    printf("  tau %f dt %f gam %f\n",tau,dt,gam);
    printf("  c1 %f c2 %f c3 %f\n",c1,c2,c3);
  }
  if (NCR > 0) {
    printf("MD parameters crowders:\n");
    printf("  Crowder mass  %f \n",mcr);
    printf("  Friction coefficient of crowders %f \n",gamcr);
    printf("  c1cr %f c2cr %f c3cr %f\n",c1cr,c2cr,c3cr);
  }
  printf("----\n\n");
  fflush(0);
}
/****************************************************************************/
void ramachan(char *fn,double b,double th,double ph) {
  FILE *fp;
  char str[100];
  
  strcpy(str,OUTDIR);
  strcat(str,fn);

  fp = fopen(str,"a");
  fprintf(fp,"%.4f %.4f %.4f\n",b,th*180/pi,ph*180/pi);
  fclose(fp);
}
/****************************************************************************/
void runtime(long it,double o[]) {
  int i;
  FILE *fp;
  char str[100];
  
  strcpy(str,OUTDIR);
  strcat(str,RT);

  fp = fopen(str,"a");

  fprintf(fp,"%li %i ",it,ind);
  for (i = 1; i < NOBS; i++)
    fprintf(fp,"%.4lf ",o[i]);
  fprintf(fp,"\n");

 fclose(fp);
} 
/****************************************************************************/
void averages(double so[][NOBS]) {
  int i,j;
  double ntot=0;
  FILE *fp;
  char str[100];
  
  strcpy(str,OUTDIR);
  strcat(str,STATS);

  fp = fopen(str,"w");
  for (i = 0; i < NTMP; i++) ntot += so[i][0];
  for (i = 0; i < NTMP; i++) 
    fprintf(fp,"temp %i %lf %li %lf\n",i, 1/beta[i],
	    (long int) so[i][0], so[i][0] / ntot);
  fprintf(fp,"nflp %li acc %li rate %lf\n",
	  nflp, accflp,
	  (double) accflp / nflp);
  fclose(fp);

  strcpy(str,OUTDIR);
  strcat(str,AVERAGES);
  
  fp = fopen(str,"w");
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"%i %lf ",i,1./beta[i]);
    for (j = 1; j < NOBS; j++)
      fprintf(fp,"%.5f ",so[i][j] / so[i][0]);
    fprintf(fp,"\n");
  }

  fclose(fp);
} 
/****************************************************************************/
void update_g(double so[][NOBS],double g[]) {
  int i;
  double g2[NTMP],ntot = 0;
  char str[100];
  FILE *fp;

  strcpy(str,OUTDIR);
  strcat(str,OUTPUTG);

  fp = fopen(str,"w");
  for (i = 0; i < NTMP; i++) ntot += so[i][0];
  for (i = 0; i < NTMP; i++) 
    g2[i] = g[i] - log( max(so[i][0] / ntot, 0.01 / NTMP ) );
  for (i = 0; i < NTMP; i++)
    fprintf(fp,"%i %lf\n",i, g2[i] - g2[0]);
  fclose(fp);

}
/****************************************************************************/
int read_checkpnt() {
  int i,err = 0;
  long orig;
  double echeck;

  if (read_data("_data",CHECKDIR,&imd0,&seed,&orig)) return 0;

  err += read_conf(0,"_conf",CHECKDIR);

  if (err > 0) {
    printf("<checkpnt> Error reading checkpnt\n");
    printf("<checkpnt> Exiting...\n");
    exit(-1);
    return 0;
  }  
  err += read_momenta("_momenta",CHECKDIR);
  err += read_forces("_forces",CHECKDIR);

  echeck = Epot;
  imd0++;
  orig_seed = orig;
  while (orig < seed) ran3n(&orig);

  cart2dof(0);

  for (i = 0; i < NCR; i++)
    fxc[i] = fyc[i] = fzc[i] = 0;

  for (i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0;

  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep=exvol(0))+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp=hp(0))+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0)); 
  
  Ekin = 0;
  for (i=0;i<N;i++)
    Ekin += vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];

  for (i = 0; i < NCR; i++)
    Ekin += vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i];
  Ekin*=0.5;

  printf("\n************************************");
  printf("\n**** Restarting from checkpoint ****");
  printf("\n************************************");
  printf("\nReading checkpoint directory %s",CHECKDIR);
  printf("\nind %i seed %li orig_seed %li",ind,seed,orig_seed);
  printf("\nEpot: calc %.14e read %.14e diff %.14e",Epot,echeck,Epot-echeck);
  printf("\nEkin: calc %.14e",Ekin);
  printf("\nRestarting from MD cycle imd %li",imd0);
  printf("\n************************************\n\n");
  
  return 1;
}
/****************************************************************************/
void write_checkpnt(void) {
  int i;
  
  cart2dof(0);

  for (i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0;
  for (i = 0; i < NCR; i++)
    fxc[i] = fyc[i] = fzc[i] = 0;

  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep=exvol(0))+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp=hp(0))+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0)); 
  
  Ekin = 0;
  for (i=0;i<N;i++)
    Ekin += vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];

  for (i = 0; i < NCR; i++)
    Ekin += vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i];
  Ekin*=0.5;

  write_data("_data",CHECKDIR,imd,seed,orig_seed);
  write_conf("_conf",CHECKDIR,"w");
  write_momenta("_momenta",CHECKDIR,"w");
  write_forces("_forces",CHECKDIR,"w");

  return ;
} 
/****************************************************************************/
int read_data(char *fname,char *fdir,long *imd0,long *seed,long *orig) {
  FILE *fp;
  char str[100];
  int n,ncr;
  
  strcpy(str,fdir);
  strcat(str,fname);

  if ( (fp = fopen(str,"r")) == NULL) {
    printf("<checkpnt> no checkpoint data found\n");
    return 1;
  }
  printf("<checkpnt> Found file %s \n",str);
  fscanf(fp,"imd %li\n",imd0);
  fscanf(fp,"seed %li\n",seed);
  fscanf(fp,"orig %li\n",orig);
  fscanf(fp,"N %i\n",&n);
  fscanf(fp,"NCR %i\n",&ncr);
  fclose(fp);

  if (n != N || ncr != NCR) {
    printf("<checkpnt> n %i (N %i) ncr %i (NCR %i) \n",n,N,ncr,NCR);
    printf("<checkpnt> unable to restart simulation\n");
    return 1;
  }

  return 0;
} 
/****************************************************************************/
void write_data(char *fname,char *fdir,long imd,long seed,long orig) {
  FILE *fp;
  char str[100];

  strcpy(str,fdir);
  strcat(str,fname);

  fp = fopen(str,"w");
  fprintf(fp,"imd %li\n",imd);
  fprintf(fp,"seed %li\n",seed);
  fprintf(fp,"orig %li\n",orig);
  fprintf(fp,"N %i\n",N);
  fprintf(fp,"NCR %i\n",NCR);
  fclose(fp);
} 
/****************************************************************************/
int read_forces(char *fname,char *fdir) {
  FILE *fp;
  char str[100];

  strcpy(str,fdir);
  strcat(str,fname);

  if ( (fp = fopen(str,"r")) == NULL) return 1;
  fread(frdx,sizeof(double),N,fp);
  fread(frdy,sizeof(double),N,fp);
  fread(frdz,sizeof(double),N,fp);
  fread(frcdx,sizeof(double),NCR,fp);
  fread(frcdy,sizeof(double),NCR,fp);
  fread(frcdz,sizeof(double),NCR,fp);
  fclose(fp);

  return 0;
} 
/****************************************************************************/
void write_forces(char *fn,char *fdir,char *fmode) {
  FILE *fp;
  char str[100];

  strcpy(str,fdir);
  strcat(str,fn);

  fp = fopen(str,fmode);
  fwrite(frdx,sizeof(double),N,fp);
  fwrite(frdy,sizeof(double),N,fp);
  fwrite(frdz,sizeof(double),N,fp);
  fwrite(frcdx,sizeof(double),NCR,fp);
  fwrite(frcdy,sizeof(double),NCR,fp);
  fwrite(frcdz,sizeof(double),NCR,fp);
  fclose(fp);
} 
/****************************************************************************/
int read_momenta(char *fname,char *fdir) {
  FILE *fp;
  char str[100];

  strcpy(str,fdir);
  strcat(str,fname);

  if ( (fp = fopen(str,"r")) == NULL) return 1;
  fread(vx,sizeof(double),N,fp);
  fread(vy,sizeof(double),N,fp);
  fread(vz,sizeof(double),N,fp);
  fread(vxc,sizeof(double),NCR,fp);
  fread(vyc,sizeof(double),NCR,fp);
  fread(vzc,sizeof(double),NCR,fp);
  fclose(fp);

  return 0;
} 
/****************************************************************************/
void write_momenta(char *fn,char *fdir,char *fmode) {
  FILE *fp;
  char str[100];

  strcpy(str,fdir);
  strcat(str,fn);

  fp = fopen(str,fmode);
  fwrite(vx,sizeof(double),N,fp);
  fwrite(vy,sizeof(double),N,fp);
  fwrite(vz,sizeof(double),N,fp);
  fwrite(vxc,sizeof(double),NCR,fp);
  fwrite(vyc,sizeof(double),NCR,fp);
  fwrite(vzc,sizeof(double),NCR,fp);
  fclose(fp);
} 
/****************************************************************************/
void write_conf(char *fn,char *dir,char *fmode) {
  FILE *fp;
  char str[100];

  strcpy(str,dir);
  strcat(str,fn);

  fp = fopen(str,fmode);
  fwrite(&ind,sizeof(int),1,fp);
  fwrite(&Epot,sizeof(double),1,fp);
  fwrite(x,sizeof(double),N,fp);
  fwrite(y,sizeof(double),N,fp);
  fwrite(z,sizeof(double),N,fp);
  fwrite(xcr,sizeof(double),NCR,fp);
  fwrite(ycr,sizeof(double),NCR,fp);
  fwrite(zcr,sizeof(double),NCR,fp);
  fclose(fp);
} 
/****************************************************************************/
int read_conf(int iflag,char *fn,char *fdir) {
  static FILE *fp = NULL;
  static char fcur[200],fname[200];

  if (iflag < 0) {
    strcpy(fcur,"");
    strcpy(fname,"");
    return 0;
  }

  strcpy(fname,fdir);
  strcat(fname,fn);

  if ( strcmp(fcur,fname) ) {
    if (fp != NULL) fclose(fp);

    if ((fp = fopen(fname,"r")) == NULL) {
      printf("<read_conf> Unable to open file %s\n",fname);
      return 1;
    } else 
      printf("<read_conf> Opening file %s\n",fname);
    strcpy(fcur,fname);
  }
  
  if (1 != fread(&ind,sizeof(int),1,fp)) return 1; 
  fread(&Epot,sizeof(double),1,fp);
  fread(x,sizeof(double),N,fp);
  fread(y,sizeof(double),N,fp);
  fread(z,sizeof(double),N,fp);
  fread(xcr,sizeof(double),NCR,fp);
  fread(ycr,sizeof(double),NCR,fp);
  fread(zcr,sizeof(double),NCR,fp);

  return 0;
} 
/****************************************************************************/
void check_cmaps(void) {
  int i,j,m,n;
  double r2;  
  int clash = 0;

  if(FF_CONT < 2) return;

  for (m = 0; m < npair; m++) {
    i = ip1[m]; j = ip2[m]; 
    r2 = ((xnat2[i]-xnat2[j]) * (xnat2[i]-xnat2[j])+
	  (ynat2[i]-ynat2[j]) * (ynat2[i]-ynat2[j])+ 
	  (znat2[i]-znat2[j]) * (znat2[i]-znat2[j]));
    if (nat2[i] && nat2[j] && r2 < distp2[m]) {
      dist_rep1[m] = r2;
      printf("<check_cmaps1> %s contact %i %i (%lf > %lf). Setting repulsive distance to %lf\n",
	     dual1[m] > 0 ? "Dual" : "Native",
	     i,j,distp[m],sqrt(r2),sqrt(dist_rep1[m]));
    }
  }

  for (n = 0; n < npair2; n++) {
    i = ip3[n]; j = ip4[n]; 
    r2 = ((xnat[i]-xnat[j]) * (xnat[i]-xnat[j])+
	  (ynat[i]-ynat[j]) * (ynat[i]-ynat[j])+ 
	  (znat[i]-znat[j]) * (znat[i]-znat[j]));
    if (nat[i] && nat[j] && r2 < distp4[n]) {
      dist_rep2[n] = r2;
      printf("<check_cmaps2> %s contact %i %i (%lf > %lf). Setting repulsive distance to %lf\n",
	     dual2[n] > 0 ? "Dual" : "Native",
	     i,j,distp3[n],sqrt(r2),sqrt(dist_rep2[n]));
    }
  }

  if (clash == 1) {
    printf("Fix clash problem!! Exiting... \n");
    exit(-1);
  }
  
  return;
}
/****************************************************************************/
int get_shared_contacts(int *mc1,int *mc2,int *dual1,int *dual2)
{
  int i,j,m,n,s;

  if (FF_CONT < 2) return 0; 

  for (m = s = 0; m < npair; m++) {
    i = ip1[m]; j = ip2[m];
    for (n = 0; n < npair2; n++) {
      if ( (i == ip3[n] && j == ip4[n]) ||
	   (j == ip3[n] && i == ip4[n]) ) {
	mc1[s] = m;
	mc2[s++] = n;
	dual1[m] = (distp[m] < distp3[n] ? 1 : 2);
	dual2[n] = (distp3[n] < distp[m] ? 1 : 2);
      }
    }
  }

  return s;
}
/****************************************************************************/
void get_natdist(double *dist,double *dist2,double *distg1,double *distg2,
		 double *dist_rep,int n,
		 int *ip1,int *ip2,int *nni1,int *nnj1,int *nni2,int *nnj2,
		 double *xn,double *yn,double *zn) {
  int i,j,i1,j1,m,q,imin=0,jmin=0;
  double r2,r2min = 0,r2max = 0;
  double gmin,dg[4];  

  gmin = dg[0] = dg[1] = dg[2] = dg[3] = 0;

  for (m = 0; m < n; m++) nni1[m] = nnj1[m] = nni2[m] = nnj2[m] = -1;
  
  for (m = 0; m < n; m++) {
    i = ip1[m]; j = ip2[m]; 
    r2 = ((xn[i]-xn[j]) * (xn[i]-xn[j])+
	  (yn[i]-yn[j]) * (yn[i]-yn[j])+ 
	  (zn[i]-zn[j]) * (zn[i]-zn[j]));
    if (r2 < r2min || m == 0) r2min = r2; 
    if (r2 > r2max || m == 0) r2max = r2; 
    dist[m] = sqrt(r2);
    dist2[m] = dist_rep[m] = r2;
    cc[i][j] = cc[j][i] = 1;

    q = 0;
    gmin = 1e6;
    for (i1=i-1; i1<=i+1; i1+=2) {
      for (j1=j-1; j1<=j+1; j1+=2) {
	if (i1 < 0 || i1 > N-1 || j1 < 0 || j1 > N-1) continue;
	dg[q] =  sqrt( (xn[i1]-xn[j1]) * (xn[i1]-xn[j1]) +
		       (yn[i1]-yn[j1]) * (yn[i1]-yn[j1]) + 
		       (zn[i1]-zn[j1]) * (zn[i1]-zn[j1]) );
	if (dg[q] < gmin) {
	  imin = i1; jmin = j1; gmin = dg[q];
	}
	q++;
      }
    }
    
    if (q == 4) {
      if (dg[0] + dg[3] < dg[1] + dg[2]) {
	nni1[m] = i-1; nnj1[m] = j-1; distg1[m] = dg[0];
	nni2[m] = i+1; nnj2[m] = j+1; distg2[m] = dg[3]; 
      } else {
	nni1[m] = i-1; nnj1[m] = j+1; distg1[m] = dg[1]; 
	nni2[m] = i+1; nnj2[m] = j-1; distg2[m] = dg[2]; 
      }
    } else {
      nni1[m] = imin; nnj1[m] = jmin; distg1[m] = gmin;
    }
  }
  
  if (n > 0)
    printf("rmin %f rmax %f\n",sqrt(r2min),sqrt(r2max));
}
/****************************************************************************/
void write_shared_dist(char *fn,int *mc1,int *mc2,int spair) {
  int i,j,l,k,m,n,s;
  char str[100];
  FILE *fp1;

  strcpy(str,TESTDIR);
  strcat(str,fn);

  printf("<write_shared_dist> Found %i shared contacts between %s and %s \n",
	 spair,NATIVE,NATIVE2);    
  printf("<write_shared_dist> Writing to %s\n",str);
  
  fp1 = fopen(str,"w");
  for (s = 0; s < spair; ++s) {
    i = ip1[(m = mc1[s])];
    j = ip2[m];

    k = ip3[(n = mc2[s])];
    l = ip4[n];
    
    if ( !( (i == k && j == l) || (j == k && i == l) ) )
      printf("<write_shared_dist> Error common contact list\n");
    
    fprintf(fp1,"%i %i %lf %lf\n",i,j,distp[m],distp3[n]);
  }
  fclose(fp1);
  return ;
}
/****************************************************************************/
void write_natdist(char *fn,double *dist,int n,int *ip1,int *ip2) {
  int i,j,m;
  char str[100];
  FILE *fp1;

  strcpy(str,TESTDIR);
  strcat(str,fn);
  
  fp1 = fopen(str,"w");

  for (m = 0; m < n; m++) {
    i = ip1[m]; j = ip2[m]; 
    fprintf(fp1,"%i %i %lf\n",i,j,dist[m]);
  }

  fclose(fp1);

  return;
}
/****************************************************************************/
void set_bonded_strength(double *kbond,double *kbend,double *ktor1,double *ktor3,
			 int *nat) {
  int i,j,k,l;

  for (i = 0; i < N-1; i++) {
    j = i + 1;
    kbond[i] = kbon;
  }
  
  for (i = 0; i < N-2; i++) {
    j = i + 1;
    k = i + 2;
    kbend[j] = kth;
  }
  
  for (i = 0; i < N-3; i++) {
    j = i + 1;
    k = i + 2;
    l = i + 3;
    if (nat[i] && nat[j] && nat[k] && nat[l]) {
      ktor1[j] = kph1;
      ktor3[j] = kph3;
    } else {
      ktor1[j] = 0;
      ktor3[j] = 0;
    }
  }
    
  return;
}
/****************************************************************************/
void get_bonded_param(double *bn,double *thn,double *phn,
		      double *xnat,double *ynat,double *znat,
		      int *nat, int *tor,char *fn) {
  int i,j,k,l;
  char str[100];
  FILE *fp1;

  for (i = 0; i < N; i++) {
    x[i] = xnat[i];
    y[i] = ynat[i];
    z[i] = znat[i];
  }

  if (1 != cart2dof(0)) printf("<get_bonded_param> Error native configuration\n");

  for (i = 0; i < N-1; i++) {
    j = i + 1;
    if (nat[i] && nat[j])
      bn[i] =  b[i];
  }
  
  for (i = 0; i < N-2; i++) {
    j = i + 1;
    k = i + 2;
    if (nat[i] && nat[j] && nat[k])
      thn[j] = th[j];
  }
  
  for (i = 0; i < N-3; i++) {
    j = i + 1;
    k = i + 2;
    l = i + 3;
    if (nat[i] && nat[j] && nat[k] && nat[l]) {
      phn[j] = ph[j];
      tor[j] = 1;
    }
  }
  
  strcpy(str,TESTDIR);
  strcat(str,"bonded_param");
  strcat(str,fn);

  fp1 = fopen(str,"w");
  for (i = 0; i < NCH; i++) {
    for (j = iBeg[i]; j <= iEnd[i]; j++) 
      fprintf(fp1,"%3i %3i %10.5f %10.5f %10.5f %3i\n",i,j,
	      bn[j],thn[j]*180./pi,phn[j]*180./pi,nat[j]);
  }
  fclose(fp1);
  
  return;
}
/****************************************************************************/
int relax_chains(int ich) {
  /* ich >= 0 : relax chain ich */
  /* ich <  0 : relax all chains */
  int i,icur,n = 0;
  double rfac = 1.0;
  double Erel,Eold,pho[N];
  double dx,dy,dz;

  if (N == 0) return 0;
  
  for (i = 1; i < N-2; i++) pho[i] = ph[i];
  
  if (ich < 0)
    printf("<init> Relaxing all chains... \n");
  else
    printf("<init> Relaxing chain %i... \n",ich);
    
  Erel = Eold = exvol(0) + cont(0);
 
  while (Erel > N * rfac * eps)  {

    icur = (ich < 0 ? NCH * ran3n(&seed) : ich );
      
    if ( ran3n(&seed) < 0.5 ) { /* turn torsion angle */

      i = iBeg[icur] + ran3n(&seed) * (iEnd[icur] - iBeg[icur] + 1);
      ph[i] =  pi * (2 * ran3n(&seed) - 1);
      dof2cart(0);
	
      if ( (Erel = exvol(0) + cont(0)) < Eold ) {
	  pho[i] = ph[i];
	  Eold = Erel;
	} else {
	  ph[i] = pho[i];
	  dof2cart(0);      
	}
	
      } else {  /* translate chain */

      dx = 10 * (2 * ran3n(&seed) - 1);
      dy = 10 * (2 * ran3n(&seed) - 1);
      dz = 10 * (2 * ran3n(&seed) - 1);
      
      trans(icur, dx, dy, dz);
      
      if ((Erel = exvol(0) + cont(0)) < Eold) {
	  Eold = Erel;
	} else {
	  trans(icur, -dx, -dy, -dz);
      }
    }
    
    if (n % 10000 == 0)  printf("Erel %lf\n",Eold);
    
    n++;    
  } 
  
  printf("<init> ...done in %i steps (Erel %lf) \n",n,Erel);

  return 0;
}
/****************************************************************************/
int relax_crowders(void) {
  int icr,n = 0;
  double rfac = 1.0;
  double Erel,Eold;
  double dx,dy,dz;
  
  if (NCR == 0) return 0;
  
  printf("<init> Relaxing crowders... \n");

  Erel = Eold = crowd_crowd(0) + crowd_bead(0);

  while (Erel > NCR * rfac * eps) {

      icr = NCR * ran3n(&seed);

      dx = 2 * (2 * ran3n(&seed) - 1) ;
      dy = 2 * (2 * ran3n(&seed) - 1) ; 
      dz = 2 * (2 * ran3n(&seed) - 1) ;
      
      trans_cr(icr, dx, dy, dz);

      if ( (Erel = crowd_crowd(0) + crowd_bead(0)) < Eold ) {
	Eold = Erel;
      } else {
	trans_cr(icr, -dx, -dy, -dz);
      }
      
      //           printf("Erel %lf\n",Eold);
      n++;
  } 
  
  printf("<init> ...done in %i steps (Erel %lf) \n",n,Erel);

  return 0;
}
/****************************************************************************/
void read_cont_param(char fn[],int ip1[],int ip2[],int npair,double kcont[]) {
  int i,j,m = 0;
  double kread;
  FILE *fp;
  
  if ( (fp = fopen(fn,"r")) != NULL ) {

    while (3 == fscanf(fp,"%i %i %lf\n",&i,&j,&kread)) {
    
      for (m = 0; m < npair; ++m) {
	if ( (i == ip1[m] && j == ip2[m]) ||
	     (i == ip2[m] && j == ip1[m]) ) break;
      }
      
      if (m == npair) {
	printf("<read_contpar> (%s) unknown contact %i %i %i\n",fn, m,i,j);
	continue;
      }
      
      kcont[m] = kread;
      printf("<read_contpar> (%s) setting strength of contact %i %i %i to %lf\n",
	     fn,m,i,j,kcont[m]);
    }
  }
  
  return ;
}
/****************************************************************************/
void read_bonded_param(char fn[],
		       double v1[],double k1[],
		       double v2[],double k2[],
		       double v3[],double k3[], double k4[]) {
  int j,n = 0;
  char str[100];
  FILE *fp;
  
  if ((fp = fopen(fn,"r")) != NULL) {
    while (1 == fscanf(fp,"%d ",&j)) {

      if (10 == fscanf(fp,"%s %lf %lf %s %lf %lf %s %lf %lf %lf\n",
		      str,v1 + j,k1 + j,
		      str,v2 + j,k2 + j,
		      str,v3 + j,k3 + j,k4 + j) && j >= 0 && j <= N-1) {
	printf("<read_bonded_param> (%s) %i bond %lf %lf bend %lf %lf tor %lf %lf %lf \n",
	       fn,j,
	       v1[j],k1[j],
	       v2[j],k2[j],
	       v3[j],k3[j],k4[j]);
      } else
	printf("<read_param> (%s) read error\n",fn);
      ++n;
    }
  } 

  return ;
}
/****************************************************************************/
/***** INITIALIZATION *******************************************************/
/****************************************************************************/
void init(int iflag) {
  int i,j,k,m,nnat1,nnat2;
  double vxsum,vysum,vzsum,vxcsum,vycsum,vzcsum;
  double c0,c0cr;
  double o[NOBS];
  char c;
  FILE *fp1;

  /* read sequence */

  fp1 = fopen(INPUT,"r");
  k = 0;
  for (j=0; j<NCH; j++) {
    iBeg[j] = k;
    while ((c=getc(fp1)) != '\n') {
      if (c >= 'A' && c <= 'z') {
	seq[k] = c;
	a2c[k] = j;
	++k;
      }
    }
    iEnd[j] = k-1;
  }
  fclose(fp1);

  /* read g parameters */

  if (IREADG == 1 && NTMP > 1) {
    fp1 = fopen(INPUTG,"r");
    if (fp1 == NULL) {
      printf("<init> No file %s\n",INPUTG);
      for (i=0; i<NTMP; i++) g[i] = 0.0;
    } else {
      for (i=0; i<NTMP; i++) fscanf(fp1,"%i %lf", &i, &g[i]);
      fclose(fp1);
    }
  }

  /* output directory */

  strcpy(OUTDIR,RESDIR);
  
  /* set seed */ 

  if (ISEED == 1){
    FILE *devrandom = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(seed), 1, devrandom);
    seed=-labs(seed%100000000);
    fclose(devrandom);
  } 

  orig_seed = seed;
  printf("orig_seed = %li\n",seed);

  beta[0] = 1. / TMIN;
  for (i = 0; i < NTMP; i++)    
    beta[i] = 1./TMAX * pow(TMAX/TMIN,(double)i/max(NTMP-1,1));

  if (N > 0) {
    printf("Chains: \n");
    for (j=0; j<NCH; j++) {
      printf("%3d %3d %3d ",j,iBeg[j],iEnd[j]);
      for (i=iBeg[j]; i<=iEnd[j]; i++) printf("%c",seq[i]); 
      printf("\n");
    }
    printf("\n");
  }
  
  printf("Index  temp  beta  g:\n");
  for (i=0; i<NTMP; i++) {
    printf("%i %lf %lf %lf\n",i,1./beta[i], beta[i], g[i]);
  }
  printf("\n");

  /* create results directory */

  printf("<init> Creating directory %s\n",RESDIR);
  mkdir(RESDIR, 0777);
  printf("<init> Creating directory %s\n",TESTDIR);
  mkdir(TESTDIR, 0777);
  printf("<init> Creating directory %s\n",CHECKDIR);
  mkdir(CHECKDIR, 0777);


  /* constants */
  
  pi=acos(-1.);
  pi2=2.*pi;
  pid2=pi/2.;
  cthmin=cos(pi/90.);
  sthmin=sin(pi/90.);

  /* energy parameters */

  kbon*=eps;
  kth*=eps;
  kph1*=eps;
  kph3*=eps;
  kcon*=eps;
  khp*=eps;
  krep*=eps;
  epsilonrep*=eps;
  eclash*=eps;
  
  kcon60 = kcon * 60;
  krep12 = krep * 12;
  kbon2 = kbon * 2;
  kth2 = kth * 2;

  /* set temperatures */
  
  ind = NTMP * ran3n(&seed);
  printf("<init> Setting temperature index to %i\n",ind);

  /* MD parameters */

  /*************** Beads  **************/
  dt=0.005*tau;
  gam=0.05/tau;
  mbd = 1.0;
  c0=gam*dt/ (2);
  c1=dt*dt/ (2.*mbd);
  c2=(1-c0)*(1-c0+c0*c0);
  c3=dt*(1-c0+c0*c0) / (2.*mbd);

  /************ Crowders **************/
  mcr = 4 * (rcrowd/8.) * (rcrowd/8.);
  gamcr = (0.025/tau) * (8./rcrowd);
  c0cr = (gamcr*dt) / (2);
  c1cr = (dt*dt) / (2.*mcr);
  c2cr = (1-c0cr)*(1-c0cr+c0cr*c0cr);
  c3cr = (dt * (1-c0cr+c0cr*c0cr)) / (2.*mcr);

  for (i=0; i<NTMP; ++i) {
    tconstbd[i] = sqrt(2*mbd*gam/dt/beta[i]);
    tconstcr[i] = sqrt(2*mcr*gamcr/dt/beta[i]);
  }

  /* Packing fractions */

  if (NCR > 0) {
    double rcrowd_effective = rcrowd + 0.5;
    double phi_cr = NCR * (4./3) * pi * pow(rcrowd_effective,3.0) / pow(BOX,3.0);
    printf("<init> Volume fraction occupied\n");
    printf("  Crowders: %lf \n",phi_cr);
  }
  
  if (iflag == 0) return ;

  /* native structures */

  if ( (nnat1 = read_native(NATIVE,xnat,ynat,znat,nat)) > 0) 
    printf("<init> NATIVE: Read %i residue positions (%s)\n",nnat1,NATIVE);
  else
    printf("<init> NATIVE: No data (%s)\n",NATIVE);

  if ( (nnat2 = read_native(NATIVE2,xnat2,ynat2,znat2,nat2)) > 0) 
    printf("<init> NATIVE2: Read %i residue positions (%s)\n",nnat2,NATIVE2);
  else
    printf("<init> NATIVE2: No data (%s)\n",NATIVE2);

  /* contact maps */
  double dtmp[MAXP];
  int ntmp[MAXP];
  
  npair = npair2 = spair = 0;
  npair3 = npair4 = npair5 = 0;
  for (m = 0; m < MAXP; m++)
    ip1[m] = ip2[m] = ip3[m] = ip4[m] = ip5[m] = ip6[m] = ip7[m] = ip8[m] = ip9[m] = ip10[m] = 0;
  for (i=0; i < N; i++) for (j = 0;j < N; j++) cc[i][j] = 0;  
  for (m=0; m < MAXP; m++) dual1[m] = dual2[m] = 0;    

  npair = read_contacts(CONTMAP,ip1,ip2);
  if (npair > 0) {
    printf("<init> CONTMAP: Read %i contacts (%s)\n",npair, CONTMAP);
    printf("<init> CONTMAP: ");
    get_natdist(distp,distp2,distg1,distg2,dist_rep1,npair,ip1,ip2,
		nni1,nnj1,nni2,nnj2,xnat,ynat,znat);
    printf("<init> CONTMAP: Writing to %s\n","cmap.out");
    write_natdist("cmap.out",distp,npair,ip1,ip2);
  } else printf("<init> CONTMAP: No data (%s)\n",CONTMAP);

  npair2 = read_contacts(CONTMAP2,ip3,ip4);
  if (npair2 > 0) {
    printf("<init> CONTMAP2: Read %i contacts  (%s)\n",npair2, CONTMAP2);
    printf("<init> CONTMAP2: ");
    get_natdist(distp3,distp4,distg3,distg4,dist_rep2,npair2,ip3,ip4,
		nni3,nnj3,nni4,nnj4,xnat2,ynat2,znat2);
    printf("<init> CONTMAP2: Writing to %s\n","cmap2.out");
    write_natdist("cmap2.out",distp3,npair2,ip3,ip4);
  } else printf("<init> CONTMAP2: No data (%s)\n",CONTMAP2);

  npair3 = read_contacts(CONTMAP3,ip5,ip6);
  if (npair3 > 0) {
    printf("<init> CONTMAP3: Read %i contacts (%s)\n",npair3, CONTMAP3);
    printf("<init> CONTMAP3: ");
    get_natdist(distp5,distp6,dtmp,dtmp,dtmp,npair3,ip5,ip6,
		ntmp,ntmp,ntmp,ntmp,xnat,ynat,znat);
    printf("<init> CONTMAP3: Writing to %s\n","cmap3.out");
    write_natdist("cmap3.out",distp5,npair3,ip5,ip6);
  } else printf("<init> CONTMAP3:  No data (%s)\n",CONTMAP3);

  npair4 = read_contacts(CONTMAP4,ip7,ip8);
  if (npair4 > 0) {
    printf("<init> CONTMAP4: Read %i contacts  (%s)\n",npair4, CONTMAP4);
    printf("<init> CONTMAP4: ");
    get_natdist(distp7,distp8,dtmp,dtmp,dtmp,npair4,ip7,ip8,
		ntmp,ntmp,ntmp,ntmp,xnat,ynat,znat);
    printf("<init> CONTMAP4: Writing to %s\n","cmap4.out");
    write_natdist("cmap4.out",distp7,npair4,ip7,ip8);
  } else printf("<init> CONTMAP4: No data (%s)\n",CONTMAP4);

  if (npair  > MAXP) {printf("npair too big\n"); exit(-1);}
  if (npair2 > MAXP) {printf("npair2 too big\n"); exit(-1);}
  if (npair3 > MAXP) {printf("npair3 too big\n"); exit(-1);}
  if (npair4 > MAXP) {printf("npair4 too big\n"); exit(-1);}

  for (m = 0; m < npair; m++) kcon_nat1[m] = kcon;
  for (m = 0; m < npair2; m++) kcon_nat2[m] = kcon;

  if (FF_CONT > 1) {
    spair = get_shared_contacts(mc1,mc2,dual1,dual2);
    write_shared_dist("cmap_corr.out",mc1,mc2,spair);
    check_cmaps();
  }

  /* contact parameters from file */  

  read_cont_param(CONTPAR,ip1,ip2,npair,kcon_nat1);
  read_cont_param(CONTPAR2,ip3,ip4,npair2,kcon_nat2);

  /* bonded interactions */

  /* default */   

  for (i = 0; i < N-1; i++) bn[i] = bn2[i] = 3.8;
  for (i = 1; i < N-1; i++) thn[i] = thn2[i] = 120.0 * pi / 180.0; 
  for (i = 1; i < N-2; i++) phn[i] = phn2[i] = 120.0 * pi / 180.0;

  /* reference values from native structure */
  
  if (nnat1 > 0) get_bonded_param(bn,thn,phn,xnat,ynat,znat,nat,tor,"_nat.out");
  if (nnat2 > 0) get_bonded_param(bn2,thn2,phn2,xnat2,ynat2,znat2,nat2,tor2,"_nat2.out");

  /* default strengths */
  
  set_bonded_strength(kbond,kbend,ktor1,ktor3,nat);
  set_bonded_strength(kbond2,kbend2,ktor1_2,ktor3_2,nat2);

  /* bonded parameters from file */  

  read_bonded_param(BONDEDPAR,bn,kbond,thn,kbend,phn,ktor1,ktor3);
  read_bonded_param(BONDEDPAR2,bn2,kbond2,thn2,kbend2,phn2,ktor1_2,ktor3_2);

  /* initialize functions */
  
  dof2cart(-1);
  crowd_crowd(-1);
  crowd_bead(-1);
  bond(-1);
  bend(-1);
  torsion(-1);
  exvol(-1);
  cont(-1);
  hp(-1);

  read_conf(-1,"","");
  histo_bond(-1);
  histo_bend(-1);
  histo_tors(-1,0);
  histoe(-1,0);
  histo_cont1(-1,0,0);
  histo_cont2(-1,0,0);

  printf("<init> npair %d npair2 %d npair3 %d npair4 %d npair5 %d \n",npair,npair2,npair3,npair4,npair5);

  /* initial chain configuration  */  

  if (ISTART == 0) {            /* native */
    printf("<init> Initializing chain(s) from NATIVE %s\n",NATIVE);
    for (i = 0; i < N; i++) {
      x[i] = xnat[i];
      y[i] = ynat[i];
      z[i] = znat[i];
      //      printf("%lf %lf %lf \n",x[i],y[i],z[i]);
    }
    if (1 != cart2dof(0)) printf("Error initial configuration");
  } else if (ISTART == 1) {     /* read */
    printf("<init> Initializing chain(s) from START %s\n",START);
    fp1 = fopen(START,"r");
    for (i = 0; i < N; i++) {
      if (4 != fscanf(fp1,"%i %lf %lf %lf",&j,&x[i],&y[i],&z[i]))
	break;
    }
    fclose(fp1);
    if (1 != cart2dof(0)) printf("Error initial configuration");
  } else if (ISTART == 2) {     /* random */
    printf("<init> Initializing chain(s) from random configuration\n");
    for (i = 0; i < N-1; i++) b[i]  =  bn[i] ;
    for (i = 1; i < N-1; i++) th[i] =  thn[i] ;
    for (i = 1; i < N-2; i++) ph[i] =  pi * (2 * ran3n(&seed) - 1);
    dof2cart(0);
    relax_chains(-1);
  } else {
    printf("\nInvalid ISTART\n");
    exit(-1);
  }

  //  printf(" Rgyr: 1-56 %lf 8-52 %lf \n",sqrt(gyr2(0,N-1)),sqrt(gyr2(7,51)));
  //  exit(-1);
  
  /* Initial conformation of crowders (random) */
  
  for (i = 0; i < NCR; i++){
    xcr[i] = BOX * ran3n(&seed);
    ycr[i] = BOX * ran3n(&seed);
    zcr[i] = BOX * ran3n(&seed);
  }
  
  relax_crowders();

  vxsum=vysum=vzsum=0;
  for (i=0;i<N;i++) {
    //    printf("%i bl %f th %f ph %f\n",i,bn[i],180*thn[i]/pi,180*phn[i]/pi);
    vx[i]=sqrt(1/beta[ind]/mbd)*gasdev2();
    vy[i]=sqrt(1/beta[ind]/mbd)*gasdev2();
    vz[i]=sqrt(1/beta[ind]/mbd)*gasdev2();
    vxsum+=vx[i];
    vysum+=vy[i];
    vzsum+=vz[i];
  }

  vxcsum=vycsum=vzcsum=0;
  for (i=0;i<NCR;i++) {
    vxc[i]=sqrt(1/beta[ind]/mcr)*gasdev2();
    vyc[i]=sqrt(1/beta[ind]/mcr)*gasdev2();
    vzc[i]=sqrt(1/beta[ind]/mcr)*gasdev2();
    vxcsum+=vxc[i];
    vycsum+=vyc[i];
    vzcsum+=vzc[i];
  }

  Ekin=0;
  for (i=0;i<N;i++) {
    vx[i]-=vxsum / max(N,1);
    vy[i]-=vysum / max(N,1);
    vz[i]-=vzsum / max(N,1);
    Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  for (i=0;i<NCR;i++) {
    vxc[i]-=vxcsum / max(NCR,1);
    vyc[i]-=vycsum / max(NCR,1);
    vzc[i]-=vzcsum / max(NCR,1);
    Ekin+=vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i];
  }
  Ekin*=0.5;

  
  /* initial energies and forces */
  
  for (i = 0; i<N; i++) fx[i]=fy[i]=fz[i]=0;
  for (i = 0; i < NCR; i++) fxc[i]=fyc[i]=fzc[i]=0;
    Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep=exvol(0))+(Etor=torsion(0))+
      (Econ=cont(0))+(Ehp=hp(0))+(Ecc = crowd_crowd(0))+(Ecb = crowd_bead(0));

  for (i=0;i<NCR;i++) {
    frcdx[i] = gasdev2() * tconstcr[ind];
    frcdy[i] = gasdev2() * tconstcr[ind];
    frcdz[i] = gasdev2() * tconstcr[ind];
  }
  
  for (i=0;i<N;i++) {
    frdx[i] = gasdev2() * tconstbd[ind];
    frdy[i] = gasdev2() * tconstbd[ind];
    frdz[i] = gasdev2() * tconstbd[ind];
  }

  if (!FANALYS && read_checkpnt()) return ;

  printf("<init> saving initial conformation to %s\n","_start.pdb");
  dumppdb("_start.pdb",o,0);

  printf("\n");
  printf("Initial conformation (%d chains, %d crowders):\n",NCH,NCR);
  printf("  Ekin %f Epot %f\n",Ekin,Epot);
  printf("  Ebon %f Eben %f Erep %f Etor %f Econ %f Econ1 %lf Econ2 %lf Ecorr %lf\n",
	 Ebon,Eben,Erep,Etor,Econ,Econ1,Econ2,Ecorr);
  printf("  Ehp %f Ecc %f Ecb %f\n",Ehp, Ecc, Ecb);

  
  /* print energy functions */
  
  bond(1);
  bend(1);
  torsion(1);
  cont_corr(1);
  cont(1);
  cont2(1);
  crowd_crowd(1);
  crowd_bead(1);  
  
  return ;
}
/****************************************************************************/
