# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <sys/types.h>
# include <sys/stat.h>
/*********************/
# include "sys.h"
# include "defs.h"
# include "global.h"
/************* miscellaneous ************************************************/
double pi,pi2,pid2;                /* pi,2*pi,pi/2                          */
double deg2rad,rad2deg;            /* conversion                            */
long seed=-11;                     /* random number seed                    */
long orig_seed;                    /* original seed                         */
char OUTDIR[100];                  /* output directory                      */
int FANALYS = 0;                   /* flag for analysis (post processing)   */
FILE *fp_log;
/****************************************************************************/
void read_cont_param(char fn[],int ip1[],int ip2[],int npair,double kcont[]);
void set_bonded_param(double *bn,double *thn,double phn[],
		      double *xnat,double *ynat,double *znat,int nnat);
void write_bonded_param(double *bn,double *thn,double *phn,char *fn);
void cont_param_salt(char fn[],int ip1[],int ip2[],int n,double kcont[]);
void read_disregs(char fn[],int *dis);
/****************************************************************************/
/***** INPUT/OUTPUT *********************************************************/
/****************************************************************************/
void printinfo(void) {
  printf("\n");
  printf("Simulation parameters: \n");
  printf("  N %i NCH %i NCR %i NPROT %i \n",N,NCH,NCR,NPROT);
  printf("  MDSTEP %li NTHERM %i\n",(long)MDSTEP,NTHERM);
  printf("  IRT %i ICONF %i ISTART %i\n",IRT,ICONF,ISTART);
  printf("  NTMP %i TMIN %lf TMAX %lf \n",NTMP, TMIN, TMAX);
  printf("  BOX %i \n",BOX);
  printf("Force field:\n");
  printf("  FF_BOND %i FF_BEND %i FF_TORS %i FF_CONT %i\n",
	 FF_BOND,FF_BEND,FF_TORS,FF_CONT);
  printf("  FF_SEQ %i FF_EL %i FF_SALT %i FF_DISULF %i\n",
	 FF_SEQ,FF_EL,FF_SALT,FF_DISULF);
  if (NCR > 0){
    printf("Interaction parameters crowders:\n");
    printf("  rcrowd %lf srefcr %lf\n",rcrowd,srefcr);
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
  fprintf(fp,"%.4f %.4f %.4f\n",b,th*rad2deg,ph*rad2deg);
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
  for (i = 3; i < NOBS; i++)
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
    for (j = 3; j < NOBS; j++)
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
int read_checkpnt(void) {
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

  cart2dof();

  for (i = 0; i < NCR; i++)
    fxc[i] = fyc[i] = fzc[i] = 0;

  for (i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0;

  exvol(0,&Erep,&Ehp,&Eel); 
 
  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep)+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp)+(Eel)+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0));

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
  
  cart2dof();

  for (i = 0; i < N; i++)
    fx[i] = fy[i] = fz[i] = 0;
  for (i = 0; i < NCR; i++)
    fxc[i] = fyc[i] = fzc[i] = 0;

  exvol(0,&Erep,&Ehp,&Eel); 
 
  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep)+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp)+(Eel)+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0));

  
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
  if (1 == fscanf(fp,"imd %li\n",imd0)) printf("<checkpnt> imd %li\n",*imd0);
  if (1 == fscanf(fp,"seed %li\n",seed)) printf("<checkpnt> seed %li\n",*seed); 
  if (1 == fscanf(fp,"orig %li\n",orig)) printf("<checkpnt> orig %li\n",*orig);
  if (1 == fscanf(fp,"N %i\n",&n)) printf("<checkpnt> N %i\n",n);
  if (1 == fscanf(fp,"NCR %i\n",&ncr)) printf("<checkpnt> NCR %i\n",ncr);
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
  int n=0;

  strcpy(str,fdir);
  strcat(str,fname);

  if ( (fp = fopen(str,"r")) == NULL) return 1;
  n+=fread(frdx,sizeof(double),N,fp);
  n+=fread(frdy,sizeof(double),N,fp);
  n+=fread(frdz,sizeof(double),N,fp);
  n+=fread(frcdx,sizeof(double),NCR,fp);
  n+=fread(frcdy,sizeof(double),NCR,fp);
  n+=fread(frcdz,sizeof(double),NCR,fp);
  fclose(fp);

  if (n != 3*N + 3*NCR)
    printf("Error reading forces\n");
  
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
  int n=0;

  strcpy(str,fdir);
  strcat(str,fname);

  if ( (fp = fopen(str,"r")) == NULL) return 1;
  n+=fread(vx,sizeof(double),N,fp);
  n+=fread(vy,sizeof(double),N,fp);
  n+=fread(vz,sizeof(double),N,fp);
  n+=fread(vxc,sizeof(double),NCR,fp);
  n+=fread(vyc,sizeof(double),NCR,fp);
  n+=fread(vzc,sizeof(double),NCR,fp);
  fclose(fp);

  if (n != 3*N + 3*NCR)
    printf("Error reading momenta\n");

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
  if (1 != fread(&Epot,sizeof(double),1,fp)) return 1;
  if (N != fread(x,sizeof(double),N,fp)) return 1;
  if (N != fread(y,sizeof(double),N,fp)) return 1;
  if (N != fread(z,sizeof(double),N,fp)) return 1;
  if (NCR != fread(xcr,sizeof(double),NCR,fp)) return 1;
  if (NCR != fread(ycr,sizeof(double),NCR,fp)) return 1;
  if (NCR != fread(zcr,sizeof(double),NCR,fp)) return 1;

  return 0;
} 
/****************************************************************************/
void check_distances(int *ip1,int *ip2,int n,double *dist2,int *dual,
		     double *xn,double *yn,double *zn) {
  /* Look for native contacts (ip1,ip2) with a distance (dist2) greater than
     the distance of that contact in another fold (xn,yn,zn). Care should be
     taken with such contacts as they may results in steric clashes. */
  int i,j,m;
  double r2;  

  if (FF_CONT < 2)
    return ;
  
  for (m = 0; m < n; m++) {
    if (dual[m] > 0) continue;

    i = ip1[m]; j = ip2[m];
    r2 = ( (xn[i] - xn[j]) * (xn[i] - xn[j]) +
	   (yn[i] - yn[j]) * (yn[i] - yn[j]) + 
	   (zn[i] - zn[j]) * (zn[i] - zn[j]) );

    if (dist2[m] > r2) 
      printf("<check_distances> cont %3i %3i (%3.5lf < %3.5lf). Caution: potential steric clash.\n",
	     i,j,sqrt(r2),sqrt(dist2[m]));

    /*    if (dist2[m] > r2) {
	  dist_rep[m] = r2;
	  fprintf(fp_log,"<correct_dist> %s cont %3i %3i (%3.5lf < %3.5lf). Setting repulsive distance to %lf\n",
	  dual[m] > 0 ? "Dual  " : "Native",i,j,sqrt(r2),sqrt(dist2[m]),sqrt(dist_rep[m]));
	  }  */
  }

  return ;
}
/****************************************************************************/
int get_shared_contacts(int *mc1,int *mc2,int *dual1,int *dual2) {
  int i,j,m,n,s = 0;

  if (FF_CONT < 2)
    return 0; 

  for (m = 0; m < npair; m++) {
    i = ip1[m]; j = ip2[m];

    for (n = 0; n < npair2; n++) {
      if ( (i == ip3[n] && j == ip4[n]) ||
	   (j == ip3[n] && i == ip4[n]) ) {
	mc1[s] = m;
	mc2[s] = n;
	dual1[m] = 1;
	dual2[n] = 1;
	s++;
      }
    }

  }

  return s;
}
/****************************************************************************/
void get_natdist(double *dist,double *dist2,int *ip1,int *ip2,int n,
		 double *xn,double *yn,double *zn) {
  /* Get native distances (dist) and distances squared (dist2) */
  int i,j,m;
  double r2;

  for (m = 0; m < n; m++) {
    i = ip1[m]; j = ip2[m]; 
    r2 = ( (xn[i]-xn[j]) * (xn[i]-xn[j]) +
	   (yn[i]-yn[j]) * (yn[i]-yn[j]) + 
	   (zn[i]-zn[j]) * (zn[i]-zn[j]) );
    dist[m] = sqrt(r2);
    dist2[m] = r2;
  }
  
  return ;
}
/****************************************************************************/
void get_nndist(double *distg1,double *distg2,int *ip1,int *ip2,int n,
		int *nni1,int *nnj1,int *nni2,int *nnj2,
		double *xn,double *yn,double *zn) {
  /* Get nearest neighbor distances */
  int i,j,ic,jc,i1,j1,m,q,imin=0,jmin=0;
  double gmin,dg[4];  

  gmin = dg[0] = dg[1] = dg[2] = dg[3] = 0;

  for (m = 0; m < n; m++) nni1[m] = nnj1[m] = nni2[m] = nnj2[m] = -1;
  
  for (m = 0; m < n; m++) {
    i = ip1[m]; j = ip2[m]; 
    ic = a2c[i]; jc = a2c[j];
    
    q = 0;
    gmin = 1e6;
    for (i1 = i-1; i1 <= i+1; i1 += 2) {
      for (j1 = j-1; j1 <= j+1; j1 += 2) {
	
	if (i1 < iBeg[ic] || i1 > iEnd[ic] || j1 < iBeg[jc] || j1 > iEnd[jc]) continue;

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

  return ;
}
/****************************************************************************/
void write_shared_dist(char *fn,int *mc1,int *mc2,int spair) {
  int i,j,l,k,m,n,s;
  char str[100];
  FILE *fp1;

  strcpy(str,TESTDIR);
  strcat(str,fn);

  printf("<write_shared_dist> Found %i shared contacts between %s and %s \n",
	 spair,CONTMAP,CONTMAP2);    
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
int read_native(char *fn,double *xr,double *yr,double *zr) {
  int j, n = 0;
  double tmpx,tmpy,tmpz;
  FILE *fp1;

  if ( NULL == (fp1 = fopen(fn,"r")) ) return 0;

  while (4 == fscanf(fp1,"%i %lf %lf %lf",&j,&tmpx,&tmpy,&tmpz)) {

    if (j < 0 || j > N-1) {
      fprintf(fp_log,"<read_native> (%s) Ignoring j = %i \n",fn,j);
      continue;
    }

    xr[j] = tmpx;
    yr[j] = tmpy;
    zr[j] = tmpz;

    ++n;
  }    

  fclose(fp1);
  
  return n;
}
/****************************************************************************/
int read_contacts(char *fn,int *ip1,int *ip2) {
  int n = 0;
  FILE *fp1;

  if ( NULL == (fp1 = fopen(fn,"r")) )
    return 0;
  
  while (2 == fscanf(fp1,"%i %i",ip1 + n,ip2 + n)) {
    if (ip1[n] < 0 || ip1[n] > N-1 || ip2[n] < 0 || ip2[n] > N-1) 
      fprintf(fp_log,"<read_contacts> (%s) Ignoring %3d, %3d\n",
	      fn,ip1[n],ip2[n]);
    else 
      n++;
  }
  
  fclose(fp1);

  return n;
}  
/****************************************************************************/
void write_bonded_param(double *bn,double *thn,double *phn,char *fn) {
  int j;
  char str[100];
  FILE *fp;

  strcpy(str,TESTDIR);
  strcat(str,fn);

  fp = fopen(str,"w");
  for (j = 0; j < N; j++)  
    fprintf(fp,"%3i BOND %8.6lf ANGLE %8.6lf TORS %8.6lf\n",
	    j,bn[j],thn[j]*rad2deg,phn[j]*rad2deg);
  fclose(fp); 
}
/****************************************************************************/
void set_bonded_param(double *bn,double *thn,double *phn,
		      double *xnat,double *ynat,double *znat,int n) {
  int i;

  printf("<set_bonded_param> Setting bonded parameters \n");

  /* default */
  
  for (i = 0; i < N-1; i++) bn[i] = 3.8; 
  for (i = 1; i < N-1; i++) thn[i] = 120.0*deg2rad;
  for (i = 1; i < N-2; i++) phn[i] = 120.0*deg2rad;

  if (n < N)
    return ;

  /* native */

  for (i = 0; i < N; i++) {x[i] = xnat[i]; y[i] = ynat[i]; z[i] = znat[i];}

  cart2dof();

  for (i = 0; i < N-1; i++)  bn[i] = b[i];
  for (i = 1; i < N-1; i++)  thn[i] = th[i];  
  for (i = 1; i < N-2; i++)  phn[i] = ph[i];

  // if (1 != check_chains()) printf("<set_bonded_param> Error native configuration\n");
}
/****************************************************************************/
void read_disregs(char fn[],int *dis) {
  int i,j;
  FILE *fp = fopen(fn,"r");
  
  if (NULL == fp) 
    return;

  while (2 == fscanf(fp,"%i %i",&i,&j) && !feof(fp)) 
    if (i >= 0 && i <= N-1 && j > 0) {
      printf("<read_dis_regions> Setting res %i as disordered (%s)\n",i,fn);
      dis[i] = 1;
    }

  fclose(fp);  
}
/****************************************************************************/
int relax_chains(int ich) {
  /* ich >= 0 : relax chain ich */
  /* ich <  0 : relax all chains */
  int i,icur,n = 0;
  double rfac = 1.0;
  double Erel,Eold,pho[N];
  double dx,dy,dz;

  if (N == 0)
    return 0;
  
  for (i = 1; i < N-2; i++) pho[i] = ph[i];
  
  if (ich < 0)
    printf("<init> Relaxing all chains... \n");
  else
    printf("<init> Relaxing chain %i... \n",ich);
    
  exvol(0,&Erep,&Ehp,&Eel); 
 
  Erel = Eold = Erep + cont(0);
 
  while (Erel > N * rfac * eps)  {

    icur = (ich < 0 ? NCH * ran3n(&seed) : ich);
      
    if ( ran3n(&seed) < 0.5 ) { /* turn torsion angle */

      i = iBeg[icur] + ran3n(&seed) * (iEnd[icur] - iBeg[icur] + 1);
      ph[i] =  pi * (2 * ran3n(&seed) - 1);
      dof2cart(0);
	    
      exvol(0,&Erep,&Ehp,&Eel); 
      
      if ( (Erel + cont(0)) < Eold ) {
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
      exvol(0,&Erep,&Ehp,&Eel); 
      if ((Erel = Erep + cont(0)) < Eold) {
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
  
  if (NCR == 0)
    return 0;
  
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
      
      if (n % 10000 == 0)  printf("Erel %lf\n",Eold);

      n++;
  } 
  
  printf("<init> ...done in %i steps (Erel %lf) \n",n,Erel);

  return 0;
}
/****************************************************************************/
//–– Convert “123.4” or “NCS1”/“NCS2” into a double
static double lookup_strength(const char *tok) {
    // numeric literal?
    if (isdigit(tok[0]) || tok[0]=='-' || tok[0]=='.') {
        return atof(tok);
    }
    // symbolic constants
    if (strcmp(tok, "NCS1")==0) return (double) NCS1;
    if (strcmp(tok, "NCS2")==0) return (double) NCS2;
    fprintf(stderr, "Unknown contact strength '%s'\n", tok);
    exit(EXIT_FAILURE);
}
/****************************************************************************/
void read_cont_param(char fn[],
                     int  ip1[],
                     int  ip2[],
                     int  npair,
                     double kcont[])
{
    int    i, j, m;
    double kread;
    char   tok[32];
    FILE  *fp;

    if ((fp = fopen(fn, "r")) != NULL) {
        while (fscanf(fp, "%i %i %31s\n", &i, &j, tok) == 3) {
            kread = lookup_strength(tok);

            for (m = 0;  m < npair;  ++m) {
                if ((i == ip1[m] && j == ip2[m]) ||
                    (i == ip2[m] && j == ip1[m])) {
                    break;
                }
            }
            if (m == npair) {
                printf("<read_contpar> (%s) unknown contact %i %i %i\n",
                       fn, m, i, j);
                continue;
            }

            kcont[m] = kread;
            printf("<read_contpar> (%s) setting strength of contact %i %i %i to %lf\n",
                   fn, m, i, j, kcont[m]);
        }
	fclose(fp);
    }
}
/****************************************************************************/
void cont_param_salt(char fn[],int ip1[],int ip2[],int n,double kcont[]) {
  int i,j,m;
  double fac;
  char str[100];
  FILE *fp;

  strcpy(str,TESTDIR);
  strcat(str,"salt_fac_");
  strcat(str,fn);

  printf("<cont_param_salt> Writing to %s\n",str);
  fp = fopen(str,"w");
  for (m = 0; m < n; m++) {
    i = ip1[m];
    j = ip2[m];
    
    if (qres[i]*qres[j] == 0) continue;

    fac = csalt_fac(csalt,qres[i],qres[j]);
    kcont[m] *= fac;

    fprintf(fp,"%3i %3i %3i %c %c qres %3lf %3lf saltfac %8.5lf kcon %8.5lf\n",
	    m,i,j,seq[i],seq[j],qres[i],qres[j],fac,kcont[m]);
  }
  fclose(fp);
}
/****************************************************************************/
static int parse_c2p_from_string(const char *spec) {
    // spec format: "0,0,0,0,1" (NCH comma-separated nonnegative ints)
    c2p = (int*)malloc(NCH * sizeof *c2p);
    int c = 0;
    const char *p = spec;
    while (*p && c < NCH) {
        int id = 0, got = 0;
        while (*p >= '0' && *p <= '9') { id = id*10 + (*p - '0'); p++; got = 1; }
        if (!got) return 0;
        c2p[c] = id;
        if (NPROT <= id) NPROT = id + 1;
        c++;
        if (*p == ',') p++;
    }
    return (c == NCH);
}

void build_protein_maps(const char *spec) {
    if (a2p) return;                        // already built

    if (spec) {
        if (!parse_c2p_from_string(spec)) {
            fprintf(stderr,"[prot] Bad c2p spec: \"%s\"\n", spec);
            exit(1);
        }
    } else {
        // fallback: each chain is its own protein
        c2p = (int*)malloc(NCH * sizeof *c2p);
        for (int c=0; c<NCH; c++) c2p[c] = c;
        NPROT = NCH;
    }

    a2p = (int*)malloc(N * sizeof *a2p);
    for (int i=0; i<N; i++) a2p[i] = c2p[a2c[i]];

    fprintf(stdout, "<prot> NCH=%d  NPROT=%d  a2p built\n", NCH, NPROT);
}
/****************************************************************************/
/***** INITIALIZATION *******************************************************/
/****************************************************************************/
void init(int iflag) {
  int i,j,k,m,nnat1,nnat2,n;
  double vxsum,vysum,vzsum,vxcsum,vycsum,vzcsum;
  double c0,c0cr;
  double o[NOBS],fdum[MAXP];
  char c,str[100];
  FILE *fp1;

  /* read sequence */

  fp1 = fopen(INPUT,"r");
  k = 0;
  for (j = 0; j < NCH; j++) {
    iBeg[j] = k;
    while ( (c = getc(fp1)) != '\n' && feof(fp1) == 0 ) {
      if (c >= 'A' && c <= 'z') {
	seq[k] = c;
	a2c[k] = j;
	switch (c) {
	case 'D' : {qres[k] = -1.0;  break;}
	case 'E' : {qres[k] = -1.0;  break;}
	case 'K' : {qres[k] = +1.0;  break;}
	case 'R' : {qres[k] = +1.0;  break;}
	default  : {qres[k] =  0.0;  break;}
	}
  switch (c) {
  case 'A' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'V' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'L' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'I' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'M' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'W' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'F' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  case 'Y' : {seqhp[k] = 1; iqkap[k] = 1.0; break;}
  default  : {seqhp[k] = 0; iqkap[k] = 0.0; break;}
  }
	++k;
      }
    }
    iEnd[j] = k-1;
  }
  fclose(fp1);
  build_protein_maps(C2P_SPEC);

  /* read g parameters */

  if (IREADG == 1 && NTMP > 1) {
    fp1 = fopen(INPUTG,"r");
    if (fp1 == NULL) {
      printf("<init> No file %s\n",INPUTG);
      for (i=0; i<NTMP; i++) g[i] = 0.0;
    } else {
      for (i=n=0; i<NTMP; i++) n += fscanf(fp1,"%i %lf", &i, &g[i]);
      printf("<init> %i lines read in %s\n",n/2,INPUTG);
      fclose(fp1);
    }
  }

  /* output directory */

  strcpy(OUTDIR,RESDIR);
  
  /* set seed */ 

  if (ISEED == 1) {
    FILE *devrandom = fopen("/dev/urandom", "r");
    n=fread(&seed, sizeof(seed), 1, devrandom);
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
    for (j = 0; j < NCH; j++) {
      printf("%3d %3d %3d ",j,iBeg[j],iEnd[j]);
      for (i = iBeg[j]; i <= iEnd[j]; i++) printf("%c",seq[i]); 
      printf("\n");
    }
    printf("\n");
  }
  
  printf("Index  temp  beta  g:\n");
  for (i = 0; i < NTMP; i++) {
    printf("%i %lf %lf %lf\n",i,1./beta[i], beta[i], g[i]);
  }
  printf("\n");

  /* create directories */

  printf("<init> Creating directory %s\n",RESDIR);
  mkdir(RESDIR, 0777);
  printf("<init> Creating directory %s\n",TESTDIR);
  mkdir(TESTDIR, 0777);
  printf("<init> Creating directory %s\n",CHECKDIR);
  mkdir(CHECKDIR, 0777);

  strcpy(str,OUTDIR);
  strcat(str,LOGFILE);
  printf("<init> Opening %s\n",str);
  fp_log = fopen(str,"a");

  /* constants */
  
  pi=acos(-1.);
  pi2=2.*pi;
  pid2=pi/2.;
  deg2rad=pi/180.;
  rad2deg=180./pi;
  cthlim=cos(pi/180.);

  /* energy parameters */

  kbon*=eps;
  kth*=eps;
  kph1*=eps;
  kph3*=eps;
  kcon*=eps;
  khp*=eps;
  kel*=eps;
  krep*=eps;
  eclash*=eps;

  thn_disa*=deg2rad;
  thn_disb*=deg2rad;
  ksi_disa*=deg2rad;
  ksi_disb*=deg2rad;
  phn_dis1*=deg2rad;
  phn_dis2*=deg2rad;
  phn_dis3*=deg2rad;
  kph_dis1*=eps;
  kph_dis2*=eps;
  kph_dis3*=eps;
  eth0*=eps;
  eph0*=eps;
  
  /* set temperatures */
  
  ind = NTMP * ran3n(&seed);
  printf("<init> Setting temperature index to %i\n",ind);

  /* MD parameters */

  /*************** Beads  **************/
  dt=0.005*tau;
  gam=0.05/tau;
  mbd=1.0;
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
  
  if (iflag == 0)
    return ;

  /* native structures */

  nnat1 = read_native(NATIVE,xnat,ynat,znat);
  printf("<init> NATIVE:  Found %3i residue positions in %s\n",nnat1,NATIVE);

  nnat2 = read_native(NATIVE2,xnat2,ynat2,znat2);
  printf("<init> NATIVE2: Found %3i residue positions in %s\n",nnat2,NATIVE2);

  /* contacts */

  npair = npair2 = spair = ndpair = 0;
  for (m = 0; m < MAXP; m++) 
    ip1[m] = ip2[m] = ip3[m] = ip4[m] = id1[m] = id2[m] = 0;
  for (i = 0; i < N; i++) for (j = 0;j < N; j++) cc[i][j] = 0;  
  for (m = 0; m < MAXP; m++) dual1[m] = dual2[m] = 0;    

  npair = read_contacts(CONTMAP,ip1,ip2);
  printf("<init> CONTMAP:  Found %i contacts in %s\n",npair,CONTMAP);
  get_natdist(distp,distp2,ip1,ip2,npair,xnat,ynat,znat);
  printf("<init> CONTMAP:  Writing to %s\n","contact_map.out");
  write_natdist("contact_map.out",distp,npair,ip1,ip2);

  npair2 = read_contacts(CONTMAP2,ip3,ip4);
  printf("<init> CONTMAP2: Found %i contacts in %s\n",npair2,CONTMAP2);
  get_natdist(distp3,distp4,ip3,ip4,npair2,xnat2,ynat2,znat2);
  printf("<init> CONTMAP2: Writing to %s\n","contact_map2.out");
  write_natdist("contact_map2.out",distp3,npair2,ip3,ip4);

  if (FF_DISULF) {
    ndpair = read_contacts(DISULFIDE,id1,id2);
    printf("<init> DISULFIDE: FF_DISULF %i \n",FF_DISULF);
    printf("<init> DISULFIDE: Found %i contacts in %s \n",ndpair,DISULFIDE);
    get_natdist(distd1,fdum,id1,id2,ndpair,xnat,ynat,znat);
    get_natdist(distd2,fdum,id1,id2,ndpair,xnat2,ynat2,znat2);
    for (n = 0; n < ndpair; n++) 
      printf("<init> DISULFIDE: %4i %4i nat %8.5lf   nat2 %8.5lf\n",id1[n],id2[n],distd1[n],distd2[n]);
  }

  if (npair  > MAXP) {printf("npair too big\n"); exit(-1);}
  if (npair2 > MAXP) {printf("npair2 too big\n"); exit(-1);}
  if (ndpair > MAXP) {printf("ndpair too big\n"); exit(-1);}

  /* NON-BONDED INTERACTIONS */

  if (FF_CONT > 0) {
    for (m = 0; m < npair; m++) kcon_nat[m] = kcon;
    get_nndist(distg1,distg2,ip1,ip2,npair,nni1,nnj1,nni2,nnj2,xnat,ynat,znat);
    for (m = 0; m < npair; m++) {
      i = ip1[m]; j = ip2[m];
      cc[i][j] = cc[j][i] = 1;
    }
  }
  
  if (FF_CONT == 2) {
    for (m = 0; m < npair2; m++) kcon_nat2[m] = kcon;
    get_nndist(distg3,distg4,ip3,ip4,npair2,nni3,nnj3,nni4,nnj4,xnat2,ynat2,znat2);
    for (m = 0; m < npair2; m++) {
      i = ip3[m]; j = ip4[m];
      cc[i][j] = cc[j][i] = 1;
    }
    spair = get_shared_contacts(mc1,mc2,dual1,dual2);
    write_shared_dist("common_contacts.out",mc1,mc2,spair);
    printf("Looking for potential clashes in NATIVE2 (%s) due to dual basin potential...\n",NATIVE2);
    check_distances(ip1,ip2,npair,distp2,dual1,xnat2,ynat2,znat2);
    printf("Looking for potential clashes in NATIVE (%s) due to dual basin potential...\n",NATIVE);
    check_distances(ip3,ip4,npair2,distp4,dual2,xnat,ynat,znat);
  }

  /* parameters from file */  

  read_cont_param(CONTPAR,ip1,ip2,npair,kcon_nat);
  read_cont_param(CONTPAR2,ip3,ip4,npair2,kcon_nat2);

  /* Effective screening of contacts between charged residues (salt effect) */
  
  if (FF_SALT) {
    printf("<init> FF_SALT %i \n",FF_SALT);
    printf("<init> csalt: %.5lf\n",csalt);
    printf("<init> Scaling contact strengths...\n");
    cont_param_salt("1",ip1,ip2,npair,kcon_nat);
    cont_param_salt("2",ip3,ip4,npair2,kcon_nat2);
  }
  
  /* BONDED INTERACTIONS */
  
  set_bonded_param(bn,thn,phn,xnat,ynat,znat,nnat1);
  set_bonded_param(bn2,thn2,phn2,xnat2,ynat2,znat2,nnat2);

  write_bonded_param(bn,thn,phn,"bonded_param_1");
  write_bonded_param(bn2,thn2,phn2,"bonded_param_2");
  
  /* disordered regions from file */

  read_disregs(DISREG,dis);
  read_disregs(DISREG2,dis2);

  /* initialize functions */
  
  double tmpf;
  
  dof2cart(-1);
  crowd_crowd(-1);
  crowd_bead(-1);
  bond(-1);
  bend(-1);
  torsion(-1);
  exvol(-1,&tmpf,&tmpf,&tmpf);
  cont(-1);
  hp(-1,0,0,0,0,0,0);
  el(-1,0,0,0,0,0,0);

  read_conf(-1,"","");
  histo_bond(-1);
  histo_bend(-1);
  histo_tors(-1,0);
  histoe(-1,0);
  histo_cont1(-1,0,0);
  histo_cont2(-1,0,0);

  /* initial chain configuration  */  

  if (ISTART == 0) {            /* native */

    printf("<init> Initializing chain(s) from NATIVE %s\n",NATIVE);
    for (i = 0; i < N; i++) {x[i] = xnat[i]; y[i] = ynat[i]; z[i] = znat[i];}
    cart2dof();
    //    if (1 != check_chains()) printf("Error initial configuration"); 

  } else if (ISTART == 1) {     /* read */

    printf("<init> Initializing chain(s) from START %s\n",START);
    read_native(START,x,y,z);
    cart2dof();
    //    if (1 != check_chains()) printf("Error initial configuration 2\n");

  } else if (ISTART == 2) {     /* random */

    printf("<init> Initializing chain(s) from random configuration\n");
    for (i = 0; i < N-1; i++) b[i]  = bn[i] ;
    for (i = 1; i < N-1; i++) th[i] = thn[i] ;
    for (i = 1; i < N-2; i++) ph[i] = pi * (2 * ran3n(&seed) - 1);
    dof2cart(0);
    relax_chains(-1);

  } else {
    printf("\nInvalid ISTART\n");
    exit(-1);
  }

  /* Initial crowder positions (random) */
  
  for (i = 0; i < NCR; i++){
    xcr[i] = BOX * ran3n(&seed);
    ycr[i] = BOX * ran3n(&seed);
    zcr[i] = BOX * ran3n(&seed);
  }
  
  relax_crowders();

  /* initialize velocities */
  
  vxsum = vysum = vzsum = 0;
  for (i = 0; i < N; i++) {
    vx[i] = sqrt(1./beta[ind]/mbd) * gasdev2();
    vy[i] = sqrt(1./beta[ind]/mbd) * gasdev2();
    vz[i] = sqrt(1./beta[ind]/mbd) * gasdev2();
    vxsum += vx[i];
    vysum += vy[i];
    vzsum += vz[i];
  }

  vxcsum = vycsum = vzcsum = 0;
  for (i = 0; i < NCR; i++) {
    vxc[i] = sqrt(1./beta[ind]/mcr) * gasdev2();
    vyc[i] = sqrt(1./beta[ind]/mcr) * gasdev2();
    vzc[i] = sqrt(1./beta[ind]/mcr) * gasdev2();
    vxcsum += vxc[i];
    vycsum += vyc[i];
    vzcsum += vzc[i];
  }

  Ekin = 0;
  for (i = 0; i < N; i++) {
    vx[i] -= vxsum / max(N,1);
    vy[i] -= vysum / max(N,1);
    vz[i] -= vzsum / max(N,1);
    Ekin += vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  for (i = 0; i < NCR; i++) {
    vxc[i] -= vxcsum / max(NCR,1);
    vyc[i] -= vycsum / max(NCR,1);
    vzc[i] -= vzcsum / max(NCR,1);
    Ekin += vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i];
  }
  Ekin *= 0.5;
  
  /* initial energies and forces */
  
  for (i = 0; i < N; i++) fx[i]=fy[i]=fz[i]=0;
  for (i = 0; i < NCR; i++) fxc[i]=fyc[i]=fzc[i]=0;

  exvol(0,&Erep,&Ehp,&Eel); 
 
  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep)+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp)+(Eel)+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0));

  for (i = 0; i < NCR; i++) {
    frcdx[i] = gasdev2() * tconstcr[ind];
    frcdy[i] = gasdev2() * tconstcr[ind];
    frcdz[i] = gasdev2() * tconstcr[ind];
  }
  
  for (i = 0; i < N; i++) {
    frdx[i] = gasdev2() * tconstbd[ind];
    frdy[i] = gasdev2() * tconstbd[ind];
    frdz[i] = gasdev2() * tconstbd[ind];
  }

  if (!FANALYS && read_checkpnt()) return ;

  printf("<init> writing initial conformation to %s\n","_start.pdb");
  dumppdb("_start.pdb",o,0);

  printf("\n");
  printf("Initial conformation (%d chains, %d crowders):\n",NCH,NCR);
  printf("  Ekin %f Epot %f\n",Ekin,Epot);
  printf("  Ebon %f Eben %f Erep %f Etor %f Econ %f Econ1 %lf Econ2 %lf Ecorr %lf\n",
	 Ebon,Eben,Erep,Etor,Econ,Econ1,Econ2,Ecorr);
  printf("  Ehp %f Eel %f Ecc %f Ecb %f\n",Ehp,Eel,Ecc,Ecb);

  /* write potential functions to file */
  
  bond(1);
  bend(1);
  torsion(1);
  cont_corr(1);
  cont(1);
  cont2(1);
  crowd_crowd(1);
  crowd_bead(1); 
 // hp(1,0,0);
 // el(1);
  
  return ;
}
/****************************************************************************/
