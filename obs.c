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
/***** OBSERVABLES **********************************************************/
/****************************************************************************/
double dist_disulf(int n) {
  return (n > ndpair-1 ? -1 : sqrt( dist2(id1[n],id2[n]) ) );
}
/****************************************************************************/
int no_cont(int ch) {
  /* Number of NATIVE contacts in chain ch (>=0) or any chain (<0) */
  int i,j,m,n = 0;

  for (m = 0; m < npair; m++) {
    i = ip1[m]; j = ip2[m];
    if ( (a2c[i] == ch && a2c[j] == ch) || ch < 0 ) 
      n += ( dist2(i,j) < 1.44 * distp2[m] );
  }
  
  return n;
}
/****************************************************************************/
int no_cont2(int ch) {
  /* Number of NATIVE2 contacts in chain ch (>=0) or any chain (<0) */
  int i,j,m,n = 0;

  for (m = 0; m < npair2; m++) {
    i = ip3[m]; j = ip4[m];
    if ( (a2c[i] == ch && a2c[j] == ch) || ch < 0) 
      n += ( dist2(i,j) < 1.44 * distp4[m] );
  }
  
  return n;
}
/****************************************************************************/
int no_cont_ch2ch(int ch1,int ch2) {
  /* Number of NATIVE contacts between chains ch1 and ch2 */
  int i,j,m,n = 0;
  
  for (m = 0; m < npair; m++) {
    i = ip1[m]; j = ip2[m];
    if ((a2c[i] == ch1 && a2c[j] == ch2) || (a2c[j] == ch1 && a2c[i] == ch2))
      n += ( dist2(i,j) < 1.44 * distp2[m] );
  }
  
  return n;
}
/****************************************************************************/
int no_cont2_ch2ch(int ch1,int ch2) {
  /* Number of NATIVE2 contacts between chains ch1 and ch2 */
  int i,j,m,n = 0;
  
  for (m = 0; m < npair2; m++) {
    i = ip3[m]; j = ip4[m];
    if ((a2c[i] == ch1 && a2c[j] == ch2) || (a2c[j] == ch1 && a2c[i] == ch2))
      n += ( dist2(i,j) < 1.44 * distp4[m] );
  }
  
  return n;
}
/****************************************************************************/
void histo_contmap(int iflag, int ind) {
  static double his[N][N],n;
  int i, j,m;
  FILE *fp;

  if (iflag<0) {
    for (m=0;m<npair;m++){
      i=ip1[m]; j=ip2[m];
      his[i][j]=n=0;
    }
    fp=fopen("results/_his_contmap","w");
    fclose(fp);
    return;
  }
  
  if (ind == NTMP-1) {  
    if (iflag==0) {   
      for (m=0;m<npair;m++) {
	i=ip1[m]; j=ip2[m];
	if ( ((x[i]-x[j])*(x[i]-x[j])+
	      (y[i]-y[j])*(y[i]-y[j])+
	      (z[i]-z[j])*(z[i]-z[j])) <= 1.44 * distp2[m]) {
	  his[i][j]++;
	  his[j][i]++;
	}
      }
      n++;
      return;
    }
  }
  if (iflag>0) {
    fp=fopen("results/_his_contmap","w");
    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
	fprintf(fp,"%i %i %lf\n", i, j, his[i][j]/n);
      }
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histo_contmap2(int iflag, int ind)
{
  static double his[N][N],n;
  int i, j,m;
  FILE *fp;
  
  if (iflag<0) {
    for (m=0;m<npair2;m++){
      i=ip3[m]; j=ip4[m];
      his[i][j]=n=0;
    }
    fp=fopen("results/_his_contmap2","w");
    fclose(fp);
    return;
  }
  
  if (ind == NTMP-1) {
    if (iflag==0) { 
      for (m=0;m<npair2;m++) {
	i=ip3[m]; j=ip4[m];
	if ( ((x[i]-x[j])*(x[i]-x[j])+
	      (y[i]-y[j])*(y[i]-y[j])+
	      (z[i]-z[j])*(z[i]-z[j])) <= 1.44 * distp4[m]) {
	  his[i][j]++;
	  his[j][i]++;
	}
      }
      n++;
      return;
    }
  }
  
  if (iflag>0) {
    fp=fopen("results/_his_contmap2","w");
    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
	fprintf(fp,"%i %i %lf\n", i, j, his[i][j]/n);
      }
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
double gyr2(int a1,int a2) {
  int i;
  double rmx,rmy,rmz,rg2;

  rmx=rmy=rmz=rg2=0;
  for (i=a1;i<=a2;i++) {
    rmx+=x[i];
    rmy+=y[i];
    rmz+=z[i];
  }
  rmx/=(a2-a1+1); rmy/=(a2-a1+1); rmz/=(a2-a1+1);
  for (i=a1;i<=a2;i++) {
    rg2+=(x[i]-rmx)*(x[i]-rmx)+
         (y[i]-rmy)*(y[i]-rmy)+
         (z[i]-rmz)*(z[i]-rmz);
  }
  return rg2/(a2-a1+1);
}
/****************************************************************************/
void center_of_mass2(int ich,double *xcm,double *ycm, double *zcm) {
  int i,n = 0;

  (*xcm) = (*ycm) = (*zcm) = 0;

  for (i = iBeg[ich]; i <= iEnd[ich]; i++) {
    (*xcm) += x[i];
    (*ycm) += y[i];
    (*zcm) += z[i];
    ++n;
  }

  (*xcm) /= n; (*ycm) /= n; (*zcm) /= n; 

  return ;
}
/****************************************************************************/
void histo_bond(int iflag) {
  static double low = 2.5,high = 4.5;
  static double his[N][NBIN],eps;
  static int out=0;
  double x,sum;
  int i,j,k;
  char str[100];
  FILE *fp;
  
  if (iflag < 0) {
    for (i=0;i<N;i++) for (j=0; j<NBIN; j++) his[i][j] = 0;
    eps = (high - low) / NBIN;

    strcpy(str,OUTDIR);
    strcat(str,"_his_bond");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag == 0 && ind == NTMP-1) {
    for (i = 0; i < N; i++) {
      x = b[i];
      if (x > low && x < high) {
	k = (x - low) / eps;
	his[i][k]++;
      } else
	out++;
    }
    return;
  }

  if (iflag > 0) {
    strcpy(str,OUTDIR);
    strcat(str,"_his_bond");

    fp = fopen(str,"w");
    fprintf(fp,"# out %i\n",out);
    for (i = 0; i < N-1; i++) {
      sum = 0;
      for (k = 0; k < NBIN; k++) sum += his[i][k];
      for (k = 0; k < NBIN; k++) {
	x = low + (k + .5) * eps;
	fprintf(fp,"%i %f %e\n",i,x,his[i][k] / sum / eps);
      }
    }
    fclose(fp);
    return;
  }

  return;
}
/****************************************************************************/
void histo_bend(int iflag) {
  static double low = 40,high = 180;
  static double his[N][NBIN],eps;
  static int out = 0;
  double x,sum;
  int i,j,k;
  char str[100];
  FILE *fp;

  if (iflag < 0) {
    for (i=0; i<N; i++) for (j=0; j<NBIN; j++) his[i][j] = 0;
    eps = (high - low) / NBIN;

    strcpy(str,OUTDIR);
    strcat(str,"_his_bend");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag == 0 && ind == NTMP-1) {
    for (i = 1; i < N-1; i++) {
      x = th[i] * 180 / pi;
      if (x > low && x < high) {
	k = (x - low) / eps;
	his[i][k]++;
      } else
	out++;
    }
    return;
  }

  if (iflag > 0) {
    strcpy(str,OUTDIR);
    strcat(str,"_his_bend");

    fp = fopen(str,"w");
    fprintf(fp,"# out %i\n",out);
    for (i = 1; i < N-1; i++) {
      sum = 0;
      for (k=0; k<NBIN; k++) sum += his[i][k];
      for (k=0; k<NBIN; k++) {
	x = low + (k + .5) * eps;
	fprintf(fp,"%i %f %e\n",i,x,his[i][k] / sum / eps);
      }
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histo_tors(int iflag,int ia) {
  static double low=-180,high=180;
  static double his[NTMP][NBIN],eps;
  static int out[NTMP];
  double x,sum;
  int i,j,k;
  char str[100];
  FILE *fp;

  if (iflag < 0) {
    for (i = 0; i < NTMP; i++) for (j=0; j<NBIN; j++) his[i][j] = out[i] = 0;
    eps = (high - low) / NBIN;

    strcpy(str,OUTDIR);
    strcat(str,"_his_tors");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag == 0) {

    x = ph[ia];
    while (x < -pi) x += pi2;
    while (x > +pi) x -= pi2;
    x *= rad2deg;

    if (x >= low && x < high) {
	k = (x - low) / eps;
	his[ind][k]++;
    } else
      out[ind]++;
    return;
  }

  if (iflag > 0) {
    strcpy(str,OUTDIR);
    strcat(str,"_his_tors");

    fp = fopen(str,"w");
    for (i = 0; i < NTMP; i++) {
      sum = 0;
      for (k = 0; k < NBIN; k++) sum += his[i][k];
      for (k = 0; k < NBIN; k++) {
	x = low + (k + .5) * eps;
	fprintf(fp,"%i %f %e\n", i, x, his[i][k] / sum / eps);
      }
    }

    for (i = 0; i < NTMP; i++)
      fprintf(fp,"# out temp %i %i\n",i, out[i]);

    fclose(fp);
    return;
  }
  
  return;
}
/****************************************************************************/
void histoe(int iflag,double x) {
  static double low=-100.0,high=200.0;
  static double his[NBIN],eps;
  int j,k;
  char str[100];
  FILE *fp;

  if (iflag<0) {
    for (j=0;j<NBIN;j++) his[j]=0;
    eps=(high-low)/NBIN;
    strcpy(str,OUTDIR);
    strcat(str,"_hise");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag==0) {
    if (x>low && x<high) {
      k=(x-low)/eps;
      his[k]++;
    } 
    return;
  }

  if (iflag>0) {
    fp=fopen("_hise","w");
    j=0; for (k=0;k<NBIN;k++) j+=his[k];
    for (k=0;k<NBIN;k++)
      if (j!=0) fprintf(fp,"%f %.15e\n",low+(k+.5)*eps,his[k]/j/eps);
    fclose(fp);
  }
}
/****************************************************************************/
void histormsd(int iflag,double x) {
  static double low=0.0,high=40.0;
  static double his[NBIN],eps;
  int j,k;
  char str[100];
  FILE *fp;

  if (iflag<0) {
    for (j=0;j<NBIN;j++) his[j]=0;
    eps=(high-low)/NBIN;
    strcpy(str,OUTDIR);
    strcat(str,"_hisrmsd");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag==0) {
    if (x>low && x<high) {
      k=(x-low)/eps;
      his[k]++;
    }
    return;
  }

  if (iflag>0) {
    strcpy(str,OUTDIR);
    strcat(str,"_hisrmsd");

    fp = fopen(str,"w");
    j=0; for (k=0;k<NBIN;k++) j+=his[k];
    for (k=0;k<NBIN;k++)
      if (j!=0) fprintf(fp,"%f %.15e\n",low+(k+.5)*eps,his[k]/j/eps);
    fclose(fp);
  }
}
/****************************************************************************/
void histoq1(int iflag, int x) {
  static int nc,his[MAXP];
  int j;

  char str[100];
  FILE *fp;

  if (iflag<0) {
    for (j=0;j<MAXP;j++) his[j]=0;
    nc=npair;
    if (nc>MAXP) {printf("histoq1: nc too big\n");exit(-1);}

    strcpy(str,OUTDIR);
    strcat(str,"_hisq1");
    fp = fopen(str,"w");
    fclose(fp);

    return;
  }
  
  if (iflag==0) {
    if (x >= 0 && x <= nc)
      his[x]++;
    return;
  }

  if (iflag>0) {
    strcpy(str,OUTDIR);
    strcat(str,"_hisq1");
    fp = fopen(str,"w");
    for (j = 0;j <= nc; j++) fprintf(fp,"%i %i\n",j,his[j]);
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void historg(int iflag,double x) {
  static double low=5.0,high=40.0;
  static double his[NBIN],eps;
  int j,k;
  char str[100];
  FILE *fp;

  if (iflag<0) {
    for (j=0;j<NBIN;j++) his[j]=0;
    eps=(high-low)/NBIN;
    strcpy(str,OUTDIR);
    strcat(str,"_hisrg");

    fp = fopen(str,"w");
    fclose(fp);
    return;
  }
  
  if (iflag==0) {
    if (x>low && x<high) {
      k=(x-low)/eps;
      his[k]++;
    }
    return ;
  }

  if (iflag>0) {
    strcpy(str,OUTDIR);
    strcat(str,"_hisrg");

    fp = fopen(str,"w");
    j=0; for (k=0;k<NBIN;k++) j+=his[k];
    for (k=0;k<NBIN;k++)
      if (j!=0) fprintf(fp,"%f %.15e\n",low+(k+.5)*eps,his[k]/j/eps);
    fclose(fp);
  }
}
/****************************************************************************/
void histo_cont1(int iflag, int ind, int n){
  static double his[NTMP][MAXP+1];
  static double ene[NTMP][MAXP+1];
  static double sum;
  int i,j;
  FILE *fp;

  if (iflag<0){
    for(i = 0;i < NTMP;i++)
      for (j = 0; j < MAXP; j++) his[i][j] = ene[i][j] = 0;
    fp = fopen("results/_his_cont1","w"); 
    fclose(fp);
    return;
  }

  if (iflag==0){
    if (n >= 0 && n <= MAXP) {
      his[ind][n]++;
      ene[ind][n]+=Epot;
    }
    return;
  }
  
  if (iflag>0){
    fp=fopen("results/_his_cont1","w");
    for(i = 0; i < NTMP; i++){
      for(j = sum = 0; j <= MAXP; j++) sum += his[i][j];
      for(j = 0; j <= MAXP; j++){
	if(his[i][j]>0) fprintf(fp,"%i %i %f %f\n",i,j,
				his[i][j]/sum,ene[i][j]/his[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
    return ;
  }

}
/****************************************************************************/
void histo_cont2(int iflag, int ind, int n){
  static double his[NTMP][MAXP+1];
  static double ene[NTMP][MAXP+1];
  static double sum;
  int i,j;
  FILE *fp;
  
  if (iflag<0){
    for(i = 0;i < NTMP;i++)
      for (j = 0; j < MAXP; j++) his[i][j] = ene[i][j] = 0;
    fp = fopen("results/_his_cont2","w"); 
    fclose(fp);
    return;
  }

  if (iflag==0){
    if (n >= 0 && n <= MAXP) {
      his[ind][n]++;
      ene[ind][n]+=Epot;
    }
    return;
  }

  if (iflag>0){
    fp=fopen("results/_his_cont2","w");
    for(i = 0; i < NTMP; i++){
      for(j = sum = 0; j <= MAXP; j++) sum += his[i][j];
      for(j = 0; j <= MAXP; j++){
	if(his[i][j]>0) fprintf(fp,"%i %i %f %f\n",i,j,his[i][j]/sum,ene[i][j]/his[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  return;
}
/****************************************************************************/
void histod1(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=40.0;
  static double his[NTMP][NBIN],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBIN;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBIN;
    fp=fopen("results/_hisd1","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("results/_hisd1","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod2(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=40.0;
  static double his[NTMP][NBIN],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBIN;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBIN;
    fp=fopen("results/_hisd2","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("results/_hisd2","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod3(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=40.0;
  static double his[NTMP][NBIN],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBIN;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBIN;
    fp=fopen("results/_hisd3","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("results/_hisd3","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i hisd3 out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
