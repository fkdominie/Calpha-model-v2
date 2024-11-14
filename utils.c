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
double xcr[NCR],ycr[NCR],zcr[NCR];       /* Crowders coordinates            */
double xcro[NCR], ycro[NCR], zcro[NCR];  /* Back up for crowders coordinates*/  
/****************************************************************************/
/***** input/output *********************************************************/
/****************************************************************************/
int read_native(char *fn,double *xr,double *yr,double *zr,int *nat) {
  int i, j, n = 0;
  FILE *fp1;

  fp1 = fopen(fn,"r");
  if (fp1 == NULL) return 0;

  for (i = 0; i < N; i++) {
    if (1 != fscanf(fp1,"%i ",&j)) break;
    if (3 != fscanf(fp1,"%lf %lf %lf",xr + j,yr + j,zr + j)) break;
    //    printf("<read_native> %i %lf %lf %lf\n",j,xr[j],yr[j],zr[j]);
    if (j < 0 || j > N-1) continue;
    nat[j] = 1;
    ++n;
  }    
  fclose(fp1);
  
  return n;
}
/****************************************************************************/
int read_native2(char *fn,double *xr,double *yr,double *zr, int a1, int a2) {
  int i, n = 0, nread;
  FILE *fp1;

  fp1 = fopen(fn,"r");
  if (fp1 == NULL) return 0;

  for (i = 0; i <= a2 ; i++) {
    if (i < a1 || i > N-1) continue;
    nread = fscanf(fp1,"%lf %lf %lf",xr + i,yr + i,zr + i);
    if (nread != 3) break;
    ++n;      
  }    
  fclose(fp1);
  
  return n;
}
/****************************************************************************/
int read_contacts(char *fn,int *ip1,int *ip2) {
  int n = 0;
  FILE *fp1;
  
  fp1 = fopen(fn,"r");

  if (fp1 == NULL) return 0;
  
  while (2 == fscanf(fp1,"%i %i",ip1 + n,ip2 + n)) {
    if (ip1[n] >= 0 && ip1[n] <= N-1 &&
	ip2[n] >= 0 && ip2[n] <= N-1)
      n++;
  }
  fclose(fp1);

  return n;
}  
/****************************************************************************/
void dumppdb(char *fn,double *o,int nobs) {
  int i,j,n;
  char str[100];
  FILE *fp;
    
  strcpy(str,OUTDIR);
  strcat(str,fn);
  
  fp = fopen(str,"w");
  
  if (fp != NULL) {

    for (i=1; i<nobs; i++)
      fprintf(fp,"# REMARK %d %lf\n",i,o[i]);
    
    for (j=0; j<NCH; j++) {
      if (CHAIN_TO_BOX) ch2box(j);
      for (i=iBeg[j]; i<=iEnd[j]; i++) 
	fprintf(fp,"ATOM  %5u  CA  ALA %c%4u    %8.3f%8.3f%8.3f  1.00\n",
		i+1,'A'+j%26,i+1,x[i],y[i],z[i]);
    }
    for (i=n=0;i<NCR;i++) {
      if (CHAIN_TO_BOX) cr2box(i);
      fprintf(fp,"HETATM %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
	      ++n," O ","HOH",i+1,xcr[i],ycr[i],zcr[i]);
    }
  }
  
  fclose(fp);
}
/****************************************************************************/
/***** Random numbers *******************************************************/
/****************************************************************************/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
double ran3n(long *idum)
/*******************************************************/
/*  returns uniform deviate r, 0<r<1                   */
/*    NR sec 7.1,  changed 1) float->double 2) r>0     */
/*******************************************************/
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
        double ret_val;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
        else {(*idum)++;}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
        ret_val = mj*FAC;
        if (mj == 0) ret_val = FAC;
	return ret_val;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/* (C) Copr. 1994 Feb 4th Jan 2:27pm  Potti-Soft  registered Trademark */
/******************************************************************/
double gasdev2(void) {
  double ran3n(long *seed);
  static int iset=0;
  static double gcos;
  double tmp1,tmp2;

  if (iset==0) { 
    tmp1=sqrt(-2*log(ran3n(&seed)));
    tmp2=pi2*ran3n(&seed);
    gcos=tmp1*cos(tmp2);
    iset=1;
    return tmp1*sin(tmp2);
  }else{
    iset=0;
    return gcos;
  }
}
/****************************************************************************/
double center_of_mass(int i1,int i2,
		      double x[],double y[],double z[]);
double correlation(int i1,int i2,
		   double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[]);
/****************************************************************************/
/* rmsd calculation                                                         */
/****************************************************************************/
double rmsd_calc(double *x1,double *y1,double *z1,
		 double *x2,double *y2,double *z2,
		 int i1,int i2) {
  int i;
  double xref1[N],yref1[N],zref1[N],g1 = 0;
  double xref2[N],yref2[N],zref2[N],g2 = 0;
  double corr;

  if (i1 < 0 || i1 > N-1) return 0;
  if (i2 < 0 || i2 > N-1) return 0;
  
  for (i = i1; i <= i2; i++) {
      xref1[i] = x1[i];
      yref1[i] = y1[i];
      zref1[i] = z1[i];
    }

  g1 = center_of_mass(i1,i2,xref1,yref1,zref1);

    for (i = i1; i <= i2; i++) {
      xref2[i] = x2[i];
      yref2[i] = y2[i];
      zref2[i] = z2[i];
    }

  g2 = center_of_mass(i1,i2,xref2,yref2,zref2);

  corr = correlation(i1,i2,xref1,yref1,zref1,xref2,yref2,zref2);
  
  return sqrt(g1 + g2 - corr);
}
/****************************************************************************/
double center_of_mass(int i1,int i2,double x[],double y[],double z[])
{
  double xcm = 0,ycm = 0,zcm = 0,gyr2 = 0;
  int j,n;

  n = 0;
  for (j = i1; j <= i2; j++) {
    xcm += x[j]; 
    ycm += y[j]; 
    zcm += z[j];
    ++n;
  }
  xcm /= n; ycm /= n; zcm /= n;
  
  for (j = i1; j <= i2; j++){
    x[j] -= xcm; 
    y[j] -= ycm; 
    z[j] -= zcm;
    gyr2 += x[j] * x[j] + y[j] * y[j] + z[j] * z[j];
  }

  return gyr2/n;
}
/****************************************************************************/
double correlation(int i1,int i2,
		   double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[])
{
  double R[3][3],RtR[3][3];
  double a,b,c,q,r,theta,detR,w[3];
  double pi2=2*acos(-1.);
  int i,j,im,n;
 
  for (i=0;i<3;i++) for (j=0;j<3;j++) R[i][j]=0;
  n=0;
  for (i=i1;i<=i2;i++) {
    R[0][0]+=x1[i]*x2[i]; R[0][1]+=x1[i]*y2[i]; R[0][2]+=x1[i]*z2[i];
    R[1][0]+=y1[i]*x2[i]; R[1][1]+=y1[i]*y2[i]; R[1][2]+=y1[i]*z2[i];
    R[2][0]+=z1[i]*x2[i]; R[2][1]+=z1[i]*y2[i]; R[2][2]+=z1[i]*z2[i];
    ++n;
  }
  detR=+R[0][0]*R[1][1]*R[2][2]
       +R[0][1]*R[1][2]*R[2][0]
       +R[0][2]*R[1][0]*R[2][1]
       -R[2][0]*R[1][1]*R[0][2]
       -R[2][1]*R[1][2]*R[0][0]
       -R[2][2]*R[1][0]*R[0][1];
  if (detR==0) {printf("****** WARNING! *******\ndetR=0\n");}

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      RtR[i][j]=R[0][i]*R[0][j]+R[1][i]*R[1][j]+R[2][i]*R[2][j];
    }
  }

  a=-(RtR[0][0]+RtR[1][1]+RtR[2][2]);
  b=RtR[0][0]*RtR[1][1]+RtR[0][0]*RtR[2][2]+RtR[1][1]*RtR[2][2]-
    RtR[0][1]*RtR[0][1]-RtR[0][2]*RtR[0][2]-RtR[1][2]*RtR[1][2];
  c=-detR*detR;

  q=(a*a-3*b)/9;
  r=(2*a*a*a-9*a*b+27*c)/54;
  if (r*r/(q*q*q)>=1) {
    printf("error r^2=%e q^3=%e %e\n",r*r,q*q*q,r*r/(q*q*q));
    exit(-1);
  }
  q=sqrt(q);
  theta=acos(r/(q*q*q));
  w[0]=sqrt(fabs(-2*q*cos(theta/3)-a/3));
  w[1]=sqrt(fabs(-2*q*cos((theta+pi2)/3)-a/3));
  w[2]=sqrt(fabs(-2*q*cos((theta-pi2)/3)-a/3));
  
  if (detR<0) {
    im = w[1] > w[0] ? 1 : 0;
    return 2.*(w[im]+fabs(w[(im+1)%3]-w[(im+2)%3]))/n;
  }else
    return 2.*(w[0]+w[1]+w[2])/n;
}
/****************************************************************************/
# undef NMAX
