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
/***** GEOMETRY *************************************************************/
/****************************************************************************/
double x[N],y[N],z[N];                   /* atom coordinates                      */
double xb[N],yb[N],zb[N];                /* periodic boundary enforced            */
double b[N];                             /* pseudo bond lenghts (0...N-1)         */
double th[N];                            /* pseudo bond angles (1...N-2)          */
double ph[N];                            /* pseudo torsion angles (1...N-3)       */
double bx[N],by[N],bz[N];                /* pseudo bond vectors                   */
double sx[N],sy[N],sz[N];                /* auxilary vector                       */
double boxhf=BOX/2.0;                    /* half box length                       */
int iBeg[NCH],iEnd[NCH];                 /* start and end indices of chains       */
int a2c[N];                              /* monomer --> chain                     */
int carterr=0,therr=0;
/****************************************************************************/
double vec2(int i,int j,double *rx,double *ry,double *rz) {
  (*rx) = x[j] - x[i]; bc(rx);
  (*ry) = y[j] - y[i]; bc(ry);
  (*rz) = z[j] - z[i]; bc(rz);
  return  (*rx)*(*rx) + (*ry)*(*ry) + (*rz)*(*rz);
}
/****************************************************************************/
void trans(int ic,double dx,double dy,double dz) {
  int i;

  for (i=iBeg[ic]; i<=iEnd[ic]; ++i) {
    x[i] += dx;
    y[i] += dy;
    z[i] += dz;
  }

  return ;
}
/****************************************************************************/
void trans_cr(int ic,double dx,double dy,double dz) {

  xcr[ic] += dx;
  ycr[ic] += dy;
  zcr[ic] += dz;

  return ;
}
/****************************************************************************/
void bc(double *x) 
{
  while ((*x)>boxhf) (*x)-=BOX;
  while ((*x)<-boxhf) (*x)+=BOX;
}
/****************************************************************************/
void in2box(void) {
  int i,j,k=0;

  for (j=0;j<NCH;j++) {
    for (i=iBeg[j];i<=iEnd[j];i++) {
      xb[k]=x[i]; yb[k]=y[i]; zb[k]=z[i];
      while (xb[k]>=BOX) xb[k]-=BOX;
      while (xb[k]<0) xb[k]+=BOX;

      while (yb[k]>=BOX) yb[k]-=BOX;
      while (yb[k]<0) yb[k]+=BOX;
      
      while (zb[k]>=BOX) zb[k]-=BOX;
      while (zb[k]<0) zb[k]+=BOX;
      k++;
    }
  }
}
/****************************************************************************/
void ch2box(int ich) {
  int i,n,nx,ny,nz;
  double xcm,ycm,zcm;
  
  xcm=ycm=zcm=n=0;
  for (i=iBeg[ich]; i<=iEnd[ich]; ++i) {
    xcm += x[i];
    ycm += y[i];
    zcm += z[i];
    ++n;
  }

  xcm /= n;
  ycm /= n;
  zcm /= n;

  nx=ny=nz=0;
  while (xcm + nx*BOX > BOX) --nx;
  while (xcm + nx*BOX < 0) ++nx; 
  while (ycm + ny*BOX > BOX) --ny;
  while (ycm + ny*BOX < 0) ++ny; 
  while (zcm + nz*BOX > BOX) --nz;
  while (zcm + nz*BOX < 0) ++nz; 

  for (i=iBeg[ich]; i<=iEnd[ich]; ++i) {
    x[i] += nx*BOX;
    y[i] += ny*BOX;
    z[i] += nz*BOX;
  }

  return;
}
/****************************************************************************/
void cr2box(int icr) {
  
  while (xcr[icr] > BOX) xcr[icr] -= BOX;
  while (xcr[icr] < 0) xcr[icr] += BOX;

  while (ycr[icr] > BOX) ycr[icr] -= BOX;
  while (ycr[icr] < 0) ycr[icr] += BOX;

  while (zcr[icr] > BOX) zcr[icr] -= BOX;
  while (zcr[icr] < 0) zcr[icr] += BOX;

  return;
}
/****************************************************************************/
void dof2cart(int iflag) {
  int i,j,k,l;
  double vx,vy,vz;
  double sth1,sth2,cth1,cth2,cph,sph;
  double renorm;
  
  if (iflag < 0) {
    printf("<dof2cart> init\n");
    for (l=0; l<NCH; ++l) {
      k = iBeg[l];
      x[k] = BOX * ran3n(&seed);
      y[k] = BOX * ran3n(&seed);
      z[k] = BOX * ran3n(&seed);
      //      printf("%d %lf %lf %lf\n",k,x[k],y[k],z[k]);
      
      x[k+1] = x[k];
      y[k+1] = y[k];
      z[k+1] = z[k] + bn[k];
      //      printf("%d %lf %lf %lf\n",k+1,x[k+1],y[k+1],z[k+1]);
      
      x[k+2] = x[k+1];
      y[k+2] = y[k+1] + bn[k+1] * sin(thn[k+1]);
      z[k+2] = z[k+1] - bn[k+1] * cos(thn[k+1]);
      //      printf("%d %lf %lf %lf\n",k+2,x[k+2],y[k+2],z[k+2]);
    }
    return ;
  }
  
  for (l=0; l<NCH; ++l) {
    k = iBeg[l];

    bx[k] = x[k+1] - x[k];
    by[k] = y[k+1] - y[k];
    bz[k] = z[k+1] - z[k];
    renorm = sqrt(bx[k]*bx[k] + by[k]*by[k] + bz[k]*bz[k]);
    bx[k] /= renorm;
    by[k] /= renorm;
    bz[k] /= renorm;
    
    k = iBeg[l] + 1;

    bx[k] = x[k+1] - x[k];
    by[k] = y[k+1] - y[k];
    bz[k] = z[k+1] - z[k];
    renorm = sqrt(bx[k]*bx[k] + by[k]*by[k] + bz[k]*bz[k]);
    bx[k] /= renorm;
    by[k] /= renorm;
    bz[k] /= renorm;
    
    for (k=iBeg[l]; k<iEnd[l]; ++k) {
      if (k > iBeg[l] + 1) {
	j=k-1;
	i=k-2;

	sth1 = sin(th[j]);
	sth2 = sin(th[k]);
	cth1 = cos(th[j]);
	cth2 = cos(th[k]);
	sph = sin(ph[j]);
	cph = cos(ph[j]);
	
	sx[j]=(by[i]*bz[j]-bz[i]*by[j])/sth1;
	sy[j]=(bz[i]*bx[j]-bx[i]*bz[j])/sth1;
	sz[j]=(bx[i]*by[j]-by[i]*bx[j])/sth1;
	
	vx=(bx[i]+bx[j]*cth1)/sth1;
	vy=(by[i]+by[j]*cth1)/sth1;
	vz=(bz[i]+bz[j]*cth1)/sth1;
	
	bx[k]=-cth2*bx[j]+sth2*(sph*sx[j]-cph*vx);
	by[k]=-cth2*by[j]+sth2*(sph*sy[j]-cph*vy);
	bz[k]=-cth2*bz[j]+sth2*(sph*sz[j]-cph*vz);

	renorm = sqrt(bx[k]*bx[k] + by[k]*by[k] + bz[k]*bz[k]);
	bx[k] /= renorm;
	by[k] /= renorm;
	bz[k] /= renorm;
      }

      x[k+1] = x[k] + b[k] * bx[k];
      y[k+1] = y[k] + b[k] * by[k];
      z[k+1] = z[k] + b[k] * bz[k];
    }
  }
}
/****************************************************************************/
int cart2dof(int err) {
  int i,j,k,ok=1;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z;
  double ux,uy,uz,u;
  double tmp1;

  /* bond vectors */

  for (k=0; k<NCH; k++) {
    for (i=iBeg[k]; i<iEnd[k]; i++) {
      j=i+1;
      b1x=x[j]-x[i];
      b1y=y[j]-y[i];
      b1z=z[j]-z[i];
      b1=sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
      
      bx[i]=b1x/b1;
      by[i]=b1y/b1;
      bz[i]=b1z/b1;
      b[i]=b1;
      
     // if (b1>5.0 || b1<2.0) {ok=0; printf("cart2dof 1: %i b1 %f\n",i,b1);}
      if (err && (b1>5.0 || b1<2.0)) {ok=0; printf("b1 %f %d \n",b1, i);}
    }    
  }

  /* bond angles */

  for (k=0; k<NCH; k++) {
    for (i=iBeg[k]+1; i<iEnd[k]; i++) {
      b1x=bx[i-1];
      b1y=by[i-1];
      b1z=bz[i-1];

      b2x=bx[i];
      b2y=by[i];
      b2z=bz[i];
    
      tmp1 = b1x*b2x + b1y*b2y + b1z*b2z;
      // if (tmp1 < -1.0) {tmp1 = -1.0; ok=0; printf("cart2dof 2: tmp1 %f\n",tmp1);}
      // if (tmp1 >  1.0) {tmp1 =  1.0; ok=0; printf("cart2dof 3: tmp1 %f\n",tmp1);}
      
      if (tmp1 < -1.0) {tmp1 = -1.0; ok=0;}
      if (tmp1 >  1.0) {tmp1 =  1.0; ok=0; }
      th[i] = pi - acos(tmp1);
      
      ux=b1y*b2z-b1z*b2y;
      uy=b1z*b2x-b1x*b2z;
      uz=b1x*b2y-b1y*b2x;
      u=sqrt(ux*ux+uy*uy+uz*uz);
      
      sx[i]=ux/u;
      sy[i]=uy/u;
      sz[i]=uz/u;
    }
  }
 
  /* torsion angles */
  
  for (k=0; k<NCH; k++) {
    for (i=iBeg[k]+1; i<iEnd[k]-1; i++) {
      j=i+1;
      
      tmp1 = sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j];
      if (tmp1 < -1.0) {tmp1 = -1.0; ok=0;}
      if (tmp1 >  1.0) {tmp1 =  1.0; ok=0;}
      // if (tmp1 < -1.0) {tmp1 = -1.0; ok=0; printf("cart2dof 5: tmp1 %f\n",tmp1);}
      // if (tmp1 >  1.0) {tmp1 =  1.0; ok=0; printf("cart2dof 6: tmp1 %f\n",tmp1);}
      
      if ((sx[i]*bx[j]+sy[i]*by[j]+sz[i]*bz[j])>0.0)
	ph[i]=acos(tmp1);
      else
	ph[i]=-acos(tmp1);
    }
  }

  return ok;
}
/****************************************************************************/
