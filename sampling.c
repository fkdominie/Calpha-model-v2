# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "sys.h"
# include "defs.h"
# include "global.h"
/************* Langevin dynamics parameters *********************************/
long imd,imd0;                     /* number of time steps                  */
const double tau=4.0;              /* time scale sqrt(m*l*l/eps0)           */
double dt;                         /* time step                             */
/* Protein: */ 
double vx[N],vy[N],vz[N];          /* velocities                            */
double gam;                        /* friction coefficient                  */
double mbd;                        /* mass                                  */
double tconstbd[NTMP];             /* random force                          */
double c1,c2,c3;                   /* mdstep parameters                     */
/* Crowders: */
double vxc[NCR],vyc[NCR],vzc[NCR]; /* velocities                            */
double gamcr;                      /* friction coefficient */
double mcr;                        /* mass */
double tconstcr[NTMP];             /* random force */
double c1cr, c2cr, c3cr;           /* mdstep parameters  */
/****************************************************************************/
double g[NTMP];                    /* g paramters                           */
double beta[NTMP];                 /* inverse temperatures                  */
int ind;                           /* temperature index                     */
long int nflp=0,accflp=0;          /* temperature flips                     */
/****************************************************************************/
void tflip(double e) {
  int nind = ind;
  
  nflp++;
  
  nind += (ran3n(&seed) < .5 ? -1 : 1);
  
  if (nind < 0  || nind > NTMP-1)
    return ;

  if (log(ran3n(&seed)) < - (beta[nind] - beta[ind]) * e + g[nind] - g[ind]) {
    ind = nind;
    accflp++;
  } 
      
  return;
}
/****************************************************************************/
/***** Langevin dynamics ****************************************************/
/****************************************************************************/
void mdstep() {
  int i;

  for (i=0;i<N;i++) {
    x[i]=x[i]+dt*vx[i]+c1*(fx[i]-mbd*gam*vx[i]+frdx[i]);
    y[i]=y[i]+dt*vy[i]+c1*(fy[i]-mbd*gam*vy[i]+frdy[i]);
    z[i]=z[i]+dt*vz[i]+c1*(fz[i]-mbd*gam*vz[i]+frdz[i]);
  }

  for (i= 0; i < NCR; i++){
    xcr[i]=xcr[i]+dt*vxc[i]+c1cr*(fxc[i]-mcr*gamcr*vxc[i]+frcdx[i]);
    ycr[i]=ycr[i]+dt*vyc[i]+c1cr*(fyc[i]-mcr*gamcr*vyc[i]+frcdy[i]);
    zcr[i]=zcr[i]+dt*vzc[i]+c1cr*(fzc[i]-mcr*gamcr*vzc[i]+frcdz[i]);
  }
  
  if (1!=cart2dof(1)) carterr++;

  for (i=0;i<N;i++) {
    frdxo[i]=frdx[i];
    frdyo[i]=frdy[i];
    frdzo[i]=frdz[i];
    fxo[i]=fx[i];
    fyo[i]=fy[i];
    fzo[i]=fz[i];
    frdx[i]=gasdev2()*tconstbd[ind];
    frdy[i]=gasdev2()*tconstbd[ind];
    frdz[i]=gasdev2()*tconstbd[ind];
  }
  
  for (i = 0; i < NCR; i++){
    frcdxo[i]=frcdx[i];
    frcdyo[i]=frcdy[i];
    frcdzo[i]=frcdz[i];
    fxco[i]=fxc[i];
    fyco[i]=fyc[i];
    fzco[i]=fzc[i];
    frcdx[i]=gasdev2()*tconstcr[ind];
    frcdy[i]=gasdev2()*tconstcr[ind];
    frcdz[i]=gasdev2()*tconstcr[ind];
  }  
  for (i = 0; i < NCR; i++) fxc[i] = fyc[i] = fzc[i] = 0;
  for (i = 0; i < N; i++) fx[i] = fy[i] = fz[i] = 0;
  
  Epot=(Ebon=bond(0))+(Eben=bend(0))+(Erep=exvol(0))+(Etor=torsion(0))+
    (Econ=cont(0))+(Ehp=hp(0))+(Ecc=crowd_crowd(0))+(Ecb=crowd_bead(0)); 
  
  Ekin=0;
  
  for (i=0;i<N;i++) {
    vx[i]=c2*vx[i]+c3*(fxo[i]+fx[i]+frdxo[i]+frdx[i]);
    vy[i]=c2*vy[i]+c3*(fyo[i]+fy[i]+frdyo[i]+frdy[i]);
    vz[i]=c2*vz[i]+c3*(fzo[i]+fz[i]+frdzo[i]+frdz[i]);
    Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  for (i = 0; i < NCR; i++) {
    vxc[i]=c2cr*vxc[i]+c3cr*(fxco[i]+fxc[i]+frcdxo[i]+frcdx[i]);
    vyc[i]=c2cr*vyc[i]+c3cr*(fyco[i]+fyc[i]+frcdyo[i]+frcdy[i]);
    vzc[i]=c2cr*vzc[i]+c3cr*(fzco[i]+fzc[i]+frcdzo[i]+frcdz[i]);
    Ekin += vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i];
  }
  
  Ekin*=0.5;

  return ;
}
/****************************************************************************/
