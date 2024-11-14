# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "param.h"
# include "sys.h"
# include "defs.h"
# include "global.h"
/*************** energy and force terms *************************************/
double Ekin,Epot,Eben,Ebon,Erep;   /* energy terms                          */ 
double Etor,Econ,Ehp, Ecc, Ecb;    /* energy terms                          */             
double Econ1,Econ2,Ecorr;          /* energy terms                          */             
double fx[N],fy[N],fz[N];          /* conformational force                  */
double fxo[N],fyo[N],fzo[N];       /* conformational force old              */
double frdx[N],frdy[N],frdz[N];    /* random force                          */
double frdxo[N],frdyo[N],frdzo[N]; /* random force old                      */
double fxc[NCR],fyc[NCR],fzc[NCR];          /* force crowders               */
double fxco[NCR],fyco[NCR],fzco[NCR];       /* force crowders, old          */
double frcdx[NCR],frcdy[NCR],frcdz[NCR];    /* random force crowders        */
double frcdxo[NCR],frcdyo[NCR],frcdzo[NCR]; /* random force crowders, old   */
/************* native structure *********************************************/
double xnat[N],ynat[N],znat[N];    /* native structure 1                    */
double xnat2[N],ynat2[N],znat2[N]; /* native structure 2                    */
double xnat3[N],ynat3[N],znat3[N]; /* native structure 3                    */
double xnat4[N],ynat4[N],znat4[N]; /* native structure 4                    */
double xnat5[N],ynat5[N],znat5[N]; /* native structure 4                    */
double bn[N],thn[N],phn[N];        /* bond lengths, bond & torsion angles   */
double bn2[N],thn2[N],phn2[N];     /* bond lengths, bond & torsion angles   */
int nat[N],nat2[N],nat3[N],nat4[N],nat5[N]; /* 1 if native coordinate read  */
int tor[N],tor2[N];                /* native coordinate read                */
int npair;                         /* # native contacts                     */
int npair2;                        /* # native contacts                     */
int npair3;
int npair4;
int npair5;
int spair;                         /* # native contacts                     */
int qpair;                         /* # native contacts                     */
int ip1[MAXP],ip2[MAXP];           /* list of contacts                      */
int ip3[MAXP],ip4[MAXP];           /* list of contacts                      */
int ip5[MAXP],ip6[MAXP];           /* list of contacts                      */
int ip7[MAXP],ip8[MAXP];           /* list of contacts                      */
int ip9[MAXP],ip10[MAXP];           /* list of contacts                     */
int iq1[MAXP],iq2[MAXP];           /* list of contacts                      */
int mc1[MAXP],mc2[MAXP];           /* common native1/native2 contacts       */
int nni1[MAXP],nnj1[MAXP];
int nni2[MAXP],nnj2[MAXP];
int nni3[MAXP],nnj3[MAXP];
int nni4[MAXP],nnj4[MAXP];
int dual1[MAXP],dual2[MAXP];       /* 1 if shared contact, 0 otherwise      */
double kcon_nat1[MAXP];            /* contact strengths                     */
double kcon_nat2[MAXP];            /* contact strengths                     */
double kbond[N],kbond2[N];         /* bond strengths                        */ 
double kbend[N],kbend2[N];         /* bend strengths                        */ 
double ktor1[N],ktor3[N];          /* tors strengths                        */ 
double ktor1_2[N],ktor3_2[N];      /* tors strengths                        */ 
double distp[MAXP];                /* distances                             */
double distp2[MAXP];               /* distances                             */
double distp3[MAXP];               /* distances                             */
double distp4[MAXP];               /* distances                             */
double distp5[MAXP];               /* distances                             */
double distp6[MAXP];               /* distances                             */
double distp7[MAXP];               /* distances                             */
double distp8[MAXP];               /* distances                             */
double dist_rep1[MAXP];           /* distances                             */
double dist_rep2[MAXP];           /* distances                             */
double distg1[MAXP];              /* distances                             */
double distg2[MAXP];              /* distances                             */
double distg3[MAXP];              /* distances                             */
double distg4[MAXP];              /* distances                             */
double iqkap[MAXP];
short cc[N][N];                    /* 1 native contact, 0 otherwise         */
/************* sequence effects *********************************************/
int seq[N];                        /* sequence                              */
int seqhp[N];                      /* sequence hp: >0 hp, 0 polar           */
double kap[N];                     /* hydrophobicity strengths              */
/************** extra *******************************************************/
double kcon60,krep12,kbon2,kth2;
double cthmin,sthmin;
/****************************************************************************/
double cont_df(double s1, double s2,double d);
double cellcalc(int inc,int nl,short list[]);
void log_sum_merge(double e1, double e2, double *e,
		 double f1, double f2, double *f,
		   double bet,double c);
void crowd_ecalc(double *e, double *f,double r2,double sigma,double rho);
void eseq_calc(int i,int j,double r2,double *es,double *fs);
void bond_ecalc(double *e,double *f,double db,double kbond);
void bend_ecalc(double *e,double *f,double dth,double kbend);
void tors_ecalc(double *e,double *f,double dph,double k1,double k3);
void cont_rep_ecalc(int iflag,double *e,double *f,double r2,double sig2,double kcont);
void cont_att_ecalc(int iflag,double kstr,double r,
		    double *eg,double *eg1,double *eg2,
		    double *fg,double *fg1,double *fg2,
		    int im1,int jm1,int im2,int jm2, 
		    double dg,double dg1,double dg2,
		    double *rg1x,double *rg1y,double *rg1z,
		    double *rg2x,double *rg2y,double *rg2z);
void cont_gcalc(int iflag,double *e, double *f,double r, double r0,double kcont,double ksi);
void cont_ecalc(int iflag,double *e,double *f,double sig2,double kcont,double r2);

double cont2(int iflag);
double cont_corr(int iflag);
void add_f(int i,int j,double fr,double rx,double ry,double rz);
/****************************************************************************/
/***** ENERGIES & FORCES ****************************************************/
/****************************************************************************/
void add_f(int i,int j,double fr,double rx,double ry,double rz) {
  if (i >= 0 && i <= N-1 && j >= 0 && j <= N-1) {
    fx[i] += fr * rx;
    fy[i] += fr * ry;
    fz[i] += fr * rz;      
    fx[j] -= fr * rx;
    fy[j] -= fr * ry;
    fz[j] -= fr * rz;
  }
}
/****************************************************************************/
void log_sum_merge(double e1, double e2, double *e,
		   double f1, double f2, double *f,
		   double bet, double c) {
  double emax,fmax,wmin,wsum,del;

  if (e1 < e2) {
    (*e) = e1;
    (*f) = f1;
    emax = e2;
    fmax = f2;
  } else {
    (*e) = e2;
    (*f) = f2;
    emax = e1;
    fmax = f1;
  }

  if ((del = bet * (emax - (*e))) > 100)
    return ;
  
  wsum = 1 + (wmin = exp(- del));
  
  (*e) = (*e) - log(wsum) / bet + c;
  (*f) = ( (*f) + fmax * wmin ) / wsum;

  return;
}
/****************************************************************************/
void crowd_ecalc(double *e, double *f,double r2,double sigma,double rho) {
  double r = sqrt(r2),rcalc,rcalc12;

  rcalc = sigma / (r - rho  + sigma);
  rcalc12 = rcalc * rcalc * rcalc; //rcalc^3
  rcalc12 = rcalc12 * rcalc12;     //rcalc^6
  rcalc12 = rcalc12 * rcalc12;     //rcalc^12

  (*e) = epsilonrep * rcalc12;
  (*f) = -12 * (*e) * rcalc / sigma / r;
}
/****************************************************************************/
double crowd_crowd(int iflag){
  int i, j;
  double d, dx, dy, dz, r2;                     /* To clculate distance between crowders */                      
  double e = 0,fcc = 0,ecc = 0;                 /* To calculate force and energy between crowders */
  static double sigma,rho,rcut,rcut2,Dlim,Dlim2;/* To apply limitations and ecalsh penalty */
  FILE *fp;
  
  
   if (NCR == 0) return 0;

  if (iflag < 0) {
    rho = rcrowd + rcrowd; // range crowder-crowder interaction
    sigma = srefcr + srefcr; // softness crowder-crowder interaction 
    rcut = rho + sigma;
    rcut2 = rcut * rcut;
    Dlim = rho - sigma * (1 - pow(eclash/epsilonrep,-1./12));
    Dlim2 = Dlim * Dlim;

    printf("<crowd_crowd> \n");
    printf("  rcrowd %lf  \n", rcrowd);
    printf("  eclash %lf  Dlim %lf \n", eclash, Dlim);
    printf("  rcut %lf  \n", rcut);

    return 0;
  }
  
  if (iflag == 0) {
    
    for (i = 0; i < NCR; i++){
      for (j = 0; j < i; j++){
	dx = xcr[i] - xcr[j]; bc(&dx);
	dy = ycr[i] - ycr[j]; bc(&dy);
	dz = zcr[i] - zcr[j]; bc(&dz);
	r2 = dx * dx + dy * dy + dz *  dz;
		
	if (r2 > rcut2) continue;

	if (r2 < Dlim2) {
	  ecc = eclash * (1 + (Dlim2 - r2) / Dlim2);
	  fcc = - 2 * eclash / Dlim2;
	} else {
	  crowd_ecalc(&ecc,&fcc,r2,sigma,rho);
	}
	
	fxc[i] -= fcc * dx;
	fyc[i] -= fcc * dy;
	fzc[i] -= fcc * dz;
	
	fxc[j] += fcc * dx;
	fyc[j] += fcc * dy;
	fzc[j] += fcc * dz;
	
	e += ecc;
      }
    }
    
    return e;    
  }

  if (iflag > 0) {
    fp = fopen("results/param_check/crowdcrowd_energy.plot","w");
    
    for (d = 0.0; d < 40.0; d += 0.1) {
      r2 = d * d;
      
      if (r2 > rcut2) continue;

      if (r2 < Dlim2) {
	ecc = eclash * (1 + (Dlim2 - r2) / Dlim2);
	fcc = - 2 * eclash / Dlim2;
      } else {
	crowd_ecalc(&ecc,&fcc,r2,sigma,rho);
      }  
      
      fprintf(fp,"%lf  %lf %lf %lf\n",d,ecc,fcc,fcc*d);
    }
    
    fclose(fp);
    
    return 0;
  }
  
  return 0;
}   
/****************************************************************************/
double crowd_bead(int iflag){
  int i, j,k;
  double d, dx, dy, dz, r2;                          /* To calculate distance between crowders */    
  double ecb = 0, e=0, fcb = 0;                      /* To calculate force and energy between crowders */
  static double rcut, rcut2, Dlim, Dlim2;            /* To apply limitations and ecalsh penalty */
  static double rho,sigma;
  FILE *fp;
  
  if (NCH == 0 || NCR == 0) return 0;
  
  if (iflag < 0) {
    rho = rcrowd + sigsa;    // range crowder-bead interaction
    sigma = srefcr + sigsa;  // softness crowder-bead interaction 
    rcut = rho + sigma;
    rcut2 = rcut * rcut;
    Dlim = rho - sigma * (1 - pow(eclash/epsilonrep,-1./12));
    Dlim2 = Dlim * Dlim;

    printf("<crowd_bead> \n");
    printf("  rcrowd %lf  \n", rcrowd);
    printf("  eclash %lf Dlim %lf \n", eclash, Dlim);
    printf("  rcut %lf  \n", rcut);

    return 0;    
  }

  if (iflag == 0) {

    for (i = 0; i < NCR; i++) {
      for (j = 0; j < NCH; j++) {
	for (k = iBeg[j]; k <= iEnd[j]; k++){ 
	  dx = xcr[i] - x[k]; bc(&dx);
	  dy = ycr[i] - y[k]; bc(&dy);
	  dz = zcr[i] - z[k]; bc(&dz);
	  r2 = dx * dx + dy * dy + dz *  dz;
	  
	  if (r2 > rcut2) continue;
	    
	  if (r2 < Dlim2) {
	    ecb = eclash * (1 + (Dlim2 - r2) / Dlim2);
	    fcb = -2 * eclash / Dlim2;
	  } else {
	    crowd_ecalc(&ecb,&fcb,r2,sigma,rho);
	  } 

	  fxc[i] -= fcb * dx;
	  fyc[i] -= fcb * dy;
	  fzc[i] -= fcb * dz;
	  fx[k] += fcb * dx;
	  fy[k] += fcb * dy;
	  fz[k] += fcb * dz;

	  e += ecb;
	}
      }
    }

    return e;    
  }


  if (iflag > 0) {
    fp = fopen("results/param_check/crowdbead_energy.plot","w");
    
    for (d = 0.0; d < 40.0; d += 0.1) {
      r2 = d * d;

      if (r2 > rcut2) continue;

      if (r2 < Dlim2) {
	ecb = eclash * (1 + (Dlim2 - r2) / Dlim2);
	fcb = -2 * eclash / Dlim2;
      } else {
	crowd_ecalc(&ecb,&fcb,r2,sigma,rho);
      }
      
      fprintf(fp,"%lf  %lf %lf %lf\n",d,ecb,fcb,fcb*d);
    }
    
    fclose(fp);
    
    return 0;
  }

  return 0;
}   
/****************************************************************************/
void bond_ecalc(double *e,double *f,double db,double kbon) {
  
  (*e) = kbon * db * db;
  (*f) = - kbon * 2 * db;

  return ;
}
/****************************************************************************/
double bond(int iflag) {
  int i,j,k;
  double fbx,fby,fbz,fb,fb1 = 0,fb2 = 0;
  double e = 0,e1 = 0,e2 = 0,et,db,db2,bb;
  static double bet=10.0;
  FILE *fp;

  if (FF_BOND == 0) return 0;

  if (iflag < 0) {
    bet=10.0;
    printf("<bond> bet %lf\n",bet);
    return 0;
  }

  if (iflag == 0) {
    for (k = 0; k < NCH; ++k) {
      for (i = iBeg[k]; i < iEnd[k]; ++i) {
	j = i + 1;
	//	printf("%i %lf %lf %lf\n",i,b[i],bn[i],kbond[i]);
	db = b[i] - bn[i];	
	bond_ecalc(&et,&fb,db,kbond[i]);
	
	if (FF_BOND == 2 && kbond2[i] != 0) {
	  db2 = b[i] - bn2[i];
	  bond_ecalc(&e2,&fb2,db2,kbond2[i]);
	  log_sum_merge(e1=et,e2,&et,fb1=fb,fb2,&fb,bet,0);
	}

	e += et;
	
	fbx = fb * bx[i];
	fby = fb * by[i];
	fbz = fb * bz[i];
	
	fx[i] -= fbx;
	fy[i] -= fby;
	fz[i] -= fbz;
	fx[j] += fbx;
	fy[j] += fby;
	fz[j] += fbz;
      }
    }

    return e;
  }

  if (iflag > 0) {
    
    fp = fopen("results/param_check/bond_energy.plot","w");
    for (k = 0; k < NCH; ++k) {
      for (i = iBeg[k]; i < iEnd[k]; ++i) {
	j = i + 1;
    
	for (bb=3.0; bb<5.0; bb+=0.01) {
	  db = bb - bn[i];
	  bond_ecalc(&et,&fb,db,kbond[i]);
	  e2 = e1 = et;
	  fb2 = fb1 = fb;
	  
	  if (FF_BOND == 2 && kbond2[i] != 0) {
	    db2 = bb - bn2[i];
	    bond_ecalc(&e2,&fb2,db2,kbond2[i]);
	    log_sum_merge(e1,e2,&et,fb1,fb2,&fb,bet,0);
	  }

	  fprintf(fp,"%i %lf  %lf %lf %lf  %lf %lf %lf\n",i,bb,e1,e2,et,fb1,fb2,fb);
	}
      }
    }
    fclose(fp);
    
    return 0;
  }

  return 0;
}
/****************************************************************************/
void bend_ecalc(double *e,double *f,double dth,double kth) {

  (*e) = kth * dth * dth;
  (*f) = - kth * 2 * dth;

  return ;
}
/****************************************************************************/
double bend(int iflag) {
  int i,j,k,l;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z,b2;
  double dix,diy,diz;
  double dkx,dky,dkz;
  double cth,sth,dth,dth2;
  double e = 0,e1,e2,et,fben,fben1,fben2,d;
  static double bet=5.0;
  FILE *fp;

  if (FF_BEND == 0) return 0;

  if (iflag < 0) {
    bet = 5.0;
    printf("<bend> bet %lf\n",bet);
    return 0;
  }
  
  if (iflag == 0) {
    for (l=0; l<NCH; ++l) {
      for (i=iBeg[l]; i<iEnd[l]-1; ++i) {
	j = i + 1;
	k = i + 2;
	dth = th[j] - thn[j];
	bend_ecalc(&et,&fben,dth,kbend[j]);
	
	if (FF_BEND == 2  && kbend2[j] != 0) {
	  dth2 = th[j]-thn2[j];
	  bend_ecalc(&e2,&fben2,dth2,kbend2[j]);
	  log_sum_merge(e1=et,e2,&et,fben1=fben,fben2,&fben,bet,0);
	}
	
	e += et;
	
	cth=max(cos(th[j]),-cthmin);
	sth=max(sin(th[j]),sthmin);
	if (sin(th[j])<sthmin) therr++;
	
	b1x=bx[i];
	b1y=by[i];
	b1z=bz[i];
	b1=b[i];
	b2x=bx[j];
	b2y=by[j];
	b2z=bz[j];
	b2=b[j];
	
	dix=-(b2x+cth*b1x)/sth/b1;
	diy=-(b2y+cth*b1y)/sth/b1;
	diz=-(b2z+cth*b1z)/sth/b1;
	dkx=(b1x+cth*b2x)/sth/b2;
	dky=(b1y+cth*b2y)/sth/b2;
	dkz=(b1z+cth*b2z)/sth/b2;
	
	fx[i]+=fben*dix;
	fy[i]+=fben*diy;
	fz[i]+=fben*diz;
	fx[j]+=fben*(-dix-dkx);
	fy[j]+=fben*(-diy-dky);
	fz[j]+=fben*(-diz-dkz);
	fx[k]+=fben*dkx;
	fy[k]+=fben*dky;
	fz[k]+=fben*dkz;
      }
    }

    return e;
  }

  if (iflag > 0)  {
    fp = fopen("results/param_check/bend_energy.plot","w");
    
    for (l=0; l<NCH; ++l) {
      for (i=iBeg[l]; i<iEnd[l]-1; ++i) {
	j=i+1;
	k=i+2;

	for (d = 0; d < pi; d += pi/180) {
	  dth = d - thn[j];
	  bend_ecalc(&et,&fben,dth,kbend[j]);
	  e2 = e1 = et;
	  fben2 = fben1 = fben;
	  
	  if (FF_BEND == 2  && kbend2[j] != 0) {
	    dth2 = d - thn2[j];
	    bend_ecalc(&e2,&fben2,dth2,kbend2[j]);
	    log_sum_merge(e1=et,e2,&et,fben1=fben,fben2,&fben,bet,0);
	  }

	  fprintf(fp,"%i %lf  %lf %lf %lf  %lf %lf %lf\n",j,d*180/pi,e1,e2,et,
		  fben1,fben2,fben);
	}
      }
    }
    
    return 0;
  }
  
  return 0;
}
/****************************************************************************/
void tors_ecalc(double *e,double *f,double dph,double k1,double k3) {
  double dph3 = 3 * dph;
  
  (*e) = k1 * (1 - cos(dph)) + k3 * (1 - cos(dph3));
  (*f) = - k1 * sin(dph) - 3 * k3 * sin(dph3); 

  return ;
}
/****************************************************************************/
double torsion(int iflag) {
  int i,j,k,l,m;
  double e=0,e1=0,e2=0,et=0;
  double fph=0,fph1=0,fph2=0;
  double b1,b2,b3;
  double dix,diy,diz;
  double djx,djy,djz;
  double dkx,dky,dkz;
  double dlx,dly,dlz;
  double cth1,cth2,sth1,sth2,dph,dph2;
  static double bet=5.0;
  
  if (FF_TORS == 0) return 0;

  if (iflag < 0) {
    bet = 5.0;
    printf("<torsion> bet %lf\n",bet);
    return 0;
  }
  
  if (iflag == 0) {
    for (m = 0; m < NCH; ++m) {
      for (i = iBeg[m]; i < iEnd[m] - 2; ++i) {
	j = i + 1;
	k = i + 2;
	l = i + 3;
	dph = ph[j] - phn[j];
	
	if (FF_TORS == 1) {
	  tors_ecalc(&et,&fph,dph,ktor1[j],ktor3[j]);
	}
	
	if (FF_TORS == 2) {
	  dph2 = ph[j] - phn2[j];
	  if (tor[j]) tors_ecalc(&e1,&fph1,dph,ktor1[j],ktor3[j]);
	  if (tor2[j]) tors_ecalc(&e2,&fph2,dph2,ktor1_2[j],ktor3_2[j]);
	  if (tor[j] && tor2[j]) log_sum_merge(e1,e2,&et,fph1,fph2,&fph,bet,0);
	  else {
	    et = (tor[j] ? e1 : e2);
	    fph = (tor[j] ? fph1 : fph2);
	  }
	  //	  printf("j=%i e1 %lf e2 %lf et %lf\n",j,e1,e2,et);
	}
	
	e += et;
	
      	cth1=max(cos(th[j]),-cthmin);
	sth1=max(sin(th[j]),sthmin);
	cth2=max(cos(th[k]),-cthmin);
	sth2=max(sin(th[k]),sthmin);
	
	b1=b[i]; 
	b2=b[j]; 
	b3=b[k]; 

	dlx=sx[k]/b3/sth2;
	dly=sy[k]/b3/sth2;
	dlz=sz[k]/b3/sth2;
	dix=-sx[j]/b1/sth1;
	diy=-sy[j]/b1/sth1;
	diz=-sz[j]/b1/sth1;
	djx=-b3/b2*cth2*dlx-(1-b1/b2*cth1)*dix;
	djy=-b3/b2*cth2*dly-(1-b1/b2*cth1)*diy;
	djz=-b3/b2*cth2*dlz-(1-b1/b2*cth1)*diz;
	dkx=-b1/b2*cth1*dix-(1-b3/b2*cth2)*dlx;
	dky=-b1/b2*cth1*diy-(1-b3/b2*cth2)*dly;
	dkz=-b1/b2*cth1*diz-(1-b3/b2*cth2)*dlz;
	
	fx[i]+=fph*dix;
	fy[i]+=fph*diy;
	fz[i]+=fph*diz;
	fx[j]+=fph*djx;
	fy[j]+=fph*djy;
	fz[j]+=fph*djz;
	fx[k]+=fph*dkx;
	fy[k]+=fph*dky;
	fz[k]+=fph*dkz; 
	fx[l]+=fph*dlx;
	fy[l]+=fph*dly;
	fz[l]+=fph*dlz;
      }

      return e;
    }
  }

  if (iflag > 0) { /* print energy function */
    FILE *fp;
    double d;
    
    fp = fopen("results/param_check/tors_energy.plot","w");
    for (m = 0; m < NCH; ++m) {
      for (i = iBeg[m]; i < iEnd[m] - 2; ++i) {
	j = i + 1;
	k = i + 2;
	l = i + 3;

	for (d = -pi; d < pi; d += pi/360) {
	  dph = d - phn[j];
	  
	  if (FF_TORS == 1) {
	    tors_ecalc(&et,&fph,dph,ktor1[j],ktor3[j]);
	    e1 = et; fph1 = fph;
	  }
	  
	  if (FF_TORS == 2) {
	    dph2 = d - phn2[j];
	    if (tor[j]) tors_ecalc(&e1,&fph1,dph,ktor1[j],ktor3[j]);
	    if (tor2[j]) tors_ecalc(&e2,&fph2,dph2,ktor1_2[j],ktor3_2[j]);
	    if (tor[j] && tor2[j]) log_sum_merge(e1,e2,&et,fph1,fph2,&fph,bet,0);
	    else {
	      et = (tor[j] ? e1 : e2);
	      fph = (tor[j] ? fph1 : fph2);
	    }
	  }
	  
	  fprintf(fp,"%i %lf %lf %lf %lf %lf %lf %lf \n",j,d*180/pi,e1,e2,et,
		  fph1,fph2,fph);
	}
      }
    }
    
    fclose(fp);
    
    return 0;
  }
  
  return 0;
}
/****************************************************************************/
double exvol(int iflag)
{
  /* iflag < 0 initialize                               */
  /* iflag > 0 calculate the full energy                */

  int i,j,ix,iy,iz,ic,nec=0;
  int lx[N],ly[N],lz[N],lc[N],in_cell,nl,cn;
  double e=0;
  short a,list[N],pnt[N]; 

  static short cell[MAXCELL];      /* division into cells */
  static double cutg;
  static int ns,h2,h3;

  if (FF_EXVOL == 0 || N == 0) return 0;

  if (iflag < 0) {
    for (i=0; i<MAXCELL; ++i) cell[i]=-1;

    ns = BOX / cut;
    h2 = ns*ns; 
    h3 = h2*ns;
    cutg = BOX/(double)ns;
    cellcalc(-1,0,list);

    if (ns * ns * ns > MAXCELL) {printf("# cells > MAXCELL\n"); exit(-1);}
    printf("<exvol> # cells %i ns %i cutg %f\n",ns*ns*ns,ns,cutg);

    return 0.0;
  }

  in2box();
  for (i=N-1; i>=0; --i) {
    ix=xb[i]/cutg; iy=yb[i]/cutg; iz=zb[i]/cutg;
    ic=ix+ns*(iy+ns*iz);
    pnt[i]=cell[ic];
    if (cell[ic]<0) {lx[nec]=ix; ly[nec]=iy; lz[nec]=iz; lc[nec++]=ic;}
    cell[ic]=i;
  }

  for (i=0;i<nec;i++) { 
    nl=0; 
    list[nl++]=a=cell[(j=lc[i])]; 
    while ((a=pnt[a])>=0) {list[nl++]=a;}
    in_cell=nl;
    cn=j+1; 
    if (lx[i]+1==ns) cn-=ns; 
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+ns; 
    if (ly[i]+1==ns) cn-=h2;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h2; 
    if (lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1+ns; 
    if (lx[i]+1==ns) cn-=ns;
    if (ly[i]+1==ns) cn-=h2;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1+h2; 
    if (lx[i]+1==ns) cn-=ns;
    if(lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+ns+h2;
    if (ly[i]+1==ns) cn-=h2;
    if (lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1-ns; 
    if (lx[i]+1==ns) cn-=ns;
    if (ly[i]==0) cn+=h2; 
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1-h2; 
    if (lx[i]+1==ns) cn-=ns;
    if (lz[i]==0) cn+=h3; 
    if ((a=cell[cn])>=0) {  
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+ns-h2;
    if (ly[i]+1==ns) cn-=h2;
    if (lz[i]==0) cn+=h3;
    if ((a=cell[cn])>=0) {  
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1+ns+h2;
    if (lx[i]+1==ns) cn-=ns;
    if (ly[i]+1==ns) cn-=h2;
    if (lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1+ns-h2;
    if (lx[i]+1==ns) cn-=ns;
    if (ly[i]+1==ns) cn-=h2;
    if (lz[i]==0) cn+=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+1-ns+h2;
    if (lx[i]+1==ns) cn-=ns;
    if (ly[i]==0) cn+=h2;
    if (lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j-1+ns+h2;
    if (lx[i]==0) cn+=ns;
    if (ly[i]+1==ns) cn-=h2;
    if (lz[i]+1==ns) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    e+=cellcalc(in_cell,nl,list);
  }

  for (i=0;i<nec;i++) cell[lc[i]]=-1;

  return krep*e;
}
/****************************************************************************/
double cellcalc(int inc,int nl,short list[])
{
  int i,j,n,m,ic;
  double r2,r6,r12,rx,ry,rz,fr,e=0;
  static double cut2,sigsa2;
  static double asa,bsa,bsa2;

  if (inc<0) {
    cut2=cut*cut;
    sigsa2=sigsa*sigsa;
    asa=-7*pow(sigsa2/cut2,6.0);
    bsa=6*pow(sigsa2/cut2,6.0)/cut2;
    bsa2=2*bsa;
    printf("<cellcalc> asa %e bsa %e\n",asa,bsa);
    return 0.0;
  }

  for (n=0;n<inc;n++) {
    i=list[n]; ic=a2c[i];
    for (m=n+1;m<nl;m++) {
      j=list[m];
      if ( (abs(i-j)<4 && ic == a2c[j]) || cc[i][j] == 1 || cc[j][i] == 1) continue;
      rx=x[j]-x[i]; bc(&rx);
      ry=y[j]-y[i]; bc(&ry);
      rz=z[j]-z[i]; bc(&rz);
      if ((r2=rx*rx+ry*ry+rz*rz)>cut2) continue;
      r6=sigsa2/r2; r6=r6*r6*r6; r12=r6*r6;
      e+=r12+asa+bsa*r2;
      fr=krep12*r12/r2-bsa2;

      fx[i]-=fr*rx;
      fy[i]-=fr*ry;
      fz[i]-=fr*rz;
      fx[j]+=fr*rx;
      fy[j]+=fr*ry;
      fz[j]+=fr*rz;
    }
  }

  return e;
}
/****************************************************************************/
double cont_df(double s1, double s2,double d) {
  /* Test: numerical derivation of cont() energy */
  double h = 0.001;
  double e1,e2,em1,em2;
  double f1,f2,fm;
  double df,bet=5.0;

  d -= h;
  cont_ecalc(0,&e1,&f1,s1,kcon,d * d);
  cont_ecalc(0,&e2,&f2,s2,kcon,d * d);
  log_sum_merge(e1,e2,&em1,f1,f2,&fm,bet,0);

  d += 2 * h;

  cont_ecalc(0,&e1,&f1,s1,kcon,d * d);
  cont_ecalc(0,&e2,&f2,s2,kcon,d * d);
  log_sum_merge(e1,e2,&em2,f1,f2,&fm,bet,0);

  df = (em2 - em1) / (2 * h);

  return df;
}

/****************************************************************************/
void cont_gcalc(int iflag,double *e, double *f, double r, double r0,
		double kcont,double ksi) {
  double dr,g;

  if (iflag < 0) {
    return ;
  }

  dr = r - r0;
  g = exp(- dr * dr / 2 / ksi);

  (*e) = kcont * g;
  (*f) = kcont * dr * g / ksi / r;
}
/****************************************************************************/
void cont_rep_ecalc(int iflag,double *e,double *f,double r2,double sig2,
		    double kcont) {
  double r4,r6;

  if (iflag < 0) {
    return ;
  }

  if (r2 >= sig2) {
    (*e) = (*f) = 0;
    return ;
  }
  
  r6 = sig2 / r2; 
  r4 = r6 * r6; 
  r6 = r4 * r6;
  
  (*e) = kcont * ( r6 * (r6 - 2) + 1 );
  (*f) = kcont * 12 * r6 * (r6 - 1) / r2 ;
  
  return ;
}
/****************************************************************************/
void cont_att_ecalc(int iflag,double kstr,double r,
		    double *eg,double *eg1,double *eg2,
		    double *fg,double *fg1,double *fg2,
		    int im1,int jm1,int im2,int jm2, 
		    double dg,double dg1,double dg2,
		    double *rg1x,double *rg1y,double *rg1z,
		    double *rg2x,double *rg2y,double *rg2z) {
  double rg1,rg2;
  
  if (iflag < 0) {
    return;
  }
  
  cont_gcalc(0,eg,fg,r,dg,kstr,ksi1);
	
  if (im1 >= 0 && jm1 >= 0) {
    rg1 = vec2(im1,jm1,rg1x,rg1y,rg1z);
    cont_gcalc(0,eg1,fg1,sqrt(rg1),dg1,1.0,ksi2);
  } else *eg1 = *fg1 = 1.0;

  if (im2 >= 0 && jm2 >= 0) {
    rg2 = vec2(im2,jm2,rg2x,rg2y,rg2z);
    cont_gcalc(0,eg2,fg2,sqrt(rg2),dg2,1.0,ksi2);
  } else *eg2 = *fg2 = 1.0;

  (*fg1) = - (*eg) * (*fg1) * (*eg2);
  (*fg2) = - (*eg) * (*eg1) * (*fg2);
  (*eg) = - (*eg) * (*eg1) * (*eg2);
  (*fg) = - (*fg) * (*eg1) * (*eg2);
}
/****************************************************************************/
void cont_ecalc(int iflag,double *e,double *f,double sig2,double kcont,
		double r2) {
  double r4,r6;

  if (iflag < 0) {
    return ;
  }

  r6 = sig2 / r2; 
  r4 = r6 * r6; 
  r6 = r4 * r6;
  
  (*e) = kcont * ( r6 * (5 * r6 - 6 * r4)  );
  (*f) = kcont * 60 * ( r6 * (r6 - r4) / r2  );
  
  return ;
}
/****************************************************************************/
double cont(int iflag) {
  int i,j,m,im1,jm1,im2,jm2;
  double r,r2,rx,ry,rz,sig2;
  double fr,er,eg,fg,e=0;
  double rg1x,rg1y,rg1z,eg1,fg1;
  double rg2x,rg2y,rg2z,eg2,fg2;
  double d;
  static double rcut_mb;
  FILE *fp1;

  if (FF_CONT == 0 || N == 0) return 0;
  
  if (iflag < 0) {
    rcut_mb = 3 * ksi1;
    
    if (FF_MULTIBODY)  {
      printf("<cont> FF_MULTIBODY %i\n",FF_MULTIBODY);
      printf("<cont> ksi1 %lf ksi2 %lf\n",ksi1,ksi2);
      printf("<cont> cutoff %lf \n",rcut_mb);
    }
    
    cont_ecalc(-1,&e,&fr,0,0,0);
    cont_gcalc(-1,&e,&fr,0,0,0,0);
    cont_rep_ecalc(-1,&e,&fr,0,0,0);
    
    if (FF_CONT == 2) {
      cont2(-1);
      cont_corr(-1);
    }

    fp1 = fopen("results/param_check/cont_param.out","w");
    for (m = 0; m < npair; m++) {
      if (dual1[m] > 0) continue;
      i = ip1[m]; j = ip2[m];
      fprintf(fp1,"%i %i %lf %lf %lf %i %i %i %i %lf %lf \n",i,j,
	      kcon_nat1[m],distp[m],sqrt(dist_rep1[m]),
	      nni1[m],nnj1[m],nni2[m],nnj2[m],
	      distg1[m],distg2[m]);
    }
    fclose(fp1);
    
    return 0;
  }

  if (iflag == 0) {
    for (m = 0; m < npair; m++) {
      if (dual1[m] > 0) continue;

      i=ip1[m]; j=ip2[m];
      r2 = vec2(i,j,&rx,&ry,&rz);

      if (!FF_MULTIBODY) {
	if ( r2 > 4 * (sig2 = distp2[m]) ) continue;
	cont_ecalc(0,&er,&fr,sig2,kcon_nat1[m],r2);
	e += er;
	add_f(i,j,-fr,rx,ry,rz);
      }
      
      if (FF_MULTIBODY) {
	if ( (r = sqrt(r2)) > distp[m] + rcut_mb ) continue;
	im1 = nni1[m]; jm1 = nnj1[m];
	im2 = nni2[m]; jm2 = nnj2[m];
	cont_rep_ecalc(0,&er,&fr,r2,dist_rep1[m],krep);
	cont_att_ecalc(0,kcon_nat1[m],r,
		       &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
		       distp[m],distg1[m],distg2[m],
		       &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);
	e += er + eg;
	fr = fr + fg;
	add_f(i,j,-fr,rx,ry,rz);
	add_f(im1,jm1,-fg1,rg1x,rg1y,rg1z);
	add_f(im2,jm2,-fg2,rg2x,rg2y,rg2z);
      }
    }

    Econ1 = e;
    
    if (FF_CONT == 2 && iflag == 0) {
      e += (Econ2 = cont2(0));
      e += (Ecorr = cont_corr(0));
    }
    
    return e;
  }
  
  if (iflag > 0) {
    fp1 = fopen("results/param_check/cont_energy1.plot","w");
    
    for (m = 0; m < npair; m++) {
      if (dual1[m] > 0) continue;
      i = ip1[m]; j = ip2[m];

      for (d = 3.0; d < 20; d += 0.01) {
	r2 = d * d;

	if (!FF_MULTIBODY) {
	  if ( r2 > 4 * (sig2 = distp2[m]) ) continue;
	  cont_ecalc(0,&er,&fr,sig2,kcon_nat1[m],r2);
	  fprintf(fp1,"%i %i %i %lf  %lf\n",m,i,j,d,er);
	}

	if (FF_MULTIBODY) {
	  if ( (r = sqrt(r2)) > distp[m] + rcut_mb ) continue;
	  cont_rep_ecalc(0,&er,&fr,r2,dist_rep1[m],krep);
	  cont_gcalc(0,&eg,&fg,r,distp[m],kcon_nat1[m],ksi1);
	  fprintf(fp1,"%i %i %i %lf  %lf %lf %lf \n",m,i,j,d,er-eg,er,eg);
	} 

      }
    }

    fclose(fp1);
    return 0;
  }

  return 0;
}
/****************************************************************************/
double cont2(int iflag) {
  int i,j,m,im1,jm1,im2,jm2;
  double r,r2,rx,ry,rz,sig2;
  double fr,er,eg,fg,e=0;
  double rg1x,rg1y,rg1z,eg1,fg1;
  double rg2x,rg2y,rg2z,eg2,fg2;
  double d;
  static double rcut_mb;
  FILE *fp1;
  
  if (iflag < 0) {
    rcut_mb = 3 * ksi1;
    if (FF_MULTIBODY == 1)  {
      printf("<cont2> cutoff %lf \n",rcut_mb);
    }

    fp1 = fopen("results/param_check/cont_param2.out","w");
    for (m = 0; m < npair2; m++) {
      if (dual2[m] > 0) continue;
      i = ip3[m]; j = ip4[m];
      fprintf(fp1,"%i %i %lf %lf %lf %i %i %i %i %lf %lf\n",i,j,
	      kcon_nat2[m],distp3[m],sqrt(dist_rep2[m]),
	      nni3[m],nnj3[m],nni4[m],nnj4[m],
	      distg3[m],distg4[m]);
    }
    fclose(fp1);

    return 0;
  }

  if (iflag == 0) {
    for (m = 0; m < npair2; m++) {
      if (dual2[m] > 0) continue;
      i=ip3[m]; j=ip4[m];
      r2 = vec2(i,j,&rx,&ry,&rz);

      if (!FF_MULTIBODY) {
	if ( r2 > 4 * (sig2 = distp4[m]) ) continue;
	cont_ecalc(0,&er,&fr,sig2,kcon_nat2[m],r2);
	add_f(i,j,-fr,rx,ry,rz);
	e += er;
      }

      if (FF_MULTIBODY) {
	if ( (r = sqrt(r2)) > distp3[m] + rcut_mb ) continue;
	im1 = nni3[m]; jm1 = nnj3[m];
	im2 = nni4[m]; jm2 = nnj4[m];
	cont_rep_ecalc(0,&er,&fr,r2,dist_rep2[m],krep);
	cont_att_ecalc(0,kcon_nat2[m],r,
		       &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
		       distp3[m],distg3[m],distg4[m],
		       &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);
	e += er + eg;
	fr = fr + fg;
	add_f(i,j,-fr,rx,ry,rz);
	add_f(im1,jm1,-fg1,rg1x,rg1y,rg1z);
	add_f(im2,jm2,-fg2,rg2x,rg2y,rg2z);
      }
    }

    return e;
  }
  
  if (iflag > 0) {
    fp1 = fopen("results/param_check/cont_energy2.plot","w");
    
    for (m = 0; m < npair2; m++) {
      if (dual2[m] > 0) continue;
      i = ip3[m]; j = ip4[m];
      
      for (d = 3.0; d < 20; d += 0.01) {
	r2 = d * d;
	
	if (!FF_MULTIBODY) {
	  if ( r2 > 4 * (sig2 = distp4[m]) ) continue;
	  cont_ecalc(0,&er,&fr,sig2,kcon_nat2[m],r2);
	  fprintf(fp1,"%i %i %i %lf  %lf\n",m,i,j,d,er);
	}

	if (FF_MULTIBODY) {
	  if ( (r = sqrt(r2)) > distp3[m] + rcut_mb ) continue;
	  cont_rep_ecalc(0,&er,&fr,r2,dist_rep2[m],krep);
	  cont_gcalc(0,&eg,&fg,sqrt(r2),distp3[m],kcon_nat2[m],ksi1);
	  fprintf(fp1,"%i %i %i %lf  %lf %lf %lf \n",m,i,j,d,er-eg,er,eg);
	} 
      }

    }
    
    fclose(fp1);
    return 0;
  }
  
  return 0;
}
/****************************************************************************/
double cont_corr(int iflag) {
  int s,m,n,i,j;
  double r,r2,rx,ry,rz;
  double e=0,er,fr;
  static double rcut_mb;
  FILE *fp1;

  int im1A,jm1A,im2A,jm2A;
  double egA,fgA;
  double rg1xA,rg1yA,rg1zA,eg1A,fg1A;
  double rg2xA,rg2yA,rg2zA,eg2A,fg2A;

  int im1B,jm1B,im2B,jm2B;
  double egB,fgB;
  double rg1xB,rg1yB,rg1zB,eg1B,fg1B;
  double rg2xB,rg2yB,rg2zB,eg2B,fg2B;

  if (!FF_MULTIBODY) {
    printf("    cont_corr() not implemented for the case FF_MULTIBODY %i\nExiting...\n",FF_MULTIBODY);
    exit(-1);
  }

  if (iflag < 0) {
    rcut_mb = 3 * ksi1;

    if (FF_MULTIBODY) {
      fp1 = fopen("results/param_check/cont_param_shared.out","w");
      for (s = 0; s < spair; ++s) {
	m = mc1[s];
	n = mc2[s];
	
	i = ip1[m];
	j = ip2[m];

	im1A = nni1[m]; jm1A = nnj1[m];
	im2A = nni2[m]; jm2A = nnj2[m];
	im1B = nni3[n]; jm1B = nnj3[n];
	im2B = nni4[n]; jm2B = nnj4[n];

	fprintf(fp1,"%i %i kcon %lf %lf dist %lf %lf %lf %lf %i %i %i %i %i %i %i %i distg %lf %lf %lf %lf\n",i,j,
		kcon_nat1[m],kcon_nat2[n],distp[m],distp3[n],
		sqrt(dist_rep1[m]),sqrt(dist_rep2[n]),
		im1A,jm1A,im2A,jm2A,im1B,jm1B,im2B,jm2B,
		distg1[m],distg2[m],
		distg3[n],distg4[n]);
      }
      fclose(fp1);
    }

    return 0;
  }

  if (iflag == 0) {
    for (s = 0; s < spair; ++s) {
      m = mc1[s];
      n = mc2[s];

      i = ip1[m];
      j = ip2[m];

      r2 = vec2(i,j,&rx,&ry,&rz);

      if (FF_MULTIBODY) {
	if ( (r=sqrt(r2)) > distp[m] + rcut_mb  && r > distp3[n] + rcut_mb ) continue;

	im1A = nni1[m]; jm1A = nnj1[m];
	im2A = nni2[m]; jm2A = nnj2[m];
	im1B = nni3[n]; jm1B = nnj3[n];
	im2B = nni4[n]; jm2B = nnj4[n];

	cont_rep_ecalc(0,&er,&fr,r2,dist_rep1[m],krep);

	cont_att_ecalc(0,kcon_nat1[m],r,
		       &egA,&eg1A,&eg2A,&fgA,&fg1A,&fg2A,im1A,jm1A,im2A,jm2A,
		       distp[m],distg1[m],distg2[m],
		       &rg1xA,&rg1yA,&rg1zA,&rg2xA,&rg2yA,&rg2zA);
	
	cont_att_ecalc(0,kcon_nat2[n],r,
		       &egB,&eg1B,&eg2B,&fgB,&fg1B,&fg2B,im1B,jm1B,im2B,jm2B,
		       distp3[n],distg3[n],distg4[n],
		       &rg1xB,&rg1yB,&rg1zB,&rg2xB,&rg2yB,&rg2zB);

	if (egA < egB) {
	  e += er + egA;
	  fr = fr + fgA;
	  add_f(i,j,-fr,rx,ry,rz);
	  add_f(im1A,jm1A,-fg1A,rg1xA,rg1yA,rg1zA);
	  add_f(im2A,jm2A,-fg2A,rg2xA,rg2yA,rg2zA);
	} else {
	  e += er + egB;
	  fr = fr + fgB;
	  add_f(i,j,-fr,rx,ry,rz);
	  add_f(im1B,jm1B,-fg1B,rg1xB,rg1yB,rg1zB);
	  add_f(im2B,jm2B,-fg2B,rg2xB,rg2yB,rg2zB);
	} 
      }
    }
    
    return e;
  }
  
  if (iflag > 0) {
    
    return 0;
  }
  
  return 0;
}
/****************************************************************************/
double hp(int iflag) {
  int i,j,m,seqhp[N];                 
  double r,r2,rx,ry,rz,dr,fr,edr,e=0;
  static double  aco,bco,cq2;

  // untested
  
  if (FF_SEQ == 0) return 0;

  if (iflag < 0) {
    for (i=0;i<N;i++) {
      seqhp[i] = 0;
      switch (seq[i]) {
      case 'A' : {seqhp[i] = 1.0; break;}
      case 'V' : {seqhp[i] = 1.0; break;}
      case 'L' : {seqhp[i] = 1.0; break;}
      case 'I' : {seqhp[i] = 1.0; break;}
      case 'M' : {seqhp[i] = 1.0; break;}
      case 'W' : {seqhp[i] = 1.0; break;}
      case 'F' : {seqhp[i] = 1.0; break;}
      case 'Y' : {seqhp[i] = 1.0; break;}
      default  : {seqhp[i] = 0.0; break;}
      }
    }
    
    for (i=0;i<N;i++) kap[i] = (seqhp[i] > 0 ? 1.0 : 0.0);
    
    qpair=0;
    for (i=0;i<N;i++) {
      if (seqhp[i] == 0) continue;
      for (j=0;j<i-3;j++) {
	if (seqhp[j] == 0) continue;
	if (cc[i][j]!=1) {
	  iq1[qpair]=i;
	  iq2[qpair]=j;
	  iqkap[qpair++]=kap[i]*kap[j];
	}
      }
    }
      
    cq2=(sighp+cuthp)*(sighp+cuthp);
    aco=-(1+cuthp*cuthp)*exp(-cuthp*cuthp/2);
    bco=cuthp*exp(-cuthp*cuthp/2);

    printf("<hp> non-native hp pairs: qpair %i\n",qpair);
    printf("<hp> sighp %f cuthp %f\n",sighp,cuthp);
    printf("<hp> aco %e bco %e\n",aco,bco);
    return 0;
  }
  
  for (m=0;m<qpair;m++) {
    i=iq1[m]; j=iq2[m];
    rx=x[j]-x[i]; bc(&rx);
    ry=y[j]-y[i]; bc(&ry);
    rz=z[j]-z[i]; bc(&rz);
    if ((r2=rx*rx+ry*ry+rz*rz)>cq2) continue;
    dr=(r=sqrt(r2))-sighp;
    edr=exp(-dr*dr/2);
    e-=iqkap[m]*(edr+aco+bco*dr);
    fr=iqkap[m]*khp*(-dr*edr+bco)/r;
    fx[i]-=fr*rx;
    fy[i]-=fr*ry;
    fz[i]-=fr*rz;
    fx[j]+=fr*rx;
    fy[j]+=fr*ry;
    fz[j]+=fr*rz; 
  }
  return khp*e;
}
/****************************************************************************/
