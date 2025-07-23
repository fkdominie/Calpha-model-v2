# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
/******************/
# include "param.h"
# include "sys.h"
# include "defs.h"
# include "global.h"
/*************** energy and force terms *************************************/
double Ekin,Epot,Eben,Ebon,Erep;   /* energy terms                          */ 
double Etor,Econ,Ecc,Ecb;          /* energy terms                          */             
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
double xnat[N],ynat[N],znat[N];    /* native 1 structure                    */
double xnat2[N],ynat2[N],znat2[N]; /* native 2 structure                    */
double bn[N],thn[N],phn[N];        /* native 1 reference values             */
double bn2[N],thn2[N],phn2[N];     /* native 2 reference values             */
int npair;                         /* # native contacts                     */
int npair2;                        /* # native contacts                     */
int ndpair;
int spair; 
int ip1[MAXP],ip2[MAXP];           /* list of contacts                      */
int ip3[MAXP],ip4[MAXP];           /* list of contacts                      */
int id1[MAXP],id2[MAXP];           /* list of contacts                      */
int mc1[MAXP],mc2[MAXP];           /* common contacts                       */
int nni1[MAXP],nnj1[MAXP];         /* nearest neighbors                     */
int nni2[MAXP],nnj2[MAXP];
int nni3[MAXP],nnj3[MAXP];
int nni4[MAXP],nnj4[MAXP];
int dual1[MAXP],dual2[MAXP];       /* shared contacts                       */
double kcon_nat[MAXP];             /* contact strengths                     */
double kcon_nat2[MAXP];            /* contact strengths                     */
int dis[N],dis2[N];                /* 1, disordered region, 0 otherwise     */
int qres[N];                       /* charges: D,E -1; K,R +1; others 0     */
/******************/
double distp[MAXP];                /* distances                             */
double distp2[MAXP];               /* distances                             */
double distp3[MAXP];               /* distances                             */
double distp4[MAXP];               /* distances                             */
double distd1[MAXP];               /* distances                             */
double distd2[MAXP];               /* distances                             */
double distg1[MAXP];               /* distances                             */
double distg2[MAXP];               /* distances                             */
double distg3[MAXP];               /* distances                             */
double distg4[MAXP];               /* distances                             */
short cc[N][N];                    /* 1 native contact, 0 otherwise         */
/************* sequence effects *********************************************/
int seq[N];                        /* amino acid sequence                   */
/************* auxiliary ****************************************************/
double cthlim;
/****************************************************************************/
double cellcalc(int inc,int nl,short list[]);
void log_sum_merge(double e1, double e2, double *e,
		   double f1, double f2, double *f,
		   double bet);
void crowd_ecalc(double *e, double *f,double r2,double sigma,double rho);
void bond_ecalc(double *e,double *f,double db);
void bend_ecalc_dis(double *e,double *f,double dtha,double dthb);
void bend_ecalc(double *e,double *f,double dth);
void tors_ecalc_dis(double *e,double *f,double phval);
void tors_ecalc(double *e,double *f,double dph);
double csalt_fac(double csalt,double qi, double qj);
void cont_rep_ecalc(double *e,double *f,double r2,double sig2,double kcont);
void cont_att_ecalc(double kstr,double r,
		    double *eg,double *eg1,double *eg2,
		    double *fg,double *fg1,double *fg2,
		    int im1,int jm1,int im2,int jm2, 
		    double dg,double dg1,double dg2,
		    double *rg1x,double *rg1y,double *rg1z,
		    double *rg2x,double *rg2y,double *rg2z);
void cont_gcalc(double *e, double *f,double r, double r0,double kcont,double ksi);
void cont_ecalc(double *e,double *f,double sig2,double kcont,double r2);
double cont(int iflag);
double cont2(int iflag);
double cont_corr(int iflag);
void add_f(int i,int j,double fr,double rx,double ry,double rz);
/****************************************************************************/
/***** ENERGIES & FORCES ****************************************************/
/****************************************************************************/
void add_f(int i,int j,double fr,double rx,double ry,double rz) {
  fx[i] -= fr * rx;
  fy[i] -= fr * ry;
  fz[i] -= fr * rz;      
  fx[j] += fr * rx;
  fy[j] += fr * ry;
  fz[j] += fr * rz;
}
/****************************************************************************/
void log_sum_merge(double e1, double e2, double *e,
		   double f1, double f2, double *f,
		   double bet) {
  double fmax,wsum,del;

  (*e) = (e1 < e2 ? e1 : e2);
  (*f) = (e1 < e2 ? f1 : f2);

  if ( (del = bet * fabs(e1 - e2)) > 100)
    return ;

  fmax = (e1 < e2 ? f2 : f1);  
  wsum = 1 + exp(-del);
  
  (*e) = (*e) - log(wsum) / bet;
  (*f) = ( (*f) + sgn(fmax) * exp(-del + log(fabs(fmax)) ) ) / wsum;
  
  return;
}
/****************************************************************************/
void crowd_ecalc(double *e, double *f,double r2,double sigma,double rho) {
  double r = sqrt(r2),rcalc,rcalc12;

  rcalc = sigma / (r - rho  + sigma);
  rcalc12 = rcalc * rcalc * rcalc; //rcalc^3
  rcalc12 = rcalc12 * rcalc12;     //rcalc^6
  rcalc12 = rcalc12 * rcalc12;     //rcalc^12

  (*e) = eps * rcalc12;
  (*f) = -12 * (*e) * rcalc / sigma / r;
}
/****************************************************************************/
double crowd_crowd(int iflag){
  int i, j;
  double d, dx, dy, dz, r2;                     /* To clculate distance between crowders */                      
  double e = 0,fcc = 0,ecc = 0;                 /* To calculate force and energy between crowders */
  static double sigma,rho,rcut,rcut2,Dlim,Dlim2;/* To apply limitations and ecalsh penalty */
  FILE *fp;
    
  if (NCR == 0)
    return 0;

  if (iflag < 0) {
    rho = rcrowd + rcrowd; // range crowder-crowder interaction
    sigma = srefcr + srefcr; // softness crowder-crowder interaction 
    rcut = rho + sigma;
    rcut2 = rcut * rcut;
    Dlim = rho - sigma * (1 - pow(eclash / eps,-1./12));
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
    Dlim = rho - sigma * (1 - pow(eclash / eps,-1./12));
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
void bond_ecalc(double *e,double *f,double db) {

  (*e) = kbon * db * db;
  (*f) = - kbon * 2 * db;
}
/****************************************************************************/
double bond(int iflag) {
  int i,j,k;
  double fb,fb1=0,fb2=0;
  double e=0,e1=0,e2=0,et,bb,db,db2,r;
  double dx,dy,dz;
  static double bet=10.0;
  FILE *fp;

  if (FF_BOND == 0)
    return 0;

  if (iflag < 0) {
    if (FF_BOND > 1) printf("<bond> bet %lf\n",bet);
    printf("<bond> kbon %f \n",kbon);    
    return 0;
  }

  if (iflag == 0) {

    for (k = 0; k < NCH; ++k) {
      for (i = iBeg[k]; i < iEnd[k]; ++i) {
	j = i + 1;

	switch (dis[i]) {
	case 1  : {bond_ecalc(&et,&fb,b[i]-bn_dis); break;}
	default : {bond_ecalc(&et,&fb,b[i]-bn[i]); break;}
	}
	
	if (FF_BOND == 2) {

	  switch (dis2[i]) {
	  case 1  : {bond_ecalc(&e2,&fb2,b[i]-bn_dis); break;}
	  default : {bond_ecalc(&e2,&fb2,b[i]-bn2[i]); break;}
	  }

	  log_sum_merge(e1=et,e2,&et,fb1=fb,fb2,&fb,bet);
	}

	e += et;

	add_f(i,j,fb,bx[i],by[i],bz[i]);
      }
    }

    if (FF_DISULF) {
    
      for (k = 0; k < ndpair; ++k) {
	i = id1[k]; j = id2[k];
	
	db = (r = sqrt(vec2(i,j,&dx,&dy,&dz))) - distd1[k];
	bond_ecalc(&et,&fb,db);
	
	if (FF_DISULF == 2) {
	  if ( db * (db2 = r - distd2[k]) < 0 ) continue;
	  if ( db * (db2 - db) < 0 ) bond_ecalc(&et,&fb,db2);
	}
	
	e += et;
	
	add_f(i,j,fb,dx,dy,dz);
      }
    }
    
    return e;
  }

  if (iflag > 0) {
    fp = fopen("results/param_check/bond_energy.plot","w");

    for (k = 0; k < NCH; ++k) {
      for (i = iBeg[k]; i < iEnd[k]; ++i) {
    
	for (bb=3.0; bb<5.0; bb+=0.01) {

	  switch (dis[i]) {
	  case 1  : {bond_ecalc(&e1,&fb1,bb-bn_dis); break;}
	  default : {bond_ecalc(&e1,&fb1,bb-bn[i]); break;}
	  }	  

	  if (FF_BOND == 2) {

	    switch (dis2[i]) {
	    case 1  : {bond_ecalc(&e2,&fb2,bb-bn_dis); break;}
	    default : {bond_ecalc(&e2,&fb2,bb-bn2[i]); break;}
	    }

	    log_sum_merge(e1,e2,&et,fb1,fb2,&fb,bet);

	    fprintf(fp,"%i %lf  %lf %lf %lf  %lf %lf %lf\n",
		    i,bb,e1,e2,et,fb1,fb2,fb);
	    continue;
	  }

	  fprintf(fp,"%i %lf  %lf %lf\n",i,bb,e1,fb1);
	}
      }
    }

    fclose(fp);
  }

  return 0;
}
/****************************************************************************/
void bend_ecalc_dis(double *e,double *f,double dtha,double dthb) {
  /* Potential V(x) = -ln[ exp(-(x-ta)^2/(2ka^2) + exp(-(x-tb)^2/(2kb^2) ]  
     where ta,tb are reference values and ka,kb are parameters. V(x) is
     shifted by additive constant eth0. Returns (*e) = V(x) and
     (*f) = - V'(x) */

  double ga,gb;
  
  dtha /= ksi_disa;
  dthb /= ksi_disb;
  
  ga = exp( - dtha * dtha / 2 );
  gb = exp( - dthb * dthb / 2 );
    
  (*e) = - eps * log(ga + gb) + eth0; 
  (*f) = - eps * (dtha / ksi_disa * ga + dthb / ksi_disb * gb) / (ga + gb);  
}
/****************************************************************************/
void bend_ecalc(double *e,double *f,double dth) {

  (*e) = kth * dth * dth;
  (*f) = - kth * 2 * dth;
}
/****************************************************************************/
double bend(int iflag) {
  int i,j,k,l;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z,b2;
  double dix,diy,diz;
  double dkx,dky,dkz;
  double cth,sth;
  double e = 0,e1,e2,et,fben,fben1,fben2,d;
  static double bet = 5.0;
  FILE *fp;

  if (FF_BEND == 0)
    return 0;

  if (iflag < 0) {
    if (FF_BEND > 1) printf("<bend> bet %lf\n",bet);
    printf("<bend> kth %lf \n",kth);
    return 0;
  }
  
  if (iflag == 0) {

    for (l = 0; l < NCH; ++l) {
      for (i = iBeg[l]; i < iEnd[l] - 1; ++i) {
	j = i + 1;
	k = i + 2;
	
	switch (dis[j]) {
	case 1  : {bend_ecalc_dis(&et,&fben,th[j]-thn_disa,th[j]-thn_disb); break;}
	default : {bend_ecalc(&et,&fben,th[j]-thn[j]); break;}
	}
	
	if (FF_BEND == 2) {

	  switch (dis2[j]) {
	  case 1  : {bend_ecalc_dis(&e2,&fben2,th[j]-thn_disa,th[j]-thn_disb); break;}
	  default : {bend_ecalc(&e2,&fben2,th[j]-thn2[j]); break;}
	  }
	  
	  log_sum_merge(e1=et,e2,&et,fben1=fben,fben2,&fben,bet);
	}
	
	e += et;
	
	cth = cos(th[j]);
	sth = sin(th[j]);
	
	b1x = bx[i];
	b1y = by[i];
	b1z = bz[i];
	b1 = b[i];

	b2x = bx[j];
	b2y = by[j];
	b2z = bz[j];
	b2 = b[j];
	
	dix = -(b2x+cth*b1x)/sth/b1;
	diy = -(b2y+cth*b1y)/sth/b1;
	diz = -(b2z+cth*b1z)/sth/b1;
	dkx =  (b1x+cth*b2x)/sth/b2;
	dky =  (b1y+cth*b2y)/sth/b2;
	dkz =  (b1z+cth*b2z)/sth/b2;
	
	fx[i] += fben*dix;
	fy[i] += fben*diy;
	fz[i] += fben*diz;
	fx[j] += fben*(-dix-dkx);
	fy[j] += fben*(-diy-dky);
	fz[j] += fben*(-diz-dkz);
	fx[k] += fben*dkx;
	fy[k] += fben*dky;
	fz[k] += fben*dkz;
      }
    }

    return e;
  }

  if (iflag > 0)  {
    fp = fopen("results/param_check/bend_energy.plot","w");
    
    for (l = 0; l < NCH; ++l) {
      for (i = iBeg[l]; i < iEnd[l]-1; ++i) {
	j = i + 1;

	for (d = 0; d < pi; d += pi/180) {

	  switch (dis[j]) {
	  case 1  : {bend_ecalc_dis(&e1,&fben1,d-thn_disa,d-thn_disb); break;}
	  default : {bend_ecalc(&e1,&fben1,d-thn[j]); break;}
	  }
	  
	  if (FF_BEND == 2) {

	    switch (dis2[j]) {
	    case 1  : {bend_ecalc_dis(&e2,&fben2,d-thn_disa,d-thn_disb); break;}
	    default : {bend_ecalc(&e2,&fben2,d-thn2[j]); break;}
	    }
	    
	    log_sum_merge(e1,e2,&et,fben1,fben2,&fben,bet);

	    fprintf(fp,"%i %lf  %lf %lf %lf  %lf %lf %lf\n",j,d*rad2deg,e1,e2,et,
		    fben1,fben2,fben);

	    continue;
	  }

	  fprintf(fp,"%i %lf  %lf %lf\n",j,d*rad2deg,e1,fben1);
	}
      }
    }
  }
  
  return 0;
}
/****************************************************************************/
void tors_ecalc_dis(double *e,double *f,double phval) {
  /* V(x) = k1(1-cos(dph)) + k2(1-cos(2dph)) + k3(1-cos(3dph)) + eth0, where 
     eth0 is a shift set so that the minimum of V is zero. Returns
     (*e) = V(x) and (*f) = -V'(x). */

  double dph1 = phval - phn_dis1;
  double dph2 = 2 * (phval - phn_dis2);
  double dph3 = 3 * (phval - phn_dis3);
  
  (*e) = ( kph_dis1 * (1 - cos(dph1)) +
	   kph_dis2 * (1 - cos(dph2)) +
	   kph_dis3 * (1 - cos(dph3)) + eph0 );
  (*f) = - ( kph_dis1 * sin(dph1) +
	     kph_dis2 * 2 * sin(dph2) +
	     kph_dis3 * 3 * sin(dph3) ); 
}
/****************************************************************************/
void tors_ecalc(double *e,double *f,double dph) {
  /* Torsion potential V(x) = kph1*(1-cos(x)) + kph3*(1-cos(3x)).
     Returns (*e) = V(x) and (*f) = -V'(x). */

  (*e) = kph1 * (1 - cos(dph)) + kph3 * (1 - cos(3 * dph));
  (*f) = - kph1 * sin(dph) - 3 * kph3 * sin(3 * dph); 
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
  double cth1,cth2,sth1,sth2;
  static double bet=5.0;
  
  if (FF_TORS == 0)
    return 0;

  if (iflag < 0) {
    if (FF_TORS == 2) printf("<torsion> bet %lf\n",bet);
    printf("<torsion> kph1 %lf kph3 %lf\n",kph1,kph3);
    return 0;
  }
  
  if (iflag == 0) {

    for (m = 0; m < NCH; ++m) {
      for (i = iBeg[m]; i < iEnd[m] - 2; ++i) {
	j = i + 1;
	k = i + 2;
	l = i + 3;

	switch (dis[j]) {
	case 1  : {tors_ecalc_dis(&et,&fph,ph[j]); break;}
	default : {tors_ecalc(&et,&fph,ph[j]-phn[j]); break;}
	}
	
	if (FF_TORS == 2) {
	  
	  switch (dis2[j]) {
	  case 1  : {tors_ecalc_dis(&e2,&fph2,ph[j]); break;}
	  default : {tors_ecalc(&e2,&fph2,ph[j]-phn2[j]); break;}
	  }
	  
	  log_sum_merge(e1=et,e2,&et,fph1=fph,fph2,&fph,bet);
	}
	
	e += et;
	
	cth1=cos(th[j]);
	sth1=sin(th[j]);
	cth2=cos(th[k]);
	sth2=sin(th[k]);
	
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
    }
    
    return e;
  }

  if (iflag > 0) { /* print potential function */
    FILE *fp;
    double d;
    
    fp = fopen("results/param_check/tors_energy.plot","w");

    for (m = 0; m < NCH; ++m) {
      for (i = iBeg[m]; i < iEnd[m] - 2; ++i) {
	j = i + 1;
	k = i + 2;
	l = i + 3;

	for (d = -pi; d < pi; d += pi/360) {

	  switch (dis[j]) {
	  case 1  : {tors_ecalc_dis(&e1,&fph1,d); break;}
	  default : {tors_ecalc(&e1,&fph1,d-phn[j]); break;}
	  }
	  
	  if (FF_TORS == 2) {

	    switch (dis2[j]) {
	    case 1  : {tors_ecalc_dis(&e2,&fph2,d); break;}
	    default : {tors_ecalc(&e2,&fph2,d-phn2[j]); break;}
	    }
	    
	    log_sum_merge(e1,e2,&et,fph1,fph2,&fph,bet);

	    fprintf(fp,"%i %lf %lf %lf %lf %lf %lf %lf \n",j,d*rad2deg,e1,e2,et,
		    fph1,fph2,fph);

	    continue;
	  }
	  
	  fprintf(fp,"%i %lf %lf %lf \n",j,d*rad2deg,e1,fph1);
	}
      }
    }
    
    fclose(fp);
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

  if (FF_EXVOL == 0 || N == 0)
    return 0;

  if (iflag < 0) {
    for (i = 0; i < MAXCELL; ++i) cell[i] = -1;

    ns = BOX / cut;
    h2 = ns*ns; 
    h3 = h2*ns;
    cutg = BOX/(double)ns;
    cellcalc(-1,0,list);

    if (ns * ns * ns > MAXCELL) {printf("# cells %i > MAXCELL %i\n",ns*ns*ns,MAXCELL); exit(-1);}
    printf("<exvol>  #cells %i ns %i cutg %f\n",ns*ns*ns,ns,cutg);
    printf("<exvol>  eps %f sigsa %f krep %f\n",eps,sigsa,krep);

    return 0;
  }

  in2box();
  for (i = N-1; i >= 0; --i) {
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
  static double krep12;
  
  if (inc<0) {
    cut2=cut*cut;
    sigsa2=sigsa*sigsa;
    asa=-7*pow(sigsa2/cut2,6.0);
    bsa=6*pow(sigsa2/cut2,6.0)/cut2;
    bsa2=2*bsa;
    krep12=12*krep;
    printf("<cellcalc> asa %e bsa %e\n",asa,bsa);
    return 0.0;
  }

  for (n=0;n<inc;n++) {
    i=list[n]; ic=a2c[i];
    for (m=n+1;m<nl;m++) {
      j=list[m];
      if ( (abs(i-j) < 4 && ic == a2c[j]) || cc[i][j] == 1 || cc[j][i] == 1) continue;
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
void cont_rep_ecalc(double *e,double *f,double r2,double sig2,double kcont) {
  /* Weeks-Chandler-Andersen-type repulsion V(r) */
  /* Returns (*e) = V(r) and (*f) = - V'(r) / r  */ 
  double r4,r6;

  if (r2 >= sig2) {
    (*e) = (*f) = 0;
    return ;
  }
  
  r6 = sig2 / r2; 
  r4 = r6 * r6; 
  r6 = r4 * r6;
  
  (*e) = kcont * ( r6 * (r6 - 2) + 1 );
  (*f) = kcont * 12 * r6 * (r6 - 1) / r2 ;
}
/****************************************************************************/
double csalt_fac(double csalt,double qi, double qj) {
  /* salt-dependent scaling factor applied to the attractive part of native
     contact interactions */
  return (csalt - 1) * qi * qj + 1;
}
/****************************************************************************/
void cont_gcalc(double *e, double *f, double r, double r0,
		double kcont,double ksi) {
  /* Gaussian function V(r) = kcont * exp(-(r-r0)^2 / 2ksi) */
  /* Returns (*e) = V(r) and (*f) = - V'(r) / r */

  double dr = r - r0;
  double g = exp(- dr * dr / 2 / ksi);

  (*e) = kcont * g;
  (*f) = kcont * dr * g / ksi / r;
}
/****************************************************************************/
void cont_att_ecalc(double kstr,double r,
		    double *eg,double *eg1,double *eg2,
		    double *fg,double *fg1,double *fg2,
		    int im1,int jm1,int im2,int jm2, 
		    double dg,double dg1,double dg2,
		    double *rg1x,double *rg1y,double *rg1z,
		    double *rg2x,double *rg2y,double *rg2z) {
  double rg1,rg2;
  
  cont_gcalc(eg,fg,r,dg,kstr,ksi1);
    
  if (im1 >= 0 && jm1 >= 0) {
    rg1 = vec2(im1,jm1,rg1x,rg1y,rg1z);
    cont_gcalc(eg1,fg1,sqrt(rg1),dg1,1.0,ksi2);
  } else *eg1 = *fg1 = 1.0;

  if (im2 >= 0 && jm2 >= 0) {
    rg2 = vec2(im2,jm2,rg2x,rg2y,rg2z);
    cont_gcalc(eg2,fg2,sqrt(rg2),dg2,1.0,ksi2);
  } else *eg2 = *fg2 = 1.0;
  
  (*fg1) = - (*eg) * (*fg1) * (*eg2);
  (*fg2) = - (*eg) * (*eg1) * (*fg2);
  (*eg) = - (*eg) * (*eg1) * (*eg2);
  (*fg) = - (*fg) * (*eg1) * (*eg2);
}
/****************************************************************************/
void cont_ecalc(double *e,double *f,double sig2,double kcont,double r2) {
  /* Lennard-Jones-type potential V(r) */
  /* Returns (*e) = V(r) and (*f) = -V'(r) / r */

  double r4,r6;

  r6 = sig2 / r2; 
  r4 = r6 * r6; 
  r6 = r4 * r6;
  
  (*e) = kcont * ( r6 * (5 * r6 - 6 * r4)  );
  (*f) = kcont * 60 * ( r6 * (r6 - r4) / r2  );
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

  if (FF_CONT == 0 || N == 0)
    return 0;
  
  if (iflag < 0) {
    rcut_mb = 3 * ksi1;
    
    if (FF_MULTIBODY)  {
      printf("<cont> FF_MULTIBODY %i\n",FF_MULTIBODY);
      printf("<cont> ksi1 %lf ksi2 %lf\n",ksi1,ksi2);
      printf("<cont> cutoff %lf \n",rcut_mb);
    }
    
    if (FF_CONT == 2) {
      cont2(-1);
      cont_corr(-1);
    }

    fp1 = fopen("results/param_check/cont_param.out","w");
    for (m = 0; m < npair; m++) {
      if (dual1[m] > 0) continue;
      i = ip1[m]; j = ip2[m];
      fprintf(fp1,"%i %i %lf %lf %i %i %i %i %lf %lf \n",i,j,
	      kcon_nat[m],distp[m],
	      nni1[m],nnj1[m],nni2[m],nnj2[m],
	      distg1[m],distg2[m]);
    }
    fclose(fp1);
    
    return 0;
  }

  if (iflag == 0) {

    for (m = 0; m < npair; m++) {
      if (dual1[m] > 0) continue;
      i = ip1[m]; j = ip2[m];

      r2 = vec2(i,j,&rx,&ry,&rz);

      if (!FF_MULTIBODY) {
	if ( r2 > 4 * (sig2 = distp2[m]) ) continue;
	cont_ecalc(&er,&fr,sig2,kcon_nat[m],r2);
	e += er;
	add_f(i,j,fr,rx,ry,rz);
      }
      
      if (FF_MULTIBODY) {
	if ( (r = sqrt(r2)) > distp[m] + rcut_mb ) continue;

	im1 = nni1[m]; jm1 = nnj1[m];
	im2 = nni2[m]; jm2 = nnj2[m];

	cont_rep_ecalc(&er,&fr,r2,distp2[m],krep);
	cont_att_ecalc(kcon_nat[m],r,
		       &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
		       distp[m],distg1[m],distg2[m],
		       &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);

	e += er + eg;
	fr = fr + fg;

	add_f(i,j,fr,rx,ry,rz);
	add_f(im1,jm1,fg1,rg1x,rg1y,rg1z);
	add_f(im2,jm2,fg2,rg2x,rg2y,rg2z);
      }

    }

    Econ1 = e;
    
    if (FF_CONT == 2) {
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
	  cont_ecalc(&er,&fr,sig2,kcon_nat[m],r2);
	  fprintf(fp1,"%i %i %i %lf  %lf\n",m,i,j,d,er);
	}

	if (FF_MULTIBODY) {
	  if ( (r = sqrt(r2)) > distp[m] + rcut_mb ) continue;
	  im1 = jm1 = im2 = jm2 = -1;
	  
	  cont_rep_ecalc(&er,&fr,r2,distp2[m],krep);
	  cont_att_ecalc(kcon_nat[m],r,
		       &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
		       distp[m],distg1[m],distg2[m],
		       &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);

	  fprintf(fp1,"%i %i %i %lf  %lf %lf %lf \n",m,i,j,d,e+eg,er,eg);
	} 

      }
    }

    fclose(fp1);
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

  if (FF_CONT < 2)
    return 0;
  
  if (iflag < 0) {
    rcut_mb = 3 * ksi1;

    if (FF_MULTIBODY == 1) 
      printf("<cont2> cutoff %lf \n",rcut_mb);

    fp1 = fopen("results/param_check/cont_param2.out","w");
    for (m = 0; m < npair2; m++) {
      if (dual2[m] > 0) continue;
      i = ip3[m]; j = ip4[m];
      fprintf(fp1,"%i %i %lf %lf %i %i %i %i %lf %lf\n",i,j,
	      kcon_nat2[m],distp3[m],
	      nni3[m],nnj3[m],nni4[m],nnj4[m],
	      distg3[m],distg4[m]);
    }
    fclose(fp1);

    return 0;
  }

  if (iflag == 0) {

    for (m = 0; m < npair2; m++) {
      if (dual2[m] > 0) continue;
      i = ip3[m]; j = ip4[m];

      r2 = vec2(i,j,&rx,&ry,&rz);

      if (!FF_MULTIBODY) {
	if ( r2 > 4 * (sig2 = distp4[m]) ) continue;
	cont_ecalc(&er,&fr,sig2,kcon_nat2[m],r2);
	e += er;
	add_f(i,j,fr,rx,ry,rz);
      }

      if (FF_MULTIBODY) {
	if ( (r = sqrt(r2)) > distp3[m] + rcut_mb ) continue;
	im1 = nni3[m]; jm1 = nnj3[m];
	im2 = nni4[m]; jm2 = nnj4[m];

	cont_rep_ecalc(&er,&fr,r2,distp4[m],krep);
	cont_att_ecalc(kcon_nat2[m],r,
		       &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
		       distp3[m],distg3[m],distg4[m],
		       &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);

	e += er + eg;
	fr = fr + fg;
	add_f(i,j,fr,rx,ry,rz);
	add_f(im1,jm1,fg1,rg1x,rg1y,rg1z);
	add_f(im2,jm2,fg2,rg2x,rg2y,rg2z);
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
	  cont_ecalc(&er,&fr,sig2,kcon_nat2[m],r2);
	  fprintf(fp1,"%i %i %i %lf  %lf\n",m,i,j,d,er);
	}

	if (FF_MULTIBODY) {
	  if ( (r = sqrt(r2)) > distp3[m] + rcut_mb ) continue;
	  im1 = jm1 = im2 = jm2 = -1;

	  cont_rep_ecalc(&er,&fr,r2,distp4[m],krep);
	  cont_att_ecalc(kcon_nat2[m],r,
			 &eg,&eg1,&eg2,&fg,&fg1,&fg2,im1,jm1,im2,jm2,
			 distp3[m],distg3[m],distg4[m],
			 &rg1x,&rg1y,&rg1z,&rg2x,&rg2y,&rg2z);

	  fprintf(fp1,"%i %i %i %lf  %lf %lf %lf \n",m,i,j,d,er+eg,er,eg);
	} 
      }

    }
    
    fclose(fp1);
  }
  
  return 0;
}
/****************************************************************************/
/*double cont_corr_old(int iflag) {
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

  if (FF_CONT < 2)
    return 0;
  
  if (!FF_MULTIBODY) {
    printf("    cont_corr() only implemented for the case FF_MULTIBODY %i\nExiting...\n",FF_MULTIBODY);
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
		kcon_nat[m],kcon_nat2[n],distp[m],distp3[n],
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

	cont_rep_ecalc(&er,&fr,r2,dist_rep1[m],krep);

	cont_att_ecalc(kcon_nat[m],r,
		       &egA,&eg1A,&eg2A,&fgA,&fg1A,&fg2A,im1A,jm1A,im2A,jm2A,
		       distp[m],distg1[m],distg2[m],
		       &rg1xA,&rg1yA,&rg1zA,&rg2xA,&rg2yA,&rg2zA);
	
	cont_att_ecalc(kcon_nat2[n],r,
		       &egB,&eg1B,&eg2B,&fgB,&fg1B,&fg2B,im1B,jm1B,im2B,jm2B,
		       distp3[n],distg3[n],distg4[n],
		       &rg1xB,&rg1yB,&rg1zB,&rg2xB,&rg2yB,&rg2zB);

	if (egA < egB) {
	  e += er + egA;
	  fr = fr + fgA;
	  add_f(i,j,fr,rx,ry,rz);
	  add_f(im1A,jm1A,fg1A,rg1xA,rg1yA,rg1zA);
	  add_f(im2A,jm2A,fg2A,rg2xA,rg2yA,rg2zA);
	} else {
	  e += er + egB;
	  fr = fr + fgB;
	  add_f(i,j,fr,rx,ry,rz);
	  add_f(im1B,jm1B,fg1B,rg1xB,rg1yB,rg1zB);
	  add_f(im2B,jm2B,fg2B,rg2xB,rg2yB,rg2zB);
	} 
      }
    }
    
    return e;
  }
  
  if (iflag > 0) {    
    return 0;
  }
  
  return 0;
  } */
/****************************************************************************/
double cont_corr(int iflag) {
  int s,m,n,i,j;
  double r,r2,rx,ry,rz;
  double e = 0,fr,erA,frA,erB,frB;
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

  if (FF_CONT < 2)
    return 0;
  
  if (!FF_MULTIBODY) {
    printf("    cont_corr() only implemented for the case FF_MULTIBODY %i\nExiting...\n",FF_MULTIBODY);
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

	fprintf(fp1,"%i %i kcon %lf %lf dist %lf %lf %i %i %i %i %i %i %i %i distg %lf %lf %lf %lf\n",i,j,
		kcon_nat[m],kcon_nat2[n],distp[m],distp3[n],
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
      r = sqrt(r2);

      if (FF_MULTIBODY) {
	if ( (r > distp[m] + rcut_mb)  && (r > distp3[n] + rcut_mb) ) continue;

	im1A = nni1[m]; jm1A = nnj1[m];
	im2A = nni2[m]; jm2A = nnj2[m];

	cont_rep_ecalc(&erA,&frA,r2,distp2[m],krep);
	cont_att_ecalc(kcon_nat[m],r,
		       &egA,&eg1A,&eg2A,&fgA,&fg1A,&fg2A,im1A,jm1A,im2A,jm2A,
		       distp[m],distg1[m],distg2[m],
		       &rg1xA,&rg1yA,&rg1zA,&rg2xA,&rg2yA,&rg2zA);
	
	im1B = nni3[n]; jm1B = nnj3[n];
	im2B = nni4[n]; jm2B = nnj4[n];

	cont_rep_ecalc(&erB,&frB,r2,distp4[n],krep);
	cont_att_ecalc(kcon_nat2[n],r,
		       &egB,&eg1B,&eg2B,&fgB,&fg1B,&fg2B,im1B,jm1B,im2B,jm2B,
		       distp3[n],distg3[n],distg4[n],
		       &rg1xB,&rg1yB,&rg1zB,&rg2xB,&rg2yB,&rg2zB);

	if (erA + egA < erB + egB) {
	  e += erA + egA;
	  frA = frA + fgA;
	  add_f(i,j,frA,rx,ry,rz);
	  add_f(im1A,jm1A,fg1A,rg1xA,rg1yA,rg1zA);
	  add_f(im2A,jm2A,fg2A,rg2xA,rg2yA,rg2zA);
	} else {
	  e += erB + egB;
	  frB = frB + fgB;
	  add_f(i,j,frB,rx,ry,rz);
	  add_f(im1B,jm1B,fg1B,rg1xB,rg1yB,rg1zB);
	  add_f(im2B,jm2B,fg2B,rg2xB,rg2yB,rg2zB);
	} 
      }
    }
    
    return e;
  }
  
  if (iflag > 0) {    
    fp1 = fopen("results/param_check/cont_energy_corr.plot","w");
    
    for (s = 0; s < spair; ++s) {
      m = mc1[s];
      n = mc2[s];

      i = ip1[m];
      j = ip2[m];

      for (r = 3.0; r < 20; r += 0.01) {
	r2 = r * r;
	
	if ( r > distp[m] + rcut_mb && r > distp3[n] + rcut_mb ) continue;
	
	if (FF_MULTIBODY) {
	  im1A = jm1A = im2A = jm2A = -1;
	  im1B = jm1B = im2B = jm2B = -1;

	  cont_rep_ecalc(&erA,&frA,r2,distp2[m],krep);
	  cont_att_ecalc(kcon_nat[m],r,
			 &egA,&eg1A,&eg2A,&fgA,&fg1A,&fg2A,im1A,jm1A,im2A,jm2A,
			 distp[m],distg1[m],distg2[m],
			 &rg1xA,&rg1yA,&rg1zA,&rg2xA,&rg2yA,&rg2zA);
	  
	  cont_rep_ecalc(&erB,&frB,r2,distp4[n],krep);
	  cont_att_ecalc(kcon_nat2[n],r,
			 &egB,&eg1B,&eg2B,&fgB,&fg1B,&fg2B,im1B,jm1B,im2B,jm2B,
			 distp3[n],distg3[n],distg4[n],
			 &rg1xB,&rg1yB,&rg1zB,&rg2xB,&rg2yB,&rg2zB);
	  
	  if (erA + egA < erB + egB) {
	    e = erA + egA;
	    fr = frA + fgA;
	  } else {
	    e = erB + egB;
	    fr = frB + fgB;
	  } 
	
	  fprintf(fp1,"%i %i %i %i %lf  %lf %lf %lf  %lf %lf %lf\n",
		  m,n,i,j,r,erA+egA,erB+egB,e,frA+fgA,frB+fgB,fr);
	}
      }
    }
    
    fclose(fp1);
  }
  
  return 0;
}
/****************************************************************************/
