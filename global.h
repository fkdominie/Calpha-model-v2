/************* geometry *****************************************************/
extern double x[],y[],z[];     
extern double xb[],yb[],zb[];
extern double xcr[],ycr[],zcr[];
extern double xcro[],ycro[],zcro[];
extern double b[];             
extern double th[];            
extern double ph[];            
extern double bx[],by[],bz[];  
extern double sx[],sy[],sz[];
extern double boxhf; 
extern int iBeg[],iEnd[];
extern int a2c[];
extern int *a2p;
extern int *c2p;
extern int NPROT;
extern const double vbox;
extern int seq[];                   
extern double kap[];
extern double charge[];
/************* energies and forces ******************************************/
extern double Ekin,Epot,Eben,Ebon,Erep;  
extern double Etor,Econ,Ehp,Eel;
extern double Econ1,Econ2,Ecorr;
extern double fx[],fy[],fz[];         
extern double fxo[],fyo[],fzo[];      
extern double frdx[],frdy[],frdz[];   
extern double frdxo[],frdyo[],frdzo[];
/************* MD parameters ************************************************/
extern long imd,imd0;
extern const double tau;          
extern double vx[],vy[],vz[];         
extern double gam;                    
extern double dt;                     
extern double tconstbd[];
extern double c1,c2,c3;
extern double mbd;
/************* MC parameters ************************************************/
extern double g[];   
extern double beta[];
extern int ind;
extern long int nflp,accflp;
/************* interactions *************************************************/
extern const double eps;          
extern double kbon;             
extern double kth;               
extern double kph1;
extern double kph3;
extern double kcon;
extern double khp;
extern double kel;               
extern double krep;               
extern double sigsa;              
extern double cut;                
extern double sighp;              
extern double cuthp;
extern double lambda_D;
extern double cutel;
extern double ksi1,ksi2;
extern double thn_disa,thn_disb,ksi_disa,ksi_disb;
extern double phn_dis1,phn_dis2,phn_dis3;
extern double kph_dis1,kph_dis2,kph_dis3;
extern double eth0,eph0,fscal;
extern int qres[];
extern double csalt;
/************* native structure *********************************************/
extern double xnat[],ynat[],znat[];
extern double xnat2[],ynat2[],znat2[];
extern double bn[],thn[],phn[];       
extern double bn2[],thn2[],phn2[];
extern int npair,npair2;
extern int ndpair;
extern int spair;
extern int qpair;
extern int zpair;                   
extern int ip1[],ip2[];        
extern int ip3[],ip4[];
extern int iq1[],iq2[];
extern int iz1[],iz2[];
extern int id1[],id2[];
extern int mc1[],mc2[];
extern int nni1[],nnj1[];
extern int nni2[],nnj2[];
extern int nni3[],nnj3[];
extern int nni4[],nnj4[];
extern int dual1[],dual2[];
extern int nat1[],nat2[],nat3[],nat4[],nat5[];
extern int link[];
extern double kcon_nat[];     
extern double kcon_nat2[];       
extern double kbond[],kbond2[];
extern double kbend[],kbend2[];
extern double ktor[][3],ktor2[][3];
extern int dis[],dis2[];
extern double distp[];           
extern double distp2[];          
extern double distp3[];           
extern double distp4[];
extern double distd1[]; 
extern double distd2[]; 
extern double distg1[],distg2[];
extern double distg3[],distg4[];
extern double iqkap[];
extern double izcharge[];
extern short cc[][N];
/************* Crowders *****************************************************/
extern double Ecc, Ecb;                        /* Energy terms related to the crowders */
extern const double rcrowd;                   /* radius of crowders  */
extern const double Phicr;                    /*  Volume fraction of crowders */
extern const double srefcr;                  /*  Coefficient of softness  */
extern double eclash;    
extern double fxc[],fyc[],fzc[];          /* conformational force crowders                 */
extern double fxco[],fyco[],fzco[];       /* conformational force old crowders             */
extern double frcdx[],frcdy[],frcdz[];    /* random force crowders                          */
extern double frcdxo[],frcdyo[],frcdzo[]; /* random force old crowders                     */
extern double vxc[],vyc[],vzc[];         /* velocity components of the crowders.  */
extern double mcr;                       /* Mass of crowders                     */
extern double gamcr;
extern double c1cr,c2cr,c3cr;
extern double tconstcr[];
/************* miscellaneous ************************************************/
extern double pi,pi2,pid2;
extern double deg2rad,rad2deg;
extern double cthlim;
extern int carterr,therr;
extern long seed,orig_seed;
extern char OUTDIR[];
extern int FANALYS;
extern FILE *fp_log;
/****************************************************************************/
/* misc.c */
void printinfo(void);
void ramachan(char *fn,double b,double th,double ph);
void runtime(long it,double o[]);
void averages(double so[][NOBS]);
void write_conf(char *fn,char *fdir,char *fmode);
void write_momenta(char *fn,char *fdir,char *fmode);
void write_forces(char *fname,char *fdir,char *fmode);
void write_data(char *fname,char *fdir,long int imd,long seed,long orig);
void write_checkpnt(void);
int read_conf(int iflag,char *fn,char *fdir);
int read_momenta(char *fn,char *fdir);
int read_forces(char *fname,char *fdir);
int read_data(char *fname,char *fdir,long *imd,long *seed,long *orig);
int read_checkpnt(void);
void check_distances(int *ip1,int *ip2,int n,double *dist2,int *dual,
		     double *xn,double *yn,double *zn);
void write_shared_dist(char *fn,int *mc1,int *mc2,int spair);
int get_shared_contacts(int *mc1,int *mc2,int *dual1,int *dual2);
void get_natdist(double *dist,double *dist2,int *ip1,int *ip2,int n,
		 double *xn,double *yn,double *zn);
void get_nndist(double *distg1,double *distg2,
		int *ip1,int *ip2,int n,
		int *nni1,int *nnj1,int *nni2,int *nnj2,
		double *xn,double *yn,double *zn);
void write_natdist(char *fn,double *dist,int n,int *ip1,int *ip2);
void write_bonded_param(double *bn,double *thn,double *phn,char *fn);
int relax_chains(int ich);
int relax_crowders(void);
void init(int iflag);
void update_g(double so[][NOBS],double g[]);
/* energy.c */
double crowd_crowd(int iflag);
double crowd_bead(int iflag);
double bond(int iflag);
double bend(int iflag);
double torsion(int iflag);
double csalt_fac(double csalt,double qi, double qj);
double cont(int iflag);
double cont2(int iflag);
double cont_corr(int iflag);
double hp(int iflag);
double el(int iflag);
double exvol(int iflag);
/* obs.c */
void center_of_mass2(int ich,double *xcm,double *ycm, double *zcm);
void msd(int iflag);
void hist_sq(int iflag,int ind);
void histo_bond(int iflag);
void histo_bend(int iflag);
void histo_tors(int iflag,int ia);
void histoe(int iflag,double x);
void histoq1(int iflag, int x);
void histormsd(int iflag,double x);
void historg(int iflag,double x);
void histo_cont1(int iflag, int ind, int n);
void histo_cont2(int iflag, int ind, int n);
void histod1(int iflag,int itmp,double x);
void histod2(int iflag,int itmp,double x);
void histod3(int iflag,int itmp,double x);
void rate(int iflag,int nc);
double gyr2(int a1,int a2);
double dcm(void);
void histo_contmap(int iflag, int ind);
void histo_contmap2(int iflag, int ind);
int contacts(int ilst[],int jlst[],double d2[],int nlst,int ch1,int ch2);
int ncont_map1(int ch);
int ncont_map2(int ch);
int ncont_map1_inter(int ch1,int ch2);
int ncont_map2_inter(int ch1,int ch2);
/* sampling */
void tflip(double e);
void mdstep(void);
/* geometry.c */
double vec2(int i,int j,double *rx,double *ry,double *rz);
double dist2(int i,int j);
void trans(int ic,double dx,double dy,double dz);
void trans_cr(int ic,double dx,double dy,double dz);
void bc(double *x);
void in2box(void);
void ch2box(int ich);
void cr2box(int icr);
void build_protein_maps(const char *spec);
//void ch2boxcr(int icr);
void dof2cart(int iflag);
void cart2dof(void);
/* utils.c */
int read_native(char *fn,double *xr,double *yr,double *zr);
int read_contacts(char *fn,int *ip1,int *ip2);
void dumppdb(char *fn,double *o,int nobs);
double rmsd_calc(double *x1,double *y1,double *z1,
		 double *x2,double *y2,double *z2,
		 int i1,int i2);
double gasdev2(void);
double ran3n(long *idum);
/****************************************************************************/
