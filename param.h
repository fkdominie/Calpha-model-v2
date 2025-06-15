/****************************************************************************/
/************* Interaction parameters ***************************************/
/****************************************************************************/
const double eps=1.0;              /* energy unit                           */
/* protein */
double kbon=100.0;                 /* bond strength                         */
double kth=20.0;                   /* angle strength                        */
double kph1=1.0;                   /* torsion strength                      */
double kph3=0.5;                   /* torsion strength                      */
double kcon=1.0;                   /* contact strength                      */
double krep=1.0;                   /* excluded volume strength              */
double sigsa=4.0;                  /* bead diameter                         */
double cut=8.0;                    /* cutoff distance excluded volume       */
double ksi1=1.0;                   /* contact parameter                     */
double ksi2=25.0;                  /* contact parameter                     */
/* disordered regions */
double bn_dis=3.8;                 /* bond length                           */
double thn_disa=93.0;              /* angle (alpha)                         */
double thn_disb=117.0;             /* angle (sheet)                         */ 
double ksi_disa=4.0;               /* angle parameter                       */
double ksi_disb=12.0;              /* angle parameter                       */
double eth0=0.13;                  /* angle energy shift                    */
double kph_dis1=0.5;              /* torsion energy strength               */ 
double kph_dis2=0.8;              /* torsion energy strength               */
double kph_dis3=0.25;             /* torsion energy strength               */
double phn_dis1=180;               /* torsion reference angle               */
double phn_dis2=55;                /* torsion reference angle               */
double phn_dis3=25;                /* torsion reference angle               */
double eph0=-0.44;                 /* torsion shift                         */
double fscal=1.16;                 /* torsion scale factor                  */
/* crowders: */
const double rcrowd = 12.0;        /* Crowder radius                        */
const double srefcr = 3.0;         /* Crowder repulsion softness (sigma)    */
double epsilonrep = 1.0;           /* Crowder repulsion strength            */
double eclash = 1e6;               /* Crowder repulsion ceiling             */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/



