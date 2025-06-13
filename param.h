/****************************************************************************/
/************* Interaction parameters ***************************************/
/****************************************************************************/
const double eps=1.0;              /* energy unit                           */
/* protein */
double kbon=100.0;                 /* bond energy strength                  */
double kth=20.0;                   /* angle energy strength                 */
double kph[3]={1.0,0,0.5};         /* torsion energy strength               */
double kcon=0.0;                   /* contact energy strength               */
double krep=1.0;                   /* excluded volume strength              */
double sigsa=4.0;                  /* bead diameter                         */
double cut=8.0;                    /* cutoff distance excluded volume       */
double ksi1=1.0;                   /* contact parameter                     */
double ksi2=25.0;                  /* contact parameter                     */
/* disordered regions */
double thn_disa=93.0;              /* reference bond angle                  */
double thn_disb=117.0;             /* reference bond angle                  */ 
double ksi_disa=4.0;               /* reference bond angle                  */
double ksi_disb=12.0;              /* reference bond angle                  */
double kph_dis[3]={0.5,0.8,0.25};  /* torsion energy strength               */
double phn_dis[3]={180,55,25};     /* reference torsion angle               */ 
/* crowders: */
const double rcrowd = 12.0;        /* Crowder radius                        */
const double srefcr = 3.0;         /* Crowder repulsion softness (sigma)    */
double epsilonrep = 1.0;           /* Crowder repulsion strength            */
double eclash = 1e6;               /* Crowder repulsion ceiling             */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/



