/************* simulation settings ******************************************/
# define NTMP 8                    /* # temperatures                        */
# define TMAX 1.00                 /* max temperature                       */
# define TMIN 0.80                 /* min temperature                       */
# define BOX 150                   /* simulation box                        */
# define ISTART 1                  /* 0 native, 1 read, 2 random            */
# define ISEED 1                   /* 1 randomize seed (/dev/urandom)       */
/************* MD parameters ************************************************/
# define MDSTEP (5000000000)       /* max # md steps                        */
# define NTHERM (1000)           /* # discarded steps                     */
# define IFLIP (100)               /* temperature flips                     */
# define ISAMP (100)               /* sample                                */
# define ICHECK (100000)           /* checkpoint                            */
# define IRT (10000)               /* runtime                               */
# define ICONF (100000)            /* configuration write                   */
# define IREADG 1                  /* read g parameters                     */
# define CHAIN_TO_BOX 1            /* translate chains/crowders periodically*/
                                   /* to original image in box [0...BOX]    */
/************* force field selection ****************************************/
# define FF_BOND 2                 /* bond() -- 1 on, 2 dual, 0 off         */
# define FF_BEND 2                 /* bend() -- 1 on, 2 dual, 0 off         */
# define FF_TORS 2                 /* tors() -- 1 on, 2 dual, 0 off         */
# define FF_CONT 2                 /* cont() -- 1 on, 2 dual, 0 off         */
# define FF_EXVOL 1                /* exvol()-- 1 on, 0 off                 */
# define FF_DISULF 2              /* disulfide bonds-- 1 on, 2 dual, 0 off */
# define FF_SEQ 0                  /* hp()   -- 1 on, 0 off                 */
# define FF_MULTIBODY 2            /* cont() -- multibody effects           */ 
/************* measurements *************************************************/
# define NBIN 200                  /* # bins                                */
# define NOBS 19                   /* # observables                         */
# define MAXCELL 25000             /* max # cells                           */
# define MAXP 1000                 /* max # contact pairs                   */
# define MAXNC 500                 /* max # contacts                        */
# define SNAP1 5000                /* write snapshots to directory SNAPDIR  */
# define SNAP2 5000                /* for interval SNAP1 < imd < SNAP2      */
# define RMSD 2                    /* 1 NATIVE, 2 NATIVE2, 0 off            */
/************* files input **************************************************/
# define NATIVE "native_2jp1_dimer_8-54_full"
# define NATIVE2 "native_1j8i_10-69_full"
# define CONTMAP "smog_2jp1_dimer_8-54_del14-18"
# define CONTMAP2 "smog_1j8i_10-69_mirror"
# define CONTMAP3 ""
# define CONTMAP4 ""
# define DISULFIDE "disulfide"
# define START "native_2jp1_dimer_8-54_full" //"native_1j8i_10-69_full"
# define INPUT "input"
# define INPUTG "inputg"
# define CONTPAR ""// "./cont_param_2KDL_smog"
# define CONTPAR2 ""// "./cont_param_2KDM_smog_0.92"
# define BONDEDPAR ""
# define BONDEDPAR2 ""
/************* files output *************************************************/
# define RT "rt"
# define PDB "current.pdb"
# define CONF "conf"
# define RAMA "rama"
# define OUTPUTG "outputg"
# define AVERAGES "averages"
# define HEATCAP "heat_capacity"
# define STATS "samp_stats"
# define RESDIR "results/"
# define ANADIR "analys/"
# define CHECKDIR "checkpnt/"
# define TESTDIR "results/param_check/"
# define SNAPDIR "snapshots/"
# define LOGFILE "logfile"
/************* functions ****************************************************/
# define max(A,B) ((A)>(B)?(A):(B))
# define min(A,B) ((A)<(B)?(A):(B))
# define sgn(A) ((A) < 0 ? -1 : 1)
/****************************************************************************/
