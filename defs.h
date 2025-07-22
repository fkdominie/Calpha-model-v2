/************* simulation settings ******************************************/
# define NTMP 1                    /* # temperatures                        */
# define TMAX 0.92                 /* max temperature                       */
# define TMIN 0.92                 /* min temperature                       */
# define BOX 150                   /* simulation box                        */
# define ISTART 0                  /* 0 native, 1 read, 2 random            */
# define ISEED 1                   /* 1 randomize seed (/dev/urandom)       */
/************* MD parameters ************************************************/
# define MDSTEP (5000000000)       /* max # md steps                        */
# define NTHERM (10000)           /* # discarded steps                     */
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
# define FF_DISULF 2               /* disulfide bonds -- 1 on, 2 dual, 0 off*/
# define FF_MULTIBODY 1            /* multibody effects -- 1 on, 0 off      */
# define FF_SALT 0                 /* screening effect (csalt) -- 1 on 0 off*/
/************* measurements *************************************************/
# define NBIN 200                  /* # bins                                */
# define NOBS 27                   /* # observables                         */
# define MAXCELL 25000             /* max # cells                           */
# define MAXP 500                  /* max # contact pairs                   */
# define SNAP1 5000                /* write snapshots to directory SNAPDIR  */
# define SNAP2 5000                /* for interval SNAP1 < imd < SNAP2      */
# define RMSD 2                    /* 1 NATIVE, 2 NATIVE2, 0 off            */
/************* files input **************************************************/
# define NATIVE "native_1j8i_model1"
# define NATIVE2 "native_2jp1_fullchain_swap"
# define CONTMAP "smog_1j8i_r9-68_mirror_extra_20-24_9-46"
# define CONTMAP2 "smog_2jp1_r8-52_sym_intrachain"
# define DISULFIDE "disulfide"
# define START "native_1j8i_model1" 
# define INPUT "input"
# define INPUTG "inputg"
# define CONTPAR "./"
# define CONTPAR2 "./"
# define DISREG "./dis_regions_1j8i"
# define DISREG2 "./dis_regions_2jp1"
# define BONDEDPAR "./"
# define BONDEDPAR2 "./"
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
