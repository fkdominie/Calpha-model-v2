PROGRAM1 = fixtemp
PROGRAM2 = simtemp
PROGRAM3 = paramfiles
SCRIPT1 = constants
SCRIPT2 = pdb2cont
SCRIPT3 = analys_conf
SCRIPT4 = conf2pdb
SCRIPT5 = analys_conf_MBAR

FILES.c = energy.c geometry.c sampling.c obs.c misc.c utils.c
FILES.h = defs.h sys.h 
FILES.o = ${FILES.c:.c=.o}

CC      = gcc
SFLAGS  = -std=c11
GFLAGS  = 
OFLAGS  = -O3
WFLAG1  = -Wall
WFLAG2  = -Wextra
## WFLAG3  = -Werror
WFLAG4  = -Wstrict-prototypes
WFLAG5  = -Wmissing-prototypes
WFLAGS  = ${WFLAG1} ${WFLAG2} ${WFLAG3} ${WFLAG4} ${WFLAG5}
UFLAGS  = # Set on command line only

CFLAGS  = ${SFLAGS} ${GFLAGS} ${OFLAGS} ${WFLAGS} ${UFLAGS}
LDFLAGS =
LDLIBS  =


${PROGRAM1}: main_fixtemp.c ${FILES.c} ${FILES.o}
	${CC} -o main ${CFLAGS} main_fixtemp.c ${FILES.o} -lm 

${PROGRAM2}: main_simtemp.c ${FILES.c} ${FILES.o}
	${CC} -o main ${CFLAGS} main_simtemp.c ${FILES.o} -lm 

${PROGRAM3}: main_paramfiles.c ${FILES.c} ${FILES.o} 
	${CC} -o main_paramfiles ${CFLAGS} main_paramfiles.c ${FILES.o} -lm 

${SCRIPT1}: tools/constants.c ${FILES.o}
	${CC} -o constants ${CFLAGS} tools/constants.c -lm

${SCRIPT2}: tools/pdb2cont.c tools/readpdb3.o
	${CC} -o pdb2cont ${CFLAGS} tools/pdb2cont.c tools/readpdb3.o -lm 

${SCRIPT3}: tools/analys_conf.c tools/analys_conf.o ${FILES.c} ${FILES.o}
	${CC} -o analys_conf ${CFLAGS} tools/analys_conf.c  ${FILES.o} -lm 

${SCRIPT4}: tools/conf2pdb.c tools/conf2pdb.o ${FILES.c} ${FILES.o}
	${CC} -o conf2pdb ${CFLAGS} tools/conf2pdb.c  ${FILES.o} -lm 

${SCRIPT5}: tools/analys_conf_MBAR.c tools/analys_conf_MBAR.o ${FILES.c} ${FILES.o}
	${CC} -o analys_conf_MBAR ${CFLAGS} tools/analys_conf_MBAR.c  ${FILES.o} -lm 



energy.o: ${FILES.h} param.h
geometry.o: ${FILES.h}
sampling.o: ${FILES.h}
misc.o: ${FILES.h}
utils.o: ${FILES.h}
obs.o: ${FILES.h}
tools/readpdb3.o: 


clean:
	rm -fr ${FILES.o} main main_paramfiles pdb2cont tools/*.o

cleandata:
	rm -fr results/*

