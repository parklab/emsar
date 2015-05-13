COMPILER = gcc
DFLAGS = 
CFLAGS = -c
OFLAGS = -lpthread -lm -lz -o
EXECNAME = emsar
EXECNAME2 = emsar-build

all:   emsar emsar-build

emsar:   bam.o sam.o alignment.o stringhash.o bool.o emsar_main.o emsar_functions.o 
	${COMPILER} ${DFLAGS} bam.o sam.o alignment.o stringhash.o bool.o emsar_main.o emsar_functions.o ${OFLAGS} ${EXECNAME}

emsar-build:   bam.o sam.o alignment.o stringhash.o bool.o emsar_build_main.o emsar_functions.o 
	${COMPILER} ${DFLAGS} bam.o sam.o alignment.o stringhash.o bool.o emsar_build_main.o emsar_functions.o ${OFLAGS} ${EXECNAME2}

emsar_functions.o:   emsar_functions.c emsar.h
	${COMPILER} ${DFLAGS} ${CFLAGS} emsar_functions.c

emsar_main.o:     emsar_main.c emsar.h
	${COMPILER} ${DFLAGS} ${CFLAGS} emsar_main.c

emsar_build_main.o:     emsar_build_main.c emsar.h
	${COMPILER} ${DFLAGS} ${CFLAGS} emsar_build_main.c

bool.o:   bool.c bool.h
	${COMPILER} ${DFLAGS} ${CFLAGS} bool.c

stringhash.o:   stringhash.c stringhash.h
	${COMPILER} ${DFLAGS} ${CFLAGS} stringhash.c

alignment.o:   alignment.c alignment.h
	${COMPILER} ${DFLAGS} ${CFLAGS} alignment.c

sam.o:   sam.c sam.h faidx.c faidx.h razf.c razf.h
	${COMPILER} ${DFLAGS} ${CFLAGS} -combine sam.c faidx.c razf.c

bam.o:   bam.c bam.h bam_endian.h kstring.c kstring.h sam_header.c sam_header.h khash.h bgzf.c bgzf.h bam_import.c kseq.h bam_aux.c
	${COMPILER} ${DFLAGS} ${CFLAGS} -combine bam.c kstring.c sam_header.c bgzf.c bam_import.c bam_aux.c



