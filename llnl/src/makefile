
#############################################################################
#
#    makefile:  Controls the build of the primary ParaDiS executable
#               and the parallel FMM correction table generator
#
#    Usage:
#        gmake           build dd3d executable and some of the
#                        associated tools and utilities.
#        gmake clean     remove executable and object files
#        gmake depend    update makefile dependencies
#
############################################################################

#
#       Include files containing the system specific macro definitions
#       as well as the settings and flags that are not system dependent.
#       Also include the source and object module lists for the executables
#       to be compiled
#

include ../makefile.sys
include ../makefile.setup
include ../makefile.srcs

INCDIR = ../include
BINDIR = ../bin

#
#	Define the executable for the application
#

PARADIS = paradis
PARADIS_BIN = $(BINDIR)/$(PARADIS)

#
#       Define the exectutable for the parallel version of the FMM
#       image correction table generator
#

CTABLEGENP = ctablegenp
CTABLEGENP_BIN = $(BINDIR)/$(CTABLEGENP)


###########################################################################
#
#	Define all the suffix rules needed to compile the
#	various modules
#
###########################################################################

.SUFFIXES : .C makefile

.c.o:		makefile
		$(CC) $(OPT) $(CCFLAG) $(INCS) -c $<

.cpp.o:		makefile
		$(CPP) $(OPT) $(CCFLAG) $(CPPFLAG) $(INCS) -c $<

.C.o:		makefile
		$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $<



###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(PARADIS) $(CTABLEGENP)

clean:
		rm -f *.o $(PARADIS_BIN) $(CTABLEGENP_BIN) 


depend:		*.c *.C $(INCDIR)/*h makefile ../makefile.setup
		makedepend -Y$(INCDIR) *.c *.C -fmakefile.dep

$(BINDIR):
		mkdir $(BINDIR)

#
# Temporarily add a target for CommSendGhosts.o. There is
# currently a bug in the compiler on uBGL which causes
# badness in ghost comm stuff if any optimizations are
# used.  So, compile this module without any optimizations
# until the compiler is fixed.
#

CommSendGhosts.o:	CommSendGhosts.c $(HEADERS)
		$(CC) -g $(CCFLAG) $(INCS) -c $<


$(PARADIS):	$(BINDIR) $(PARADIS_BIN)
$(PARADIS_BIN):	Main.o $(PARADIS_OBJS) $(HEADERS)
		$(CPP) Main.o $(OPT) $(OPENMP_FLAG) $(PARADIS_OBJS) \
		-o $@  $(LIB)


$(CTABLEGENP):	$(BINDIR) $(CTABLEGENP_BIN)
$(CTABLEGENP_BIN): 	$(CTABLEGEN_SRCS) $(CTABLEGEN_OBJS)
		$(CPP) $(OPT) $(OPENMP_FLAG) $(CTABLEGEN_OBJS) -o $@  $(LIB)

purify:		Main.o $(PARADIS_OBJS) $(HEADERS)
		purify ./gcc Main.o $(OPT) $(PARADIS_OBJS) -o $(PARADIS) \
			$(LIB) 

prof:		Main.o $(PARADIS_OBJS) $(HEADERS)
		$(CPP) Main.o $(OPT) $(PARADIS_OBJS) -p  -o $(PARADIS) $(LIB) 


include ./makefile.dep
