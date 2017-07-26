############################################################################
#
#    makefile: controls the build of various ParaDiS support utilities
#
#    Builds the following utilities:
#
#        paradisgen     --  problem generator
#        paradisrepart  --  creates a new problem decomposition with a new
#                           domain geometry from an previous nodal data file
#        paradisconvert --
#        ctablegen      --
#
############################################################################
include ../makefile.sys
include ../makefile.setup

DEFS += -D_THINFILM
#DEFS += -D_NOYOFFESTRESS
#DEFS += -D_NOVIRTUALSEG
DEFS += -D_TFIMGSTRESS
#DEFS += -D_BENDING
#DEFS += -D_STACKINGFAULT
# we need the following line to turn off multi node splitting for partials
#DEFS += -DMULTI_NODE_SPLIT_FREQ=100000

## The following flags for data output
#DEFS += -D_PRINTSTRESS
#DEFS += -D_WRITENODEFORCE

SRCDIR = ../src
INCDIR = ../include
BINDIR = ../bin

#
#	The utilities use various source modules from the parallel
#       code located in the parent directory.  Maintain a list of
#	these files.
#
#	These modules are compiled in the parent directory with a
#	different set of preprocessor definitions than are needed
#	here, so we need to create links in this directory back to
#	the source modlues and create separate object modules for
#	these sources.
#

EXTERN_C_SRCS = CellCharge.c   \
      Collision.c              \
      CommSendGhosts.c         \
      CommSendGhostPlanes.c    \
      CommSendMirrorNodes.c    \
      CommSendRemesh.c         \
      CommSendSecondaryGhosts.c\
      CommSendSegments.c       \
      CommSendVelocity.c       \
      CorrectionTable.c        \
      ParadisInit.c            \
      ParadisFinish.c          \
      DebugFunctions.c         \
      Decomp.c                 \
      DisableUnneededParams.c  \
      DLBfreeOld.c             \
      deWitInteraction.c       \
      FindPreciseGlidePlane.c  \
      FixRemesh.c              \
      FMComm.c                 \
      FMSigma2.c               \
      FMSupport.c              \
      FreeInitArrays.c         \
      GenerateOutput.c         \
      GetDensityDelta.c        \
      GetNewNativeNode.c       \
      GetNewGhostNode.c        \
      Gnuplot.c                \
      Heap.c                   \
      InitCellDomains.c        \
      InitCellNatives.c        \
      InitCellNeighbors.c      \
      InitHome.c               \
      InitRemoteDomains.c      \
      InitSendDomains.c        \
      LoadCurve.c              \
      Matrix.c                 \
      Meminfo.c                \
      MemCheck.c               \
      Migrate.c                \
      MobilityLaw_FCC_0b.c     \
      MobilityLaw_FCC_climb.c  \
      NodeVelocity.c           \
      OsmoticForce.c           \
      ParadisThread.c          \
      Parse.c                  \
      PickScrewGlidePlane.c    \
      QueueOps.c               \
      ReadRestart.c            \
      ReadBinaryRestart.c      \
      RBDecomp.c               \
      RSDecomp.c               \
      RemapInitialTags.c       \
      RemeshRule_3.c           \
      RemoteSegForces.c        \
      RemoveNode.c             \
      ResetGlidePlanes.c       \
      SortNativeNodes.c        \
      SortNodesForCollision.c  \
      SplitSurfaceNodes.c      \
      Tecplot.c                \
      Timer.c                  \
      WriteArms.c              \
      WriteAtomEye.c           \
      WriteBinaryRestart.c     \
      WriteDensFlux.c          \
      WriteDensityField.c      \
      WriteForce.c             \
      WriteFragments.c         \
      WritePoleFig.c           \
      WritePovray.c            \
      WriteProp.c              \
      WriteRestart.c           \
      WriteVelocity.c          \
      WriteVisit.c

EXTERN_CPP_SRCS = DisplayC.C   \
                  display.C 

EXTERN_SRCS = $(EXTERN_C_SRCS) $(EXTERN_CPP_SRCS)
EXTERN_OBJS = $(EXTERN_C_SRCS:.c=.o) $(EXTERN_CPP_SRCS:.C=.o)

#
#       Define the sources in the partial/ directory
#

# Files that need to be merged with partial/ directory
#                   Initialize.c          (1 block)
#                   LocalSegForces.c      (3 blocks) 
#                   NodeForce.c           (6 blocks) 
#                   Param.c               (1 block) 
#                   Topology.c            (RemoveDoubleLinks, SplitNode, MergeNode) 

#
#	Define the exectutable, source and object modules for
#	the problem generator.
#

PARADISTF     = paradistf
PARADISTF_BIN = $(BINDIR)/$(PARADISTF)

PARADISTF_C_SRCS = TF_Main.c                \
                   TF_Util.c                \
                   ABCcoeff.c               \
                   AllSegmentStress.c       \
                   AllYoffeStress.c         \
                   CrossSlipFCC.c           \
                   CrossSlipBCC.c           \
				   PredictiveCollision.c    \
				   ProximityCollision.c     \
                   Compute.c                \
          		   CrossSlip.c              \
                   DeltaPlasticStrain.c     \
                   DeltaPlasticStrain_BCC.c \
                   DeltaPlasticStrain_FCC.c \
                   DispStress.c             \
				   ForwardEulerIntegrator.c \
                   Fourier_transforms.c     \
                   Initialize.c             \
	    	       InputSanity.c            \
                   LocalSegForces.c         \
                   Minvmatrix.c             \
                   MobilityLaw_BCC_0.c      \
                   MobilityLaw_FCC_0.c      \
                   MobilityLaw_BCC_0b.c     \
                   MobilityLaw_BCC_glide.c  \
                   MobilityLaw_Relax.c      \
                   NodeForce.c              \
                   ParadisStep.c            \
                   Param.c                  \
                   Plot.c                   \
                   PrintStress.c            \
                   RemeshRule_2.c           \
                   ThinFilm_Remesh.c        \
                   thinfilm.c               \
                   Remesh.c                 \
                   SegmentStress.c          \
                   SemiInfiniteSegSegForce.c\
                   TrapezoidIntegrator.c    \
                   Topology.c               \
				   Util.c		 		    \
                   Yoffe.c                  \
                   Yoffe_corr.c             

PARADISTF_INCS.linux = -I/usr/local/include -I/home/ryu/gsl/include
PARADISTF_LIBS.linux = -L/usr/local/lib/ -lfftw3 -lm 

#PARADISTF_INCS.linux = -I/opt/apps/intel17/impi17_0/fftw3/3.3.6/include -I/opt/apps/intel17/gsl/2.3/include
#PARADISTF_LIBS.linux = -L/opt/apps/intel17/impi17_0/fftw3/3.3.6/lib -lfftw3 -L/opt/apps/intel17/gsl/2.3/lib -lgsl -lgslcblas -lm 

#PARADISTF_INCS.linux = -I/usr/include
#PARADISTF_LIBS.linux = -L/usr/lib/ -lfftw3 -lm

PARADISTF_INCS.COMET = -I/opt/fftw/3.3.4/intel/openmpi_ib/include -I/opt/gsl/2.1/intel/include
PARADISTF_LIBS.COMET = -L/opt/fftw/3.3.4/intel/openmpi_ib/lib -lfftw3 -L/opt/gsl/2.1/intel/lib -lgsl -lgslcblas -lm 

PARADISTF_INCS.STAMPEDE2 = -I/opt/apps/intel17/impi17_0/fftw3/3.3.6/include -I/opt/apps/intel17/gsl/2.3/include
PARADISTF_LIBS.STAMPEDE2 = -L/opt/apps/intel17/impi17_0/fftw3/3.3.6/lib -lfftw3 -L/opt/apps/intel17/gsl/2.3/lib -lgsl -lgslcblas -lm 

PARADISTF_INCS.OfficeServer = -I/usr/local/include -I/home/ryu/gsl/include
PARADISTF_LIBS.OfficeServer = -L/usr/local/lib/ -lfftw3 -L/home/ryu/gsl/lib -lgsl -lgslcblas -L/opt/intel/compilers_and_libraries/linux/lib/intel64


PARADISTF_INCS.mc2 = -I/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/include
PARADISTF_LIBS.mc2 = -L/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/lib/ -lfftw3 -lm

PARADISTF_INCS.linux.pc = -I/usr/include
PARADISTF_LIBS.linux.pc = -L/usr/lib/ -lfftw3 -lm


PARADISTF_INCS.wcr = -I/opt/fftw-3.1.2/intel/include
PARADISTF_LIBS.wcr = -L/opt/fftw-3.1.2/intel/lib/ -lfftw3  -lm

PARADISTF_INCS.vip = 
PARADISTF_LIBS.vip =  -L/u/system/Power/libs/fftw-3.1.2/lib -lfftw3 -lm

PARADISTF_INCS.su-ahpcrc = -I/lustre/home/mpotts/fftw/include
PARADISTF_LIBS.su-ahpcrc = -L/lustre/home/mpotts/fftw/lib -lfftw3 -lm  

# For mjm
PARADISTF_INCS.linux.opteron = -I/usr/cta/pet/MATH/include
PARADISTF_LIBS.linux.opteron = -L/usr/cta/pet/MATH/lib/ -lfftw3 -lm

PARADISTF_INCS.cygwin = 
PARADISTF_LIBS.cygwin = -lfftw3

PARADISTF_INCS.mc-cc = 
PARADISTF_LIBS.mc-cc = -lfftw3

PARADISTF_INCS.cray =-I/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/include
PARADISTF_LIBS.cray =-L/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/lib/ -lfftw3


PARADISTF_INCS.Harold = -I/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/include
PARADISTF_LIBS.Harold = -L/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/lib -lfftw3


PARADISTF_INCS = -I Include $(PARADISTF_INCS.$(SYS))
PARADISTF_LIBS = $(PARADISTF_LIBS.$(SYS))

PARADISTF_SRCS = $(PARADISTF_C_SRCS) $(PARADISTF_CPP_SRCS)
PARADISTF_OBJS = $(PARADISTF_C_SRCS:.c=.o) $(PARADISTF_CPP_SRCS:.C=.o)



###########################################################################
#
#	Define a rule for converting .c files to object modules.
#	All modules are compile serially in this directory
#
###########################################################################

.c.o:		makefile ../makefile.sys ../makefile.setup
		$(CC) $(OPT) $(CCFLAG) $(PARADISTF_INCS) $(INCS) -c $<

.C.o:		makefile ../makefile.sys ../makefile.setup
		$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $<


###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(EXTERN_OBJS) $(PARADISTF) 

clean:
		rm -f *.o $(EXTERN_SRCS) $(PARADISTF_BIN) test_YoffeInfStress

depend:		 *.c $(SRCDIR)/*.c $(INCDIR)/*.h makefile
		makedepend -Y$(INCDIR) *.c  -fmakefile.dep

#
#	Create any necessary links in the current directory to source
#	modules located in the SRCDIR directory
#

$(EXTERN_SRCS): $(SRCDIR)/$@
		- @ ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1

# For vip
#$(EXTERN_SRCS): $(SRCDIR)/$@
#                ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1


$(PARADISTF):	$(PARADISTF_BIN)
$(PARADISTF_BIN): $(PARADISTF_SRCS) $(PARADISTF_OBJS) $(EXTERN_OBJS) $(HEADERS)
	echo $(PARADISTF_OBJS)
		$(CPP) $(OPT) $(PARADISTF_OBJS) $(EXTERN_OBJS) -o $@  $(LIB) $(PARADISTF_LIBS)

test_YoffeInfStress: test_YoffeInfStress.o AllYoffeStress.o Yoffe_corr.o Yoffe.o Util.o Heap.o QueueOps.o
	$(CC) -o $@ $^ -lm

