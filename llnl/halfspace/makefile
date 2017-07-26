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

DEFS += -D_HALFSPACE
#DEFS += -D_NOYOFFESTRESS
#DEFS += -D_NOVIRTUALSEG
DEFS += -D_HSIMGSTRESS
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

EXTERN_C_SRCS = CellCharge.c             \
      Collision.c              \
      CommSendGhosts.c         \
      CommSendGhostPlanes.c    \
      CommSendMirrorNodes.c    \
      CommSendRemesh.c         \
      CommSendSecondaryGhosts.c \
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
      ReadRestart.c	       \
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

EXTERN_CPP_SRCS = DisplayC.C       \
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

PARADISHS     = paradishs
PARADISHS_BIN = $(BINDIR)/$(PARADISHS)

PARADISHS_C_SRCS = HS_Main.c                \
                   HS_Util.c                \
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
                   DeltaPlasticStrain_BCC.c     \
                   DeltaPlasticStrain_FCC.c     \
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
                   HalfSpace_Remesh.c        \
                   halfspace.c               \
                   Remesh.c                 \
                   SegmentStress.c          \
                   SemiInfiniteSegSegForce.c          \
                   TrapezoidIntegrator.c    \
                   Topology.c               \
		   Util.c		    \
                   Yoffe.c                  \
                   Yoffe_corr.c             

PARADISHS_INCS.linux = -I/usr/include
PARADISHS_LIBS.linux = -L/usr/lib/ -lfftw3 -lm

PARADISHS_INCS.mc2 = -I/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/include
PARADISHS_LIBS.mc2 = -L/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/lib/ -lfftw3 -lm

PARADISHS_INCS.linux.pc = -I/usr/include
PARADISHS_LIBS.linux.pc = -L/usr/lib/ -lfftw3 -lm


PARADISHS_INCS.wcr = -I/opt/fftw-3.1.2/intel/include
PARADISHS_LIBS.wcr = -L/opt/fftw-3.1.2/intel/lib/ -lfftw3  -lm

PARADISHS_INCS.vip = 
PARADISHS_LIBS.vip =  -L/u/system/Power/libs/fftw-3.1.2/lib -lfftw3 -lm

PARADISHS_INCS.su-ahpcrc = -I/lustre/home/mpotts/fftw/include
PARADISHS_LIBS.su-ahpcrc = -L/lustre/home/mpotts/fftw/lib -lfftw3 -lm  

# For mjm
PARADISHS_INCS.linux.opteron = -I/usr/cta/pet/MATH/include
PARADISHS_LIBS.linux.opteron = -L/usr/cta/pet/MATH/lib/ -lfftw3 -lm

PARADISHS_INCS.cygwin = 
PARADISHS_LIBS.cygwin = -lfftw3

PARADISHS_INCS.mc-cc = 
PARADISHS_LIBS.mc-cc = -lfftw3

PARADISHS_INCS.cray =-I/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/include
PARADISHS_LIBS.cray =-L/mnt/cfs/pkgs/PTOOLS/pkgs/fftw-2.1.5/lib/ -lfftw3


PARADISHS_INCS.Harold = -I/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/include
PARADISHS_LIBS.Harold = -L/mnt/lustre/usrcta/pet/pkgs/fftw/3.2.2/lib -lfftw3


PARADISHS_INCS = -I Include $(PARADISHS_INCS.$(SYS))
PARADISHS_LIBS = $(PARADISHS_LIBS.$(SYS))

PARADISHS_SRCS = $(PARADISHS_C_SRCS) $(PARADISHS_CPP_SRCS)
PARADISHS_OBJS = $(PARADISHS_C_SRCS:.c=.o) $(PARADISHS_CPP_SRCS:.C=.o)



###########################################################################
#
#	Define a rule for converting .c files to object modules.
#	All modules are compile serially in this directory
#
###########################################################################

.c.o:		makefile ../makefile.sys ../makefile.setup
		$(CC) $(OPT) $(CCFLAG) $(PARADISHS_INCS) $(INCS) -c $<

.C.o:		makefile ../makefile.sys ../makefile.setup
		$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $<


###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(EXTERN_OBJS) $(PARADISHS) 

clean:
		rm -f *.o $(EXTERN_SRCS) $(PARADISHS_BIN) test_YoffeInfStress

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


$(PARADISHS):	$(PARADISHS_BIN)
$(PARADISHS_BIN): $(PARADISHS_SRCS) $(PARADISHS_OBJS) $(EXTERN_OBJS) $(HEADERS)
	echo $(PARADISHS_OBJS)
		$(CPP) $(OPT) $(PARADISHS_OBJS) $(EXTERN_OBJS) -o $@  $(LIB) $(PARADISHS_LIBS)

test_YoffeInfStress: test_YoffeInfStress.o AllYoffeStress.o Yoffe_corr.o Yoffe.o Util.o Heap.o QueueOps.o
	$(CC) -o $@ $^ -lm

