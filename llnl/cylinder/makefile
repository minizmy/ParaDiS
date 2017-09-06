#############################################################################
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

DEFS += -D_CYLINDER
DEFS += -D_NOYOFFESTRESS
DEFS += -D_NOVIRTUALSEG
#DEFS += -D_CYLIMGSTRESS
#DEFS += -D_CYLMETHOD1
#DEFS += -D_STACKINGFAULT
#DEFS += -D_TORSION
#DEFS += -D_TENSIONafterTORSION
#DEFS += -D_TORSIONDEBUG

## The following flags for data output
#DEFS += -D_PRINTSTRESS
#DEFS += -D_WRITENODEFORCE

## The following flags for test case for cylinder code 
#DEFS += -D_CYL_TEST1
#DEFS += -D_CYL_TEST23

## The following flags for surface dislocation loop nucleation
DEFS += -D_NUCLEATION

## The following flags for boundary
#DEFS += -D_BOUNDARY

## The following flags for RELAX step to stop 
## when the density less than this.
## Not workiing very well since it is not clear that there is 
## a critical plastic strain when the structure is stabilized.
#DEFS += -D_RELAX

## The following flags for BCC surface cross slip
#DEFS += -D_BCC_CROSS

## The following flags to write essential output only
DEFS += -D_WRITEPROP_SHORT

## The following flags for nonproportional loading : Surface Coating
#DEFS += -D_SurfCoat

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
      Decomp.c                 \
      DisableUnneededParams.c  \
      DLBfreeOld.c             \
      deWitInteraction.c       \
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
      Matrix.c                 \
      Meminfo.c                \
      MemCheck.c               \
      Migrate.c                \
      MobilityLaw_FCC_0b.c     \
      MobilityLaw_FCC_climb.c  \
      MobilityLaw_BCC_0b.c     \
      MobilityLaw_Relax.c      \
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
      WriteBinaryRestart.c     \
      WriteDensFlux.c          \
      WriteDensityField.c      \
      WriteForce.c             \
      WriteFragments.c         \
      WritePoleFig.c           \
      WritePovray.c            \
      WriteRestart.c           \
      WriteVelocity.c          \
      WriteVisit.c

EXTERN_CPP_SRCS = DisplayC.C   \
                  display.C 

EXTERN_SRCS = $(EXTERN_C_SRCS) $(EXTERN_CPP_SRCS)
EXTERN_OBJS = $(EXTERN_C_SRCS:.c=.o) $(EXTERN_CPP_SRCS:.C=.o)

# Files that need to be merged with partial/ directory
#                   Initialize.c          (1 block)
#                   LocalSegForces.c      (3 blocks) 
#                   NodeForce.c           (6 blocks) 

#
#	Define the exectutable, source and object modules for
#	the problem generator.
#

PARADISCYL     = paradiscyl
PARADISCYL_BIN = $(BINDIR)/$(PARADISCYL)
PARADISCYL_C_SRCS= CYL_Main.c               \
                   CYL_Util.c               \
                   CYL_Funcs.c              \
                   ABCcoeff.c               \
                   AllSegmentStress.c       \
                   AllYoffeStress.c         \
                   CrossSlipFCC.c           \
                   CrossSlipBCC.c           \
	      		   DebugFunctions.c         \
				   PredictiveCollision.c    \
				   ProximityCollision.c     \
				   FindPreciseGlidePlane.c  \
				   Collision.c              \
                   Compute.c                \
          		   CrossSlip.c              \
                   DeltaPlasticStrain.c     \
                   DeltaPlasticStrain_BCC.c \
                   DeltaPlasticStrain_FCC.c \
	        	   ForwardEulerIntegrator.c \
                   Fourier_transforms.c     \
                   Initialize.c             \
	    	       InputSanity.c            \
				   LoadCurve.c              \
                   LocalSegForces.c         \
                   LoopGenerate.c	    	\
                   LoopGenerate_BCC.c	    \
                   LoopGenerate_FCC.c	    \
				   MobilityLaw_BCC_0.c      \
                   MobilityLaw_FCC_0.c      \
                   MobilityLaw_BCC_glide.c  \
                   NodeForce.c              \
                   ParadisStep.c            \
                   Param.c                  \
                   Plot.c                   \
                   PrintStress.c            \
                   RemeshRule_2.c           \
                   Cylinder_Remesh.c        \
                   cylinder.c               \
                   Remesh.c                 \
                   SegmentStress.c          \
                   SemiInfiniteSegSegForce.c\
                   TrapezoidIntegrator.c    \
                   Topology.c               \
				   Boundary.c   		    \
				   WriteAtomEye.c           \
				   WriteProp.c              \
				   Util.c		 		    \
                   Yoffe.c                  \
                   Yoffe_corr.c             \
                   final_matrices.c         \
                   greenstress_a.c          \
                   gridstress.c             \
                   m2matrix.c               \
                   n2invmatrix.c            \
                   n2matrix.c               \
                   m2invmatrix.c            \
                   minvmatrix.c             \
                   mmatrix.c                \
                   ninvmatrix.c             \
                   nmatrix.c                  

PARADISCYL_INCS.linux = -I/usr/local/include -I ~/gsl/include
PARADISCYL_LIBS.linux = -L/usr/local/lib/ -lfftw3 -L ~/gsl/lib -lgsl -lgslcblas -L/opt/intel/compilers_and_libraries/linux/lib/intel64

PARADISCYL_INCS.OfficeServer = -I/usr/local/include -I ~/gsl/include
PARADISCYL_LIBS.OfficeServer = -L/usr/local/lib/ -lfftw3 -L ~/gsl/lib -lgsl -lgslcblas -L/opt/intel/compilers_and_libraries/linux/lib/intel64

PARADISCYL_INCS.wcr = -I/opt/fftw-3.1.2/intel/include
PARADISCYL_LIBS.wcr = -L/opt/fftw-3.1.2/intel/lib/ -lfftw3  -lg2c -lm

PARADISCYL_INCS.mc2 = -I/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/include
PARADISCYL_LIBS.mc2 = -L/export/apps/fftw3/3.3alpha1/mvapich/1.2rc1-intel-11-dell-gen2/lib/ -lfftw3 -L$(HOME)/usr/lib -lgsl -lgslcblas -lm

PARADISCYL_INCS.vip = 
PARADISCYL_LIBS.vip =  -L/u/system/Power/libs/fftw-3.1.2/lib -lfftw3 -lm
PARADISCYL_LIBS.vip =  -L/u/system/Power/libs/fftw-3.1.2/lib -lfftw3 -lm

PARADISCYL_INCS.su-ahpcrc = -I/lustre/home/mpotts/fftw/include
PARADISCYL_LIBS.su-ahpcrc = -L/lustre/home/mpotts/fftw/lib -lfftw3 $(HOME)/usr/lib/libgsl.a -lm  

PARADISCYL_INCS.cygwin = 
PARADISCYL_LIBS.cygwin = -lfftw3

PARADISCYL_INCS.mc-cc = 
PARADISCYL_LIBS.mc-cc = -lfftw3

PARADISCYL_INCS = -I Include $(PARADISCYL_INCS.$(SYS))
PARADISCYL_LIBS = $(PARADISCYL_LIBS.$(SYS))

PARADISCYL_SRCS = $(PARADISCYL_C_SRCS) $(PARADISCYL_CPP_SRCS)
PARADISCYL_OBJS = $(PARADISCYL_C_SRCS:.c=.o) $(PARADISCYL_CPP_SRCS:.C=.o)



###########################################################################
#
#	Define a rule for converting .c files to object modules.
#	All modules are compile serially in this directory
#
###########################################################################

.c.o:		makefile ../makefile.sys ../makefile.setup
		$(CC) $(OPT) $(CCFLAG) $(PARADISCYL_INCS) $(INCS) -c $<

.C.o:		makefile ../makefile.sys ../makefile.setup
		$(CPP) $(OPT) $(CPPFLAG) $(INCS) -c $<


###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(EXTERN_OBJS) $(PARADISCYL) 

clean:
		rm -f *.o $(EXTERN_SRCS) $(PARADISCYL_BIN)

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

$(PARADISCYL):	$(PARADISCYL_BIN)
$(PARADISCYL_BIN): $(PARADISCYL_SRCS) $(PARADISCYL_OBJS) $(EXTERN_OBJS) $(HEADERS)
		$(CPP) $(OPT) $(PARADISCYL_OBJS) $(EXTERN_OBJS) -o $@  $(LIB) $(PARADISCYL_LIBS)

