
#############################################################################
#
#    makefile:  This is the primary make file controlling the build
#               of the ParaDiS parallel code and associated utilities.
#
#    Usage:
#        gmake           build paradis executable and some of the
#                        associated tools and utilities.
#        gmake clean     remove executable and object files
#        gmake depend    update makefile dependencies
#
############################################################################

#
#       Include files containing the system specific macro definitions
#       as well as the settings and flags that are not system dependent.
#

include makefile.sys
include makefile.setup


###########################################################################
#
#	Define the object modules for the application
#
###########################################################################

PARADIS    = paradis
CTABLEGENP = ctablegenp
UTILS      = utils

DIRS    = src utilities

BINDIR = ./bin

###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(PARADIS) $(CTABLEGENP) $(UTILS)

clean:
		@ ( for i in $(DIRS) ; do	\
			cd $$i ;		\
			$(MAKE) $@ ;		\
			cd .. ;			\
		done ;				\
		)

$(BINDIR):
		mkdir $(BINDIR)

$(PARADIS):	$(BINDIR)
		@ ( cd src ;             \
		    $(MAKE) $@ ;           \
		    cd .. ;              \
		  )

$(CTABLEGENP): 	$(BINDIR)
		@ ( cd src ;             \
		    $(MAKE) $@ ;           \
		    cd .. ;              \
		  )

$(UTILS):	$(BINDIR)
		@ ( cd utilities ;       \
		    $(MAKE) ;              \
		    cd .. ;              \
		  )

depend: 	
		@ ( for i in $(DIRS) ; do	\
			cd $$i ;		\
			$(MAKE) $@ ;		\
			cd .. ;			\
		done ;				\
		)

purify:		$(BINDIR)	
		@ ( cd src ;             \
		    $(MAKE) $@ ;           \
		    cd .. ;              \
		  )

prof:		$(BINDIR)
		@ ( cd src ;             \
		    $(MAKE) $@ ;           \
		    cd .. ;              \
		  )

