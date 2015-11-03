#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Ledaps.
#----------------------------------------------------------------------------
.PHONY: all install clean all-ledaps install-ledaps clean-ledaps all-ledaps-aux install-ledaps-aux clean-ledaps-aux

DIR_LEDAPS = ledapsSrc
DIR_LEDAPS_AUX = ledapsAncSrc

#-----------------------------------------------------------------------------
all: all-ledaps

install: install-ledaps

clean: clean-ledaps

#-----------------------------------------------------------------------------
all-ledaps:
	echo "make all in $(DIR_LEDAPS)..."; \
        (cd $(DIR_LEDAPS); $(MAKE))

install-ledaps: all-ledaps
	echo "make install in $(DIR_LEDAPS)..."; \
        (cd $(DIR_LEDAPS); $(MAKE) install)

clean-ledaps:
	echo "make clean in $(DIR_LEDAPS)..."; \
        (cd $(DIR_LEDAPS); $(MAKE) clean)

#-----------------------------------------------------------------------------
all-ledaps-aux:
	echo "make all in $(DIR_LEDAPS_AUX)..."; \
        (cd $(DIR_LEDAPS_AUX); $(MAKE))

install-ledaps-aux: all-ledaps-aux
	echo "make install in $(DIR_LEDAPS_AUX)..."; \
        (cd $(DIR_LEDAPS_AUX); $(MAKE) install)

clean-ledaps-aux:
	echo "make clean in $(DIR_LEDAPS_AUX)..."; \
        (cd $(DIR_LEDAPS_AUX); $(MAKE) clean)

