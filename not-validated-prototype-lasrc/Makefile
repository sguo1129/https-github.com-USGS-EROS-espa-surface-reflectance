#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing LaSRC.
#----------------------------------------------------------------------------
.PHONY: all install clean all-lasrc install-lasrc clean-lasrc all-lasrc-aux install-lasrc-aux clean-lasrc-aux

DIR_LaSRC = c_version
DIR_LaSRC_AUX = landsat_aux

#-----------------------------------------------------------------------------
all: all-lasrc

install: install-lasrc

clean: clean-lasrc

#-----------------------------------------------------------------------------
all-lasrc:
	echo "make all in $(DIR_LaSRC)..."; \
        (cd $(DIR_LaSRC); $(MAKE))

install-lasrc: all-lasrc
	echo "make install in $(DIR_LaSRC)..."; \
        (cd $(DIR_LaSRC); $(MAKE) install)

clean-lasrc:
	echo "make clean in $(DIR_LaSRC)..."; \
        (cd $(DIR_LaSRC); $(MAKE) clean)

#-----------------------------------------------------------------------------
all-lasrc-aux:
	echo "make all in $(DIR_LaSRC_AUX)..."; \
        (cd $(DIR_LaSRC_AUX); $(MAKE))

install-lasrc-aux: all-lasrc-aux
	echo "make install in $(DIR_LaSRC_AUX)..."; \
        (cd $(DIR_LaSRC_AUX); $(MAKE) install)

clean-lasrc-aux:
	echo "make clean in $(DIR_LaSRC_AUX)..."; \
        (cd $(DIR_LaSRC_AUX); $(MAKE) clean)

