#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing L8 SR.
#----------------------------------------------------------------------------
.PHONY: all install clean all-l8-sr install-l8-sr clean-l8-sr all-l8-sr-aux install-l8-sr-aux clean-l8-sr-aux

DIR_L8_SR = c_version
DIR_L8_SR_AUX = landsat_aux

#-----------------------------------------------------------------------------
all: all-l8-sr

install: install-l8-sr

clean: clean-l8-sr

#-----------------------------------------------------------------------------
all-l8-sr:
	echo "make all in $(DIR_L8_SR)..."; \
        (cd $(DIR_L8_SR); $(MAKE))

install-l8-sr: all-l8-sr
	echo "make install in $(DIR_L8_SR)..."; \
        (cd $(DIR_L8_SR); $(MAKE) install)

clean-l8-sr:
	echo "make clean in $(DIR_L8_SR)..."; \
        (cd $(DIR_L8_SR); $(MAKE) clean)

#-----------------------------------------------------------------------------
all-l8-sr-aux:
	echo "make all in $(DIR_L8_SR_AUX)..."; \
        (cd $(DIR_L8_SR_AUX); $(MAKE))

install-l8-sr-aux: all-l8-sr-aux
	echo "make install in $(DIR_L8_SR_AUX)..."; \
        (cd $(DIR_L8_SR_AUX); $(MAKE) install)

clean-l8-sr-aux:
	echo "make clean in $(DIR_L8_SR_AUX)..."; \
        (cd $(DIR_L8_SR_AUX); $(MAKE) clean)

