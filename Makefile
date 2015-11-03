#-----------------------------------------------------------------------------
# Makefile
#
# Project Name: surface reflectance
#-----------------------------------------------------------------------------
.PHONY: check-environment all install clean all-script install-script clean-script all-ledaps install-ledaps clean-ledaps all-ledaps-aux install-ledaps-aux clean-ledaps-aux all-l8-sr install-l8-sr clean-l8-sr

include make.config

DIR_LEDAPS = ledaps
DIR_L8 = not-validated-prototype-l8_sr

#-----------------------------------------------------------------------------
all: all-ledaps all-l8-sr

install: check-environment install-ledaps install-l8-sr

clean: clean-ledaps clean-l8-sr

#-----------------------------------------------------------------------------
all-script:
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all);

install-script:
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install);

clean-script:
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean);

#-----------------------------------------------------------------------------
all-l8-sr:
	echo "make all in $(DIR_L8)"; \
        (cd $(DIR_L8); $(MAKE) all);

install-l8-sr:
	echo "make install in $(DIR_L8)"; \
        (cd $(DIR_L8); $(MAKE) install);

clean-l8-sr:
	echo "make clean in $(DIR_L8)"; \
        (cd $(DIR_L8); $(MAKE) clean);

#-----------------------------------------------------------------------------
all-ledaps:
	echo "make all in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) all-ledaps);

install-ledaps:
	echo "make install in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) install-ledaps);

clean-ledaps:
	echo "make clean in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) clean-ledaps);

#-----------------------------------------------------------------------------
all-ledaps-aux:
	echo "make all in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) all-ledaps-aux);

install-ledaps-aux:
	echo "make install in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) install-ledaps-aux);

clean-ledaps-aux:
	echo "make clean in $(DIR_LEDAPS)"; \
        (cd $(DIR_LEDAPS); $(MAKE) clean-ledaps-aux);

#-----------------------------------------------------------------------------
check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

