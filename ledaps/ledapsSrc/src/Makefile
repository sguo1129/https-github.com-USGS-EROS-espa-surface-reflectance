#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Ledaps.
#-----------------------------------------------------------------------------
.PHONY: all install clean

MODULES = lndpm lndcal lndsr lndsrbm 6sV-1.0B

all:
	@for module in $(MODULES); do \
	    (cd $$module; $(MAKE)) \
	done

install:
	@for module in $(MODULES); do \
	    (cd $$module; $(MAKE) install) \
	done

clean:
	@for module in $(MODULES); do \
	    (cd $$module; $(MAKE) clean) \
	done
