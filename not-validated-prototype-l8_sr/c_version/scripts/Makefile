#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing L8 SR.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

# Define scripts
SCRIPTS = do_l8_sr.py

#-----------------------------------------------------------------------------
all:

install:
	install -d $(link_path)
	install -d $(l8_bin_install_path)
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(l8_bin_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(l8_link_source_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done;

clean:

