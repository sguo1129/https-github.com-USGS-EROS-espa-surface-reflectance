#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

# Define the script installation path
script_install_path = $(espa_project_dir)/bin
script_source_link_path = ../$(project_name)/bin

SCRIPTS = surface_reflectance.py

all:

install:
	install -d $(link_path)
	install -d $(script_install_path)
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(script_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(script_source_link_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done;

clean:

