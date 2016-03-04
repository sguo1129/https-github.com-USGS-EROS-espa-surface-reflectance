#-----------------------------------------------------------------------------
# Makefile
#
# For building lndcal.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the source code and object files
SRC1 = lndsrbm.c
OBJ1 = $(SRC1:.c=.o)
SRC2 = comptemp.c
OBJ2 = $(SRC2:.c=.o)
SRC3 = dump_meta.c
OBJ3 = $(SRC3:.c=.o)

SRC4 = LS_geoloc_driver.c
OBJ4 = LS_geoloc_driver_inv.o
SRC5 = LS_geoloc_driver.c
OBJ5 = LS_geoloc_driver_for.o
SRC6 = LS_geoloc.c
OBJ6 = LS_geoloc_inv.o
SRC7 = LS_geoloc.c
OBJ7 = LS_geoloc_for.o

SRC8 = SDSreader3.0.c
OBJ8 = $(SRC8:.c=.o)

ALL_OBJ = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8)

# Define include paths
INCDIR  = -I. -I$(ESPAINC) -I$(XML2INC)
SDS_INCDIR = -I$(HDFINC)
GEOLOC_INCDIR = -I$(HDFEOS_GCTPINC)
NCFLAGS = $(EXTRA) $(INCDIR) $(SDS_INCDIR) $(GEOLOC_INCDIR)

# Define the object libraries and paths
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz
SDS_EXLIB = -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg
GEOLOC_EXLIB  = -L$(HDFEOS_GCTPLIB) -lGctp
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define scripts
SCRIPTS = lndsrbm.ksh

# Define C executables
EXE1 = lndsrbm
EXE2 = comptemp
EXE3 = dump_meta
EXE6 = xy2geo
EXE7 = geo2xy
EXE8 = SDSreader3.0
ALL_EXE = $(EXE1) $(EXE2) $(EXE3) $(EXE6) $(EXE7) $(EXE8)

#-----------------------------------------------------------------------------
all: $(ALL_OBJ) $(ALL_EXE)

$(EXE1): $(OBJ1)
	$(CC) $(EXTRA) -o $(EXE1) $(OBJ1) $(GEOLOC_EXLIB) $(LOADLIB)

$(EXE2): $(OBJ2)
	$(CC) $(EXTRA) -o $(EXE2) $(OBJ2) $(GEOLOC_EXLIB) $(LOADLIB)

$(EXE3): $(OBJ3)
	$(CC) $(EXTRA) -o $(EXE3) $(OBJ3) $(GEOLOC_EXLIB) $(LOADLIB)

$(EXE6): $(OBJ4) $(OBJ6)
	$(CC) $(EXTRA) -o $(EXE6) $(OBJ4) $(OBJ6) $(GEOLOC_EXLIB) $(LOADLIB)

$(EXE7): $(OBJ5) $(OBJ7)
	$(CC) $(EXTRA) -o $(EXE7) $(OBJ5) $(OBJ7) $(GEOLOC_EXLIB) $(LOADLIB)

$(EXE8): $(OBJ8)
	$(CC) $(EXTRA) -o $(EXE8) $(OBJ8) $(SDS_EXLIB) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	@for executable in $(ALL_EXE); do \
            cmd="install -m 755 $$executable $(ledaps_bin_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(ledaps_link_source_path)/$$executable $(link_path)/$$executable"; \
            echo "$$cmd"; $$cmd; \
        done
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(ledaps_bin_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(ledaps_link_source_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done

#-----------------------------------------------------------------------------
clean:
	$(RM) -f *.o $(ALL_EXE)

#-----------------------------------------------------------------------------
$(OBJ1): $(SRC1)
$(OBJ2): $(SRC2)
$(OBJ3): $(SRC3)
$(OBJ8): $(SRC8)
.c.o:
	$(CC) $(NCFLAGS) -c $<

$(OBJ4): $(SRC4)
	$(CC) $(NCFLAGS) -DINV -c $(SRC4) -o $(OBJ4)
$(OBJ5): $(SRC5)
	$(CC) $(NCFLAGS) -DFOR -c $(SRC5) -o $(OBJ5)

$(OBJ6): $(SRC6)
	$(CC) $(NCFLAGS) -DINV -c $(SRC6) -o $(OBJ6)
$(OBJ7): $(SRC7)
	$(CC) $(NCFLAGS) -DFOR -c $(SRC7) -o $(OBJ7)

