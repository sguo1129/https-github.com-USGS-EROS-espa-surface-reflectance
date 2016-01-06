#-----------------------------------------------------------------------------
# Makefile for L8 SR code
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

#-----------------------------------------------------------------------------
# Set up compile options
CC    = gcc
RM    = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the include files
INC = common.h date.h input.h output.h lut_subr.h l8_sr.h

# Define the source code and object files
SRC = compute_refl.c      \
      date.c              \
      get_args.c          \
      input.c             \
      lut_subr.c          \
      output.c            \
      subaeroret.c        \
      l8_sr.c
OBJ = $(SRC:.c=.o)

# Define include paths
INCDIR = -I. -I$(ESPAINC) -I$(XML2INC)
HDF_INCDIR = -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
NCFLAGS  = $(EXTRA) $(INCDIR) $(HDF_INCDIR)

# Define the object libraries and paths
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -l_espa_format_conversion \
        -L$(HDFEOS_LIB) -lhdfeos \
        -L$(HDFLIB) -lmfhdf -ldf \
        -L$(HDFEOS_GCTPLIB) -lGctp \
        -L$(XML2LIB) -lxml2 \
        -L$(JPEGLIB) -ljpeg \
        -L$(LZMALIB) -llzma \
        -L$(SZIPLIB) -lsz \
        -L$(ZLIBLIB) -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define C executables
EXE = l8_sr

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(l8_bin_install_path)
	install -m 755 $(EXE) $(l8_bin_install_path)
	ln -sf $(l8_link_source_path)/$(EXE) $(link_path)/$(EXE)

#-----------------------------------------------------------------------------
clean:
	$(RM) -f *.o $(EXE)

#-----------------------------------------------------------------------------
$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<
