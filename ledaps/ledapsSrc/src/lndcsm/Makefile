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

# Define the include files
INC = bool.h const.h csm.h date.h error.h input.h keyvalue.h lndcsm.h lut.h \
      myhdf.h myproj_const.h myproj.h mystring.h names.h output.h param.h \
      space.h tiff.h util.h virbuf.h

# Define the source code and object files
SRC = \
      csm.c      \
      date.c     \
      degdms.c   \
      error.c    \
      input.c    \
      lndcsm.c   \
      lut.c      \
      myhdf.c    \
      mystring.c \
      names.c    \
      output.c   \
      param.c    \
      space.c    \
      tiff.c     \
      util.c     \
      virbuf.c
OBJ = $(SRC:.c=.o)

# Define include paths
INCDIR  = -I. -I$(ESPAINC) -I$(XML2INC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the object libraries and paths
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
        -L$(HDFLIB) -lmfhdf -ldf \
        -L$(HDFEOS_LIB) -lhdfeos \
        -L$(HDFEOS_GCTPLIB) -lGctp \
        -L$(JPEGLIB) -ljpeg \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define C executables
EXE = lndcsm

#-----------------------------------------------------------------------------
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(ledaps_bin_install_path)
	install -m 755 $(EXE) $(ledaps_bin_install_path)
	ln -sf $(ledaps_link_source_path)/$(EXE) $(link_path)/$(EXE)

#-----------------------------------------------------------------------------
clean:
	$(RM) -f *.o $(TARGET1)

#-----------------------------------------------------------------------------
$(OBJ): $(INC)

.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
