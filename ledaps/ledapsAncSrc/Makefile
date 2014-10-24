NCDFINC = /usr/include
NCDFLIB = /usr/lib
HDF5INC = /usr/include
HDF5LIB = /usr/lib



CC = gcc
EXTRA   = -ansi -Wall

INCDIR  = -I. -I$(NCDFINC) -I$(HDF5INC) -I$(HDFINC) -I$(JPEGINC) 
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

NCDF_LIBS = \
-L$(NCDFLIB) -lnetcdf -L$(HDF5LIB) -lhdf5_hl -lhdf5

LIB1 = -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg -lz -lm
LIB2 = $(NCDF_LIBS) -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg -lz -lm

TARGET1 = convert_ozone
OBJ1    = convert_ozone.o

TARGET2 = ncep_repackage
OBJ2    = ncep_reanalysis_surface_repackage.o

all: $(TARGET1) $(TARGET2)

$(OBJ1): $(INC1)

$(TARGET1): $(OBJ1)
	$(CC) $(EXTRA) $(OBJ1) -o $(TARGET1) $(LIB1)

$(OBJ2): $(INC2)

$(TARGET2): $(OBJ2)
	$(CC) $(EXTRA) $(OBJ2) -o $(TARGET2) $(LIB2)

clean:
	rm -f *.o $(TARGET1) $(TARGET2)

install:
	cp $(TARGET1) $(BIN)
	cp $(TARGET2) $(BIN)
#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -c $< -o $@
