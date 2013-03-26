EXTRA   = -D_BSD_SOURCE -Wall -O2

INCDIR  = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

EXLIB	= -L$(GEOTIFF_LIB) -lgeotiff -L$(TIFFLIB) -ltiff \
          -L$(HDFEOS_LIB) -lhdfeos -L$(HDFEOS_GCTPLIB) -lGctp \
          -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1	= lndcal
OBJ1    = lndcal.o param.o degdms.o space.o lut.o output.o cal.o util.o\
          names.o date.o myhdf.o mystring.o error.o input.o
INC1    = lndcal.h keyvalue.h param.h input.h lut.h output.h cal.h \
          names.h date.h myhdf.h mystring.h bool.h const.h error.h \
          space.h myproj.h myproj_const.h util.h

all: $(TARGET1)

x: $(TARGET1)

$(OBJ1): $(INC1)

$(TARGET1): $(OBJ1)
	$(CC) $(EXTRA) -o $(TARGET1) $(OBJ1) $(LOADLIB)

clean:
	rm -f *.o $(TARGET1)

install:
	cp $(TARGET1) $(BIN)
#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
