EXTRA   = -g -D_BSD_SOURCE -Wall -O2

INCDIR  = -I. -I$(HDFEOS_GCTPINC) -I$(XML2INC) -I$(ESPAINC)
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

EXLIB	= -L$(HDFEOS_GCTPLIB) -lGctp \
          -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
          -L$(XML2LIB) -lxml2 -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1	= lndcal
OBJ1    = lndcal.o param.o lut.o output.o cal.o util.o \
          date.o mystring.o error.o input.o
INC1    = lndcal.h keyvalue.h param.h input.h lut.h output.h cal.h \
          date.h mystring.h bool.h const.h error.h \
          myproj.h myproj_const.h util.h

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
