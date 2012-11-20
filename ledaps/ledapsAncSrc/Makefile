EXTRA   = -Wall

INCDIR  = -I. -I$(HDFINC)
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

EXLIB	= -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1	= convert_ozone
OBJ1    = convert_ozone.o

TARGET2	= ncep_reanalysis_surface_repackage
OBJ2    = ncep_reanalysis_surface_repackage.o

all: $(TARGET1) $(TARGET2)

$(OBJ1): $(INC1)

$(TARGET1): $(OBJ1)
	$(CC) $(EXTRA) -o $(TARGET1) $(OBJ1) $(LOADLIB)

$(OBJ2): $(INC2)

$(TARGET2): $(OBJ2)
	$(CC) $(EXTRA) -o $(TARGET2) $(OBJ2) $(LOADLIB)

clean:
	rm -f *.o $(TARGET1) $(TARGET2)

install:
	cp $(TARGET1) $(BIN)
	cp $(TARGET2) $(BIN)
#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
