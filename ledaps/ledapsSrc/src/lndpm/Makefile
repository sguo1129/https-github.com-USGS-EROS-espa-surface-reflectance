#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
EXTRA  = -D_BSD_SOURCE -Wall -O2

TARGET = lndpm

INC = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC)

LIB = -L$(GEOTIFF_LIB) -lgeotiff -L$(TIFFLIB) -ltiff -L$(JPEGLIB) -ljpeg -lz -lm

OBJ =  lndpm.c

all: $(TARGET)

# Make the program
$(TARGET) : $(OBJ)
	$(CC) $(EXTRA) $(INC) $(OBJ) $(LIB)  -o $@

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) $(BIN)

#******************* End of make file *******************************


