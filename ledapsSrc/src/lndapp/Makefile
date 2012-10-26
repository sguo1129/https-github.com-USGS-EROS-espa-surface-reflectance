#############################################################################
# !make
#
# makefile name: lndapp.mk (for linux)
#
##!END
#############################################################################
CC = gcc
CFLAGS = -Wall

TARGET = lndapp

INC = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)

LIB = -L$(HDFEOS_LIB) -lhdfeos -L$(HDFEOS_GCTPLIB) -lGctp \
      -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm

OBJ = lndapp.o

all: $(TARGET)

$(OBJ) : $(INC_FILES)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LIB)

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) $(BIN)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

#******************* End of make file *******************************

