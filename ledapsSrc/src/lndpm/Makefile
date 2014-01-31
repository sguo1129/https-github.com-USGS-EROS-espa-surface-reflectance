#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
EXTRA  = -g -D_BSD_SOURCE -Wall -O2

TARGET = lndpm

INC = -I. -I$(ESPAINC) -I$(XML2INC)

LIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -L$(XML2LIB) -lxml2 \
      -lz -lm

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


