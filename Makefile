##############################################################################
# Compile commands from original implementation
##############################################################################
# cc -O0 -c qa_index_routine.c
#
# gfortran -ffixed-line-length-132 -O -o $1 $1.f sub_deg2utm.f sub_utm2deg.f newLUTldcm_subr.f subaeroret.f qa_index_routine.o -I $HDFINC -I JPEGINC -L $HDFLIB -L JPEGLIB -lmfhdf -ldf -ljpeg -lz -ldl -lm -lsz
##############################################################################

CC = cc
FC = gfortran

INCDIR  = -I. -I$(JPEGINC) -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
LIBDIR  = -L. -L$(JPEGLIB) -L$(TIFFLIB) -L$(GEOTIFF_LIB) -L$(HDFLIB) -L$(HDFEOS_LIB) -L$(HDFEOS_GCTPLIB)

C_EXTRA   = -D_BSD_SOURCE -Wall -O2
NCFLAGS   = $(CFLAGS) $(C_EXTRA) $(INCDIR)

F_EXTRA   = -Wall -O2 -ffixed-line-length-132

C_OBJECTS = qa_index_routine.o

F_OBJECTS = sub_deg2utm.o sub_utm2deg.o subaeroret.o newLUTldcm_subr.o LDCMSR-v1.0.o

TARGET = LDCMSR-v1.0

all: $(TARGET)

$(TARGET): $(C_OBJECTS) $(F_OBJECTS)
	$(FC) $(F_EXTRA) -o $(TARGET) $(F_OBJECTS) $(C_OBJECTS) $(INCDIR) $(LIBDIR) -lmfhdf -ldf -ljpeg -lz -ldl -lm

clean:
	rm $(TARGET) *.o

#
# Rules
#
.c.o:
	$(CC) $(NCFLAGS) -I$(HDFINC) -c $< -o $@

.f.o:
	$(FC) $(F_EXTRA) -c $< -o $@

