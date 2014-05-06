#cc -O0 -c qa_index_routine.c

#gfortran -ffixed-line-length-132 -O -o $1 $1.f sub_deg2utm.f sub_utm2deg.f newLUTldcm_subr.f subaeroret.f qa_index_routine.o -I $HDFINC -I JPEGINC -L $HDFLIB -L JPEGLIB -lmfhdf -ldf -ljpeg -lz -ldl -lm -lsz

C_OBJECTS = qa_index_routine.o



LDCMSR: qa_index_routine.o
	touch LDCMSR-v1.0
	touch LDCMSR

#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@

