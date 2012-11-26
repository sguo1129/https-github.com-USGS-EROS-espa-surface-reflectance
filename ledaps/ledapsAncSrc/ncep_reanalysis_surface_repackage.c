#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "hdf4_netcdf.h"
#include "hdf.h"
#include "mfhdf.h"

int verbose;
int copy_sds(int32 sdin_id,char *sds_name,int32 first_time_index,int32 sdout_id);

int main(int argc,char **argv) {
	int32 sdin_id,sdout_id,sds_id;
	int32 nb_datasets,nb_globattrs;
	int32 index,count,start[MAX_VAR_DIMS],edge[MAX_VAR_DIMS],stride[MAX_VAR_DIMS];
	int32 rank,data_type,attributes;
	int32 dim_sizes[MAX_VAR_DIMS];
	int32 first_time_index;
	char name[MAX_NC_NAME];
	VOIDP buffer;
	double *timebuff;
	char title[200];
	int i,write_metadata;
	int16 doy,base_date[3];

	if (argc!=4) {
		fprintf(stderr,"usage: %s <input> <output> <doy>\n",argv[0]);
		exit(-1);
	}
	verbose=1;
/****
	open input
****/
	if ((sdin_id=SDstart(argv[1], DFACC_RDONLY))<0) {
   	fprintf(stderr,"can't open file %s\n",argv[1]);
		exit(-1);
	}
/****
	check and open output
****/
	write_metadata=0;
	if ((sdout_id=SDstart(argv[2], DFACC_RDONLY))<0) {
		if ((sdout_id=SDstart(argv[2], DFACC_CREATE))<0) {
   		fprintf(stderr,"can't create output %s\n",argv[2]);
			exit(-1);
		}
		write_metadata=1;
	} else {
		SDend(sdout_id);
		if ((sdout_id=SDstart(argv[2], DFACC_WRITE))<0) {
   		fprintf(stderr,"can't open output %s\n",argv[2]);
			exit(-1);
		}
	}
	doy=(int16)atoi(argv[3]);		
	
/****
	Determine the contents of the file
****/
   SDfileinfo(sdin_id, &nb_datasets, &nb_globattrs);
	if (verbose)
		printf ("Nb of global attributes = %d\n",nb_globattrs);
 	if (write_metadata) {
/****
		Copy global attributes to output if creating a new file
****/
  		for (index = 0; index < nb_globattrs; index++) {
     			SDattrinfo(sdin_id, index, name, &data_type, &count);
			if (verbose)
				printf ("attribute #%2d = %s (%d) size = %d\n",index,name,data_type,count);
			/* Allocate a buffer to hold the attribute data. */
			buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
			/* Read the attribute data. */
			SDreadattr(sdin_id, index, buffer);
			/* Write the attribute data. */
			SDsetattr(sdout_id, name, data_type,count,buffer);
	/* Save a copy of the title to get the year if needed (when base_date is absent) */
			if (!strcmp(name,"title")) {
				strncpy(title,buffer,count);
				title[count]='\0';
			}
			/* Free buffer */
			free(buffer);
  		}
/****
		if base_date is not available, create it
****/

/* check the content of the file */
   	for (index = 0; index < nb_datasets; index++) {
		sds_id = SDselect(sdin_id, index);
		SDgetinfo(sds_id, name, &rank, dim_sizes, \
                         &data_type, &attributes);
      				printf("SDS #%2d = %s  Rank = %d Type=%d  Nb Attrs=%d\n", index,name,rank,data_type,attributes);
				printf("\tDims = ");
				for (i=0;i<rank;i++) {
					if (i==0)
						printf("%d",dim_sizes[i]);
					else
						printf("x%d",dim_sizes[i]);
				}

                }
		printf("\n");
		index=2;
		sds_id = SDselect(sdin_id, index);
		SDgetinfo(sds_id, name, &rank, dim_sizes, \
                         &data_type, &attributes);
 		timebuff = (double*)malloc(dim_sizes[0]* DFKNTsize(data_type));
		start[0]=0;
		stride[0]=1;
		edge[0]=dim_sizes[0];
		SDreaddata(sds_id,&start[0],&stride[0],&edge[0],timebuff);
/*		printf("%f\n",timebuff[0]);
		printf("%f\n",timebuff[1]);*/
		base_date[0]=(int16)(timebuff[0]/(double)8765.81277)+1;
		printf("year %d\n",base_date[0]);
	        base_date[1]=1;
	        base_date[2]=1;
		SDsetattr(sdout_id, "base_date", DFNT_INT16,3,&base_date);

/****
		Write Day Of Year to output
****/
		SDsetattr(sdout_id, "Day Of Year", DFNT_INT16,1,&doy);
		
/****
		Copy lat/lon SDS
****/
		if (copy_sds(sdin_id,"lat",-1,sdout_id)) {
			fprintf(stderr,"ERROR: couldn't copy SDS lat ... ABORT\n");
			exit(-1);
		}
		if (copy_sds(sdin_id,"lon",-1,sdout_id)) {
			fprintf(stderr,"ERROR: couldn't copy SDS lon ... ABORT\n");
			exit(-1);
		}		
	}


	if (verbose)  	
		printf ("Nb of SDS = %d\n",nb_datasets);
/****
	Loop over SDS, get name and nb of attributes for each one
****/
   	for (index = 0; index < nb_datasets; index++) {
		sds_id = SDselect(sdin_id, index);
		SDgetinfo(sds_id, name, &rank, dim_sizes, \
                         &data_type, &attributes);

		if ((!strcmp(name,"pres"))||(!strcmp(name,"pr_wtr"))||(!strcmp(name,"slp"))||(!strcmp(name,"air"))) {
			if (verbose) {
      				printf("SDS #%2d = %s  Rank = %d Type=%d  Nb Attrs=%d\n", index,name,rank,data_type,attributes);
				printf("\tDims = ");
				for (i=0;i<rank;i++) {
					if (i==0)
						printf("%d",dim_sizes[i]);
					else
						printf("x%d",dim_sizes[i]);
				}
			}
			first_time_index=(doy-1)*4;
			if (copy_sds(sdin_id,name,first_time_index,sdout_id)) {
				fprintf(stderr,"ERROR: couldn't copy SDS %s ... ABORT\n",name);
				exit(-1);
			}		


		}
	}
/****
		Close input & output
****/
	SDend(sdin_id);
	SDend(sdout_id);
    return 0;
}   


int copy_sds(int32 sdin_id,char *sds_name,int32 first_time_index,int32 sdout_id) {
	int32 index,sds_id,sdsout_id,start[MAX_VAR_DIMS],edge[MAX_VAR_DIMS];
	char name[MAX_NC_NAME];
	int32 rank,data_type,attributes,count,dim_id,dimout_id;
	int32 dim_sizes[MAX_VAR_DIMS];
	char *buffer,*oneline;
	int buf_size,i,il;

			
        if ((index=SDnametoindex(sdin_id,sds_name))<0)
		return -1;
	sds_id = SDselect(sdin_id, index);
	SDgetinfo(sds_id, name, &rank, dim_sizes,&data_type, &attributes);
	if ((first_time_index >= 0) && (rank != 3))
		return -1;

	if ( rank==3 ) dim_sizes[0]=4; 

	if ((sdsout_id=SDcreate(sdout_id,name,data_type,rank,dim_sizes))<0) 
		return -1;
	
	if (first_time_index < 0) {
		buf_size=1;
		for (i=0;i<rank;i++) {
			buf_size *= dim_sizes[i];
			start[i]=0;
			edge[i]=dim_sizes[i];
		}
		buf_size *= DFKNTsize(data_type);
		if ((buffer=(VOIDP)malloc(buf_size))==NULL)
			return -1;
		if (SDreaddata(sds_id,start,NULL,edge,buffer)<0)
			return -1;
		if (SDwritedata(sdsout_id,start,NULL,edge,buffer)<0)
			return -1;
		free(buffer);
	} else {
		start[1]=0;
		edge[1]=dim_sizes[1];
		start[2]=0;
		edge[2]=dim_sizes[2];

		buf_size=dim_sizes[1]*dim_sizes[2]*DFKNTsize(data_type);

		if ((buffer=(VOIDP)malloc(buf_size))==NULL)
			return -1;
		if ((oneline=(VOIDP)malloc(dim_sizes[2]*DFKNTsize(data_type)))==NULL)
			return -1;
		edge[0]=1;
		for (i=0;i<4;i++) {
			start[0]=first_time_index+i;
			if (SDreaddata(sds_id,start,NULL,edge,buffer)<0)
				return -1;
			for (il=0;il<dim_sizes[1];il++) {
				memcpy(oneline,&buffer[(il*dim_sizes[2]+dim_sizes[2]/2)*DFKNTsize(data_type)],(dim_sizes[2]/2)*DFKNTsize(data_type));
				memcpy(&oneline[(dim_sizes[2]/2)*DFKNTsize(data_type)],&buffer[il*dim_sizes[2]*DFKNTsize(data_type)],(dim_sizes[2]-dim_sizes[2]/2)*DFKNTsize(data_type));
				memcpy(&buffer[il*dim_sizes[2]*DFKNTsize(data_type)],oneline,dim_sizes[2]*DFKNTsize(data_type));

			}
			start[0]=i;
			if (SDwritedata(sdsout_id,start,NULL,edge,buffer)<0)
				return -1;
		}
		free(buffer);
	}

/****
	Loop over SDS attributes and copy
****/
      	for (index=0;index<attributes;index++) {
		/* Get information about the data set attribute. */
		SDattrinfo(sds_id, index, name, &data_type, &count);
		/* Allocate a buffer to hold the attribute data. */
		buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
		/* Read attribute  */
		SDreadattr(sds_id, index, buffer);
		/* Write attribute */
		SDsetattr(sdsout_id, name, data_type, count, buffer);
		/* Free buffer */
		free(buffer);
	}
		
/****
	Loop over SDS dimensions
****/
	for (index=0;index<rank;index++) {
		dim_id = SDgetdimid(sds_id, index);
		SDdiminfo(dim_id, name, &count, &data_type, &attributes);
		dimout_id=SDgetdimid(sdsout_id, index);
		SDsetdimname(dimout_id, name); 
		if (verbose)
			printf ("    Dim = %s  size = %d	data type = %d  nb attrs = %d\n",name,count,data_type,attributes);
		if ((data_type != 0)&&(count != 0)) {
			/* Allocate a buffer to hold scale data. */
			buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
			/* Read scale data. */
			SDgetdimscale(dim_id, buffer);
			/* Write scale data. */
			SDsetdimscale(dimout_id, count, data_type, buffer); 
			/* Free buffer */
			free(buffer); 
		}
	}
/****
	Close SDS
****/
	SDendaccess(sds_id);
	SDendaccess(sdsout_id);
	return 0;
}

