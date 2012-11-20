#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "hdf.h"
#include "mfhdf.h"
#define PRINT_FLAG (0)
#define MAX_STR_LEN (1024)
/* NLINE changed for OMI (1 degree res.) */
/*#define NLINE (12)*/
int NLINE, longline, shortline;

int verbose;
int GetLine(FILE *fp, char *s);
int read_ozone(char* fname, short int** data, int* doy, int* year, int* nlats,
               int* nlons, float* minlat, float* minlon, float* maxlat,
               float* maxlon, float* latsteps, float* lonsteps, float* lat_array, float* lon_array);

int main(int argc,char **argv) {
  int32 sdsout_id;
  int32 dimout_id;
	int32 sdin_id,sdout_id,sds_id,dim_id,hdf_ind;
	int32 nb_datasets,nb_globattrs;
	int32 status,index,count,dims[MAX_VAR_DIMS],start[MAX_VAR_DIMS];
        int32 edge[MAX_VAR_DIMS];
	int32 attr_index;
	int32 rank,data_type,attributes;
	int32 dim_index,dim_sizes[MAX_VAR_DIMS];
	int32 sds_r1,sds_r2,sds_r3,sds_t3,sds_t4,sds_qc;
	int32 sds_lat,sds_lon,sds_ts,sds_tv,sds_phi;
	float32 ul_lat,ul_lon,lr_lat,lr_lon;
	char name[MAX_NC_NAME];
	char names[2][MAX_NC_NAME];
	VOIDP buffer;
	char label[200], unit[200], format[200], coordsys[200];
	int write_metadata;
	float32 scalef=1.0;
	float32 addoff=0.0;
        char* oz_units= "Dobson";
	int16 doy;

   int32  num_lones,vgroup_id;  
   int32 *ref_array;    
   char   vgroup_name[VGNAMELENMAX], vgroup_class[VGNAMELENMAX];	
	char list_of_sds[4][10]={"pres","pr_wtr"};
   int idoy,year,nlats,nlons;
   float minlat, minlon, maxlat, maxlon, latsteps,lonsteps;
   short int* data;
   int i,ival;
   int16 base_date[3];
   float lat_array[1024], lon_array[1024];
   char* platform= argc>3 ? argv[3] : (char*) "Earthprobe";


	if (argc<3) {
		fprintf(stderr,"usage: %s <input> <output> <platform>\n",argv[0]);
		exit(-1);
	}
	verbose=0;
/****
	open input
****/
/* determine NLINE by sensor - Feng */
	if(strstr(argv[1], "omi")) {
	  NLINE = 15;
	  longline = 25;
	  shortline = 10;
	}
	else {
	  NLINE =12;
	  longline =25; 
	  shortline = 13;
	}
	printf("%d %d %d\n", NLINE, longline, shortline);
  read_ozone(argv[1], &data, &idoy, &year, &nlats,&nlons, &minlat, &minlon, &maxlat,
                      &maxlon, &latsteps, &lonsteps, lat_array, lon_array);

  if ( nlats!=180 || fabs(minlat+ 89.500)>0.0001 || fabs(maxlat- 89.500)>0.0001 || 
       nlons!=288 || fabs(minlon+179.375)>0.0001 || fabs(maxlon-179.375)>0.0001 ||
       fabs(latsteps-1.0)>0.0001 || fabs(lonsteps-1.25)>0.0001 )printf(
   "*** unexpected values ***\n nlats=%d nlons=%d minlat=%f maxlat=%f minlon=%f maxlon=%f\n"
   ,nlats,nlons,minlat,maxlat,minlon,maxlon);



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
	doy=(int16)idoy;
	
/****
	Determine the contents of the file
****/
 	if (write_metadata) {
/****
		Copy global attributes to output if creating a new file
		Write Day Of Year to output
****/
                base_date[0]= year; base_date[1]=1; base_date[2]=1;
		status = SDsetattr(sdout_id, "base_date", DFNT_INT16,3,base_date);
		status = SDsetattr(sdout_id, "Day Of Year", DFNT_INT16,1,&doy);
		
                count=strlen(platform)+1;
		printf("*** platform=(%s) len=%d ***\n",platform,count);
 		status=SDsetattr(sdout_id,"Platform",DFNT_CHAR8,count,(void*)platform);
/****
		Copy lat/lon SDS
****/

/************************/
/**** Latitude (lat) ****/
/************************/

	dim_sizes[0]= nlats;
        strcpy(name,"lat");
        strcpy(names[0],"lat");
        rank=1;
	if ((sdsout_id=SDcreate(sdout_id,name,DFNT_FLOAT32,rank,&dim_sizes[0]))<0) 
		return -1;

	start[0]=0;
        edge[0]=dim_sizes[0];
	if (SDwritedata(sdsout_id,start,NULL,edge,lat_array)<0)
		return -1;

        dimout_id=SDgetdimid(sdsout_id, 0);
	SDsetdimname(dimout_id, name); 

 
/*************************/
/**** Longitude (lon) ****/
/*************************/


	dim_sizes[0]= nlons;
        strcpy(name,"lon");
        strcpy(names[1],"lon");
        rank=1;
	if ((sdsout_id=SDcreate(sdout_id,name,DFNT_FLOAT32,rank,&dim_sizes[0]))<0) 
		return -1;


	start[0]=0;
        edge[0]=dim_sizes[0];
	if (SDwritedata(sdsout_id,start,NULL,edge,lon_array)<0)
		return -1;

        dimout_id=SDgetdimid(sdsout_id, 0);
	SDsetdimname(dimout_id, name); 

/************************/
/**** Ozone (ozone)  ****/
/************************/

	dim_sizes[0]= nlats;
	dim_sizes[1]= nlons;
        strcpy(name,"ozone");
        rank=2;
	if ((sdsout_id=SDcreate(sdout_id,name,DFNT_INT16,rank,&dim_sizes[0]))<0) 
		return -1;

	start[0]=0;
	start[1]=0;
        edge[0]=dim_sizes[0];
        edge[1]=dim_sizes[1];
	if (SDwritedata(sdsout_id,start,NULL,edge,data)<0)
		return -1;

	for (index=0;index<rank;index++) {
          dimout_id=SDgetdimid(sdsout_id, index);
          SDsetdimname(dimout_id, names[index]); 
	}

        strcpy(name,"scale_factor");
        count=1;
        data_type= DFNT_FLOAT32;
        SDsetattr(sdsout_id, name, data_type, count, (void*)&scalef);

        strcpy(name,"add_offset");
        count=1;
        data_type= DFNT_FLOAT32;
        SDsetattr(sdsout_id, name, data_type, count, (void*)&addoff);

        strcpy(name,"units");
        count=strlen(oz_units)+1;
        data_type= DFNT_CHAR8;
        SDsetattr(sdsout_id, name, data_type, count, (void*)oz_units);

	}


/****
		Close input & output
****/
	status = SDend(sdout_id);
}



int read_ozone(char* fname, short int** out_data, int* doy, int* year, int* nlats,
               int* nlons, float* minlat, float* minlon, float* maxlat,
               float* maxlon, float* latsteps, float* lonsteps,
               float* lat_array, float* lon_array)
{ 
 FILE *fp;
 char line[MAX_STR_LEN+1];
 char dtype[10],month[4],valbuf[4],cyear[5];
 char c_day_of_month[3];
 int day_of_month;
 int irow,icol,i,ival;
 int answer,ilat;
 int len,number;
 short int* data;
 float mylat;
 float comp_lat,comp_lon;
 char hemisphere1, hemisphere2;
 /*
 Day: 260 Sep 17, 2001    EP/TOMS    NRT OZONE    GEN:01.271 Asc LECT: 11:09 AM 
012345678 1 2345678 2 2345678 3 2345678 4 2345678 5 2345678 6 2345678 7 2345678
 Longitudes:  288 bins centered on 179.375 W to 179.375 E  (1.25 degree steps)  
 Latitudes :  180 bins centered on  89.5   S to  89.5   N  (1.00 degree steps)  
 */

 fp = fopen(fname, "r");
 len = GetLine(fp,line);
 sscanf(line," Day: %3d ",doy);
 strncpy(month,&line[10],3); month[3]='\0';
 strncpy(c_day_of_month,&line[14],2); c_day_of_month[2]='\0';
 strncpy(dtype,&line[26],7); dtype[7]='\0';
 day_of_month= atoi( c_day_of_month );
 *year= atoi( &line[18] );
 printf("year=%d month=%s day_of_month=%d\n",*year,month,day_of_month);

 len = GetLine(fp,line);
 sscanf(line," Longitudes:  %3d bins centered on %7f W to %7f E  (%5f ", 
        nlons, minlon, maxlon, lonsteps);
 /* OMI data used different index - Feng */
 if(strstr(fname, "omi")) {
   hemisphere1= line[42];  hemisphere2= line[55];  
 }
 else {
 hemisphere1= line[43];  hemisphere2= line[56];
 } 
 *minlon *= (hemisphere1=='S'||hemisphere1=='W'?-1.0:1.0);
 *maxlon *= (hemisphere2=='S'||hemisphere2=='W'?-1.0:1.0);

 len = GetLine(fp,line);
 sscanf(line," Latitudes :  %3d bins centered on %7f S to %7f N  (%4f ",
        nlats, minlat, maxlat, latsteps);
 /* OMI data used different index - Feng */
 if(strstr(fname, "omi")) {
   hemisphere1= line[42];  hemisphere2= line[55];  
 }
 else {
 hemisphere1= line[43];  hemisphere2= line[56]; 
 }
 *minlat *= (hemisphere1=='S'||hemisphere1=='W'?-1.0:1.0);
 *maxlat *= (hemisphere2=='S'||hemisphere2=='W'?-1.0:1.0);

 comp_lat= *minlat;
 comp_lon= *minlon;

 data = (short int *)calloc(((*nlats)*(*nlons)), sizeof(short int));

 for (icol=0; icol<*nlons; icol++)
   {
   lon_array[icol]= comp_lon;
   comp_lon += *lonsteps;
   } 

 ilat=0;
 for (irow=0; irow<*nlats; irow++)
 {
   ival=0;
   for (icol=0; icol<NLINE; icol++)
     {
     len= GetLine(fp,line);
     for ( i=0; i<(icol<(NLINE-1)?longline:shortline); i++ )
       {
       
       strncpy(valbuf,&line[i*3+1],3);
       valbuf[3]='\0';
       sscanf(valbuf,"%3d",&number); 
       data[(*nlats-ilat-1)*(*nlons)+ival++]= number;
       /*
       sscanf(&line[i*3+1],"%3d",&data[ival++]);
       */
       }
     }
   sscanf(&line[i*3+1],"   lat =   %f",&mylat);
   if (0 )   printf("*** in=(%s) mylat=(%f) complat=(%f) ***\n",&line[i*3+1],mylat,comp_lat);
   lat_array[irow]= comp_lat;
   comp_lat += *latsteps;
   ilat++;
 }

 *out_data= data;
 return 0;
}
int GetLine(FILE *fp, char *s)
{
  int i=0, c;
  while ((c = fgetc(fp)) != EOF)
  {
    s[i++]= c;
    if (c == '\n' || i >= MAX_STR_LEN ) break;
  }
  s[i-1]= '\0';
  return i;
}
