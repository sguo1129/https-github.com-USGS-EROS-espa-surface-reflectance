#include "util.h"
#include "error.h"
/*
!C****************************************************************************

!File: util.c
  
!Description: Utility Functions (these read radiometric tables)

!Revision History:
 Revision 1.0 2004/06/16
 Jonathan Kutler
 Original Version.


!Team Unique Header:
  This software was written for the LEDAPS project and is based on the software 
  written by the MODIS Land Science Team Support Group for the Laboratory for 
  Terrestrial Physics (Code 922) at the   National Aeronautics and Space Administration,
  Goddard Space Flight Center

 ! References and Credits:

      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
      J. Kutler 
      LEDAPS,                           Raytheon ITSS
      jonathan_l_kutler@raytheon.com    7501 Forbes Blvd, ste 103
      phone: 301-352-2152               Seabrook MD, 20706
  
 ! functions included : 
 
   mygetline   - C line read

!END****************************************************************************
*/

int mygetline(char* l,int m,FILE* f)
{
  if ( fgets(l,m,f) == NULL )return -1;
  else { if ( l[strlen(l)-1]=='\n' )l[strlen(l)-1]='\0'; return strlen(l);  }
}


int skipline(FILE* f, int n)
{
  char l[1024];
  int i;
  for ( i=0; i<n; i++)if ( fgets(l,1024,f) == NULL )break;
  return i+1;
}


/* this program zooms an integer array */
void zoomIt(int16 *o,int16 *i,int n,int z)
  {int j,c=n*z;for(j=0;j<c;j++)o[j]=i[j/z];}


/* this program zooms a char array */
void zoomIt8(unsigned char*o,unsigned char*i,int n,int z)
  {int j,c=n*z;for(j=0;j<c;j++)o[j]=i[j/z];}


int big_endian()
  {long int i=1,j=0; char b=1; memcpy(&j,&b,1); return i!=j;}
