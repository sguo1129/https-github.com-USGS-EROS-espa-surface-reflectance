#ifndef UTIL_HPP
#define UTIL_HPP
/****************************************************************************/
#include <stdio.h>      /* for true, false                                  */
#include <string.h>     /* for strlen                                       */
#include <math.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <time.h>
#define max(A,B)(A>B?A:B)
#define min(A,B)(A>B?B:A)
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))
#define bounded(A,B,C)(A>B?(A<C?A:C):B)
#define FABS(A)(A>0?A:-1.0*A)
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "lndcal.h"
#include "date.h"
/*************************************************************************/
int mygetline(char* line, int maxlen, FILE* fin);
int skipline(FILE* f, int n);
void zoomIt(int16*o,int16*i,int n,int z);
void zoomIt8(unsigned char*o,unsigned char*i,int n,int z);
int big_endian();
/*************************************************************************/
#endif /* UTIL_HPP */
