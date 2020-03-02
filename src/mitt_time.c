/* Function to provide Unix time to Fortran
   The function input is an array with 6 elements
   following [year,month,day,hour,min,sec], and 
   the output is the number of days since 1970-01-01 00:00:00
   according to Unix epoch.

   (c) 2018 Pablo Saavedra G. (pablosaa@uni-bonn.de)
   Geophysical Institute, University of Bergen
   SEE LICENCE.TXT
 */

#include<stdio.h>
#include<time.h>
#include<stdlib.h>

double* F2UnixTime(int **dato, int ntime) {
  struct tm tv;
  time_t tt;
  double *NumberOfDays;

  NumberOfDays = (double *) malloc(ntime);
  for(int i=0; i<ntime; i++){
    tt = time(NULL);
    /* Filling tv with dummy local time */
    tv = *localtime(&tt);
  
    /* Assigning input Date to the struct time  */
    tv.tm_year = dato[i][0]-1900;  /* struct tm year starts on 1900 */
    tv.tm_mon  = dato[i][1];
    tv.tm_mday = dato[i][2];
    tv.tm_hour = dato[i][3];
    tv.tm_min  = dato[i][4];
    tv.tm_sec  = dato[i][5];

    /* setting Date as time_t (No time zone considered) */
    tt = mktime(&tv) - timezone;

    /* Converting seconds since Epoch to days*/
    NumberOfDays[i] = ((double) tt)/60/60/24;
  }
  
  return(NumberOfDays);
}

/* End of function */
