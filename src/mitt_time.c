/* --------------------------------------------------------------
   Function to provide Unix time to Fortran
   The function input is an array with 6 elements
   following [year,month,day,hour,min,sec], and 
   the output is the number of days since 1970-01-01 00:00:00
   according to Unix epoch.

   ---
   (c) 2018 Pablo Saavedra G. (pablosaa@uni-bonn.de)
   Geophysical Institute, University of Bergen
   SEE LICENCE.TXT
   --------------------------------------------------------------*/

#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<string.h>

void F2UnixTime(int ntime, int dato[6][ntime], double NumberOfDays[ntime]) {

  for(int i=0; i<ntime; i++){
    struct tm tv;
    time_t tt;

    tt = time(NULL);

    /* Filling tv with dummy local time */
    tv = *localtime(&tt);

    /* Assigning input Date to the struct time  */
    tv.tm_year = dato[0][i]-1900;  /* struct tm year starts on 1900 */
    tv.tm_mon  = dato[1][i]-1;
    tv.tm_mday = dato[2][i];
    tv.tm_hour = dato[3][i];
    tv.tm_min  = dato[4][i];
    tv.tm_sec  = dato[5][i];

    /* setting Date as time_t (No time zone considered) */
    tt = mktime(&tv) - timezone;

    /* Converting seconds since Epoch to days*/
    NumberOfDays[i] = ((double) tt)/60/60/24;
  }
  
  return;
}


void CH2UnixTime(int ntime, int nlen, char dato[ntime], double NumberOfDays[ntime]) {
  
  for(int i=0; i<ntime; i++){
    struct tm tv;
    time_t tt;
    int yy,mm,dd,hh,mi,se;
    char temp[nlen];
    strncpy(temp,  &dato[i*nlen], nlen);
    sscanf(temp,"%4d%*c%2d%*c%2d%*c%2d%*c%2d%*c%2d",&yy,&mm,&dd,&hh,&mi,&se);

    tt = time(NULL);
    /* Filling tv with dummy local time */
    tv = *localtime(&tt);
  
    /* Assigning input Date to the struct time */
    tv.tm_year = yy - 1900;  /* struct tm year starts on 1900 */
    tv.tm_mon  = mm - 1;
    tv.tm_mday = dd;
    tv.tm_hour = hh;
    tv.tm_min  = mi;
    tv.tm_sec  = se;

    /* setting Date as time_t (No time zone considered) */
    tt = mktime(&tv) - timezone;

    /* Converting seconds since Epoch to days */
    NumberOfDays[i] = ((double) tt)/60/60/24;
  }
  return;
}
/* End of function */
