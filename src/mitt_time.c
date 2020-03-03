/* --------------------------------------------------------------
   Functions to provide Unix time to Fortran.
   There are 2 functions which can be called from Fortran:

   1) F2UnixTime: The function input is an array [6 x N], with 
   6 integers for the date following [year,month,day,hour,min,sec]
   and N the number inputs of dates.
   The output is a vector [N] with the number of days
   since 1970-01-01 00:00:00 according to Unix epoch.

   2) CH2UnixTime: The function input is a vector of chars, with
   every element of the vector should be a 19 length char with the
   form of "YYYY-MM-DD_hh:mi:se". 
   The output is a vector [N] with the number of days
   since 1970-01-01 00:00:00 according to Unix epoch.


   ---
   (c) 2018 Pablo Saavedra G. (pablo.saa@uib.no)
   Geophysical Institute, University of Bergen
   SEE LICENCE.TXT
   --------------------------------------------------------------*/

#include<stdio.h>
#include<time.h>


/*
  Main Function to convert an 6 element array to Number of Days.
  This function is called by the other variants depending on the 
  input variable they have: array of char.
*/
double Date2UnixTime(int datum[6]){

  struct tm tv;
  time_t tt;

  tt = time(NULL);

  /* Filling tv with dummy local time */
  tv = *localtime(&tt);

  /* Assigning input Date to the struct time  */
  tv.tm_year = datum[0] - 1900;  /* struct tm year starts on 1900 */
  tv.tm_mon  = datum[1] - 1;
  tv.tm_mday = datum[2];
  tv.tm_hour = datum[3];
  tv.tm_min  = datum[4];
  tv.tm_sec  = datum[5];

  /* setting Date as time_t (No time zone considered) */
  tt = mktime(&tv) - timezone;

  /* Converting seconds since Epoch to days*/
  double NumberOfDays = ((double) tt)/60/60/24;

  return(NumberOfDays);
}


/*
  Wrapper function for Date2UnixTime, when the input is an array type [6][N]:
*/
void F2UnixTime(int ntime, int dato[6][ntime], double NumberOfDays[ntime]) {

  for(int i=0; i<ntime; i++){

    int datum[6] = {dato[0][i], dato[1][i], dato[2][i], dato[3][i], dato[4][i], dato[5][i]};        
    NumberOfDays[i] = Date2UnixTime(datum);
    
  }
  
  return;
}


/*
  Wrapper function for Date2UnixTime, when the inut is a char array [N, "YYYY-MM-DD_hh:mi:se"]
*/
void CH2UnixTime(int ntime, int nlen, char dato[ntime], double NumberOfDays[ntime]) {
  
  for(int i=0, datum[6]; i<ntime; i++){
    
    sscanf(&dato[i*nlen], "%4d%*c%2d%*c%2d%*c%2d%*c%2d%*c%2d", &datum[0], &datum[1], &datum[2], &datum[3], &datum[4], &datum[5]);
    NumberOfDays[i] = Date2UnixTime(datum);

  }
  return;
}


void UnixTime2Date(double unixt, int datum[6]){

  time_t tt = unixt*60*60*24;  /* seconds since epoch */
  struct tm dd = *gmtime( &tt);
  datum[0] = dd.tm_year + 1900;
  datum[1] = dd.tm_mon + 1;
  datum[2] = dd.tm_mday;
  datum[3] = dd.tm_hour;
  datum[4] = dd.tm_min;
  datum[5] = dd.tm_sec;
  return;
}

/* End of functions
   ________________________________________________________ */


//#include<stdlib.h>
//#include<string.h>

/* char temp[nlen]; */
/* strncpy(temp,  &dato[i*nlen], nlen); */


