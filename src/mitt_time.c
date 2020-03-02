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
#include<string.h>

int F2UnixTime(int ntime, int dato[6][ntime], double NumberOfDays[ntime]) {
  //double *NumberOfDays;

  //NumberOfDays = (double *) malloc(sizeof(double)*ntime);

  for(int i=0; i<ntime; i++){
    struct tm tv;
    time_t tt;

    tt = time(NULL);

    /* Filling tv with dummy local time */
    tv = *localtime(&tt);

    /* Assigning input Date to the struct time  */
    tv.tm_year = dato[0][i]-1900;  /* struct tm year starts on 1900 */
    printf("Number of data: %d out of %d\n", i, ntime);  
    tv.tm_mon  = dato[1][i]-1;
    tv.tm_mday = dato[2][i];
    tv.tm_hour = dato[3][i];
    tv.tm_min  = dato[4][i];
    tv.tm_sec  = dato[5][i];

    printf("time is : %d %d %d %d\n", dato[0][i], dato[1][i], dato[2][i], dato[3][i]);
    /* setting Date as time_t (No time zone considered) */
    tt = mktime(&tv) - timezone;

    printf("seconds : %f\n", (float) tt);
    /* Converting seconds since Epoch to days*/
    NumberOfDays[i] = ((double) tt)/60/60/24;
    printf("days : %f\n", NumberOfDays[i]);
  }
  
  return(1);
}


int CH2UnixTime(int ntime, int nlen, char dato[ntime], double NumberOfDays[ntime]) {
  struct tm tv;
  time_t tt;
  char temp[nlen];
  
  for(int i=0; i<ntime; i++){
/*     tt = time(NULL); */
/*     /\* Filling tv with dummy local time *\/ */
/*     tv = *localtime(&tt); */
  
/*     /\* Assigning input Date to the struct time  *\/ */
/*     tv.tm_year = dato[i][0]-1900;  /\* struct tm year starts on 1900 *\/ */
/*     tv.tm_mon  = dato[i][1]; */
/*     tv.tm_mday = dato[i][2]; */
/*     tv.tm_hour = dato[i][3]; */
/*     tv.tm_min  = dato[i][4]; */
/*     tv.tm_sec  = dato[i][5]; */

/*     /\* setting Date as time_t (No time zone considered) *\/ */
/*     tt = mktime(&tv) - timezone; */

/*     /\* Converting seconds since Epoch to days*\/ */
/*     NumberOfDays[i] = ((double) tt)/60/60/24; */
    strncpy(temp,  &dato[i*nlen], nlen);
    printf("%d->  %s\n",i, temp);
  }

    
  return(1);
}
/* End of function */
