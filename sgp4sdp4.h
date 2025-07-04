/*
 *  Header file for kelso
 */

#ifndef KELSO_H
#define KELSO_H 1


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>

/** Type definitions **/

/* Two-line-element satellite orbital data */
typedef struct
{
  double
	epoch, xndt2o, xndd6o, bstar,
	xincl, xnodeo, eo, omegao, xmo, xno;
  int
	catnr,  /* Catalogue Number  */
	elset,  /* Element Set       */
	revnum; /* Revolution Number */
  char
	sat_name[128], /* Satellite name string    */
	idesg[9];     /* International Designator */
} tle_t;

/* Geodetic position structure */
typedef struct
{
  double
	lat, lon, alt, theta;
} geodetic_t;

/* General three-dimensional vector structure */
typedef struct
{
  double
	x, y, z, w;
} vector_t;

/* Common arguments between deep-space functions */
typedef struct
{
  /* Used by dpinit part of Deep() */
  double
	eosq,sinio,cosio,betao,aodp,theta2,sing,cosg,
	betao2,xmdot,omgdot,xnodot,xnodp;
  /* Used by dpsec and dpper parts of Deep() */
  double
	xll,omgadf,xnode,em,xinc,xn,t;
  /* Used by thetg and Deep() */
  double
	ds50;
} deep_arg_t;


/** Table of constant values **/
#define de2ra    1.74532925E-2   /* Degrees to Radians */
#define pi       3.1415926535898 /* Pi */
#define pio2     1.5707963267949 /* Pi/2 */
#define x3pio2   4.71238898      /* 3*Pi/2 */
#define twopi    6.2831853071796 /* 2*Pi  */
#define e6a      1.0E-6
#define tothrd   6.6666667E-1    /* 2/3 */
#define xj2      1.0826158E-3    /* J2 Harmonic */
#define xj3     -2.53881E-6      /* J3 Harmonic */
#define xj4     -1.65597E-6      /* J4 Harmonic */
#define xke      7.43669161E-2
#define xkmper   6.378135E3      /* Earth radius km */
#define xmnpda   1.44E3          /* Minutes per day */
#define ae       1.0
#define ck2      5.413079E-4
#define ck4      6.209887E-7
#define f        3.352779E-3
#define ge       3.986008E5
#define s        1.012229
#define qoms2t   1.880279E-09
#define secday   8.6400E4        /* Seconds per day */
#define omega_E  1.0027379
#define omega_ER 6.3003879
#define zns      1.19459E-5
#define c1ss     2.9864797E-6
#define zes      1.675E-2
#define znl      1.5835218E-4
#define c1l      4.7968065E-7
#define zel      5.490E-2
#define zcosis   9.1744867E-1
#define zsinis   3.9785416E-1
#define zsings  -9.8088458E-1
#define zcosgs   1.945905E-1
#define zcoshs   1
#define zsinhs   0
#define q22      1.7891679E-6
#define q31      2.1460748E-6
#define q33      2.2123015E-7
#define g22      5.7686396
#define g32      9.5240898E-1
#define g44      1.8014998
#define g52      1.0508330
#define g54      4.4108898
#define root22   1.7891679E-6
#define root32   3.7393792E-7
#define root44   7.3636953E-9
#define root52   1.1428639E-7
#define root54   2.1765803E-9
#define thdt     4.3752691E-3
#define rho      1.5696615E-1
#define mfactor  7.292115E-5
#define sr       6.96000E5      /*Solar radius - kilometers (IAU 76)*/
#define AU       1.49597870E8   /*Astronomical unit - kilometers (IAU 76)*/

/* Entry points of Deep() */
#define dpinit   1 /* Deep-space initialization code */
#define dpsec    2 /* Deep-space secular code        */
#define dpper    3 /* Deep-space periodic code       */

/* Carriage return and line feed */
#define CR  0x0A
#define LF  0x0D


/* Flow control flag definitions */
#define ALL_FLAGS              -1
#define  SGP_INITIALIZED_FLAG  0x000001
#define SGP4_INITIALIZED_FLAG  0x000002
#define SDP4_INITIALIZED_FLAG  0x000004
#define SGP8_INITIALIZED_FLAG  0x000008
#define SDP8_INITIALIZED_FLAG  0x000010
#define SIMPLE_FLAG            0x000020
#define DEEP_SPACE_EPHEM_FLAG  0x000040
#define LUNAR_TERMS_DONE_FLAG  0x000080
#define NEW_EPHEMERIS_FLAG     0x000100
#define DO_LOOP_FLAG           0x000200
#define RESONANCE_FLAG         0x000400
#define SYNCHRONOUS_FLAG       0x000800
#define EPOCH_RESTART_FLAG     0x001000
#define VISIBLE_FLAG           0x002000
#define SAT_ECLIPSED_FLAG      0x004000


/* Funtion prototypes produced by cproto */
/* main.c */
int main(void);
/* sgp4sdp4.c */
void SGP4(double tsince, tle_t *tle, vector_t *pos, vector_t *vel);
void SDP4(double tsince, tle_t *tle, vector_t *pos, vector_t *vel);
void Deep(int ientry, tle_t *tle, deep_arg_t *deep_arg);
int isFlagSet(int flag);
int isFlagClear(int flag);
void SetFlag(int flag);
void ClearFlag(int flag);
/* sgp_in.c */
int Checksum_Good(char *tle_set);
int Good_Elements(char *tle_set);
void Convert_Satellite_Data(char *tle_set, tle_t *tle);
int Input_Tle_Set(char *tle_file, tle_t *tle);
void select_ephemeris(tle_t *tle);
/* sgp_math.c */
int Sign(double arg);
double Sqr(double arg);
double Cube(double arg);
double Radians(double arg);
double Degrees(double arg);
double ArcSin(double arg);
double ArcCos(double arg);
void Magnitude(vector_t *v);
void Vec_Add(vector_t *v1, vector_t *v2, vector_t *v3);
void Vec_Sub(vector_t *v1, vector_t *v2, vector_t *v3);
void Scalar_Multiply(double k, vector_t *v1, vector_t *v2);
void Scale_Vector(double k, vector_t *v);
double Dot(vector_t *v1, vector_t *v2);
double Angle(vector_t *v1, vector_t *v2);
void Cross(vector_t *v1, vector_t *v2, vector_t *v3);
void Normalize(vector_t *v);
double AcTan(double sinx, double cosx);
double FMod2p(double x);
double Modulus(double arg1, double arg2);
double Frac(double arg);
int Round(double arg);
double Int(double arg);
void Convert_Sat_State(vector_t *pos, vector_t *vel);
/* sgp_obs.c */
void Calculate_User_PosVel(double time, geodetic_t *geodetic, vector_t *obs_pos, vector_t *obs_vel);
void Calculate_LatLonAlt(double time, vector_t *pos, geodetic_t *geodetic);
void Calculate_Obs(double time, vector_t *pos, vector_t *vel, geodetic_t *geodetic, vector_t *obs_set);
void Calculate_RADec(double time, vector_t *pos, vector_t *vel, geodetic_t *geodetic, vector_t *obs_set);
/* sgp_time.c */
double Julian_Date_of_Epoch(double epoch);
double Epoch_Time(double jd);
int DOY(int yr, int mo, int dy);
double Fraction_of_Day(int hr, int mi, int se, suseconds_t usec);
void Calendar_Date(double jd, struct tm *cdate);
void Time_of_Day(double jd, struct tm *cdate);
double Julian_Date(struct tm *cdate, struct timeval *tv);
void Date_Time(double julian_date, struct tm *cdate);
int Check_Date(struct tm *cdate);
struct tm Time_to_UTC(struct tm *cdate);
struct tm Time_from_UTC(struct tm *cdate);
double JD_to_UTC(double jt);
double JD_from_UTC(double jt);
double Delta_ET(double year);
double Julian_Date_of_Year(double year);
double ThetaG(double epoch, deep_arg_t *deep_arg);
double ThetaG_JD(double jd);
void UTC_Calendar_Now(struct tm *cdate, struct timeval *tv);
/* solar.c */
void Calculate_Solar_Position(double time, vector_t *solar_vector);
int Sat_Eclipsed(vector_t *pos, vector_t *sol, double *depth);

  int
Checksum_Good( char *tle_set )
{
  int i, check_digit, value, checksum = 0;

  for(i = 0; i < 68; i++)
  {
	if( (tle_set[i] >= '0') && (tle_set[i] <= '9') )
	  value = tle_set[i] - '0';
	else if( tle_set[i] == '-' )
	  value = 1;
	else
	  value = 0;

	checksum += value;
  } /* End for(i = 0; i < 68; i++) */

  checksum %= 10;
  check_digit = tle_set[68] - '0';

  return( checksum == check_digit );
} /* Function Checksums_Good */

/*------------------------------------------------------------------*/

/* Carries out various checks on a TLE set to verify its validity */
/* tle_set is a character string holding the two lines read    */
/* from a text file containing NASA format Keplerian elements. */
  int
Good_Elements( char *tle_set )
{
  /* Verify checksum of both lines of a TLE set */
  if( !Checksum_Good(&tle_set[0]) || !Checksum_Good(&tle_set[69]) )
	return (0);
  /* Check the line number of each line */
  if( (tle_set[0] != '1') || (tle_set[69] != '2') )
	return (0);
  /* Verify that Satellite Number is same in both lines */
  if( strncmp( &tle_set[2], &tle_set[71], 5 ) != 0 )
	return (0);
  /* Check that various elements are in the right place */
  if( (tle_set[ 23] != '.') ||
	  (tle_set[ 34] != '.') ||
	  (tle_set[ 80] != '.') ||
	  (tle_set[ 89] != '.') ||
	  (tle_set[106] != '.') ||
	  (tle_set[115] != '.') ||
	  (tle_set[123] != '.') ||
	  (strncmp(&tle_set[61], " 0 ", 3) != 0)
	)
	return (0);

  return(1);
}  /* Function Good_Elements */

/*------------------------------------------------------------------*/

/* Converts the strings in a raw two-line element set  */
/* to their intended numerical values. No processing   */
/* of these values is done, e.g. from deg to rads etc. */
/* This is done in the select_ephemeris() function.    */
  void
Convert_Satellite_Data( char *tle_set, tle_t *tle )
{
  char buff[15];
  int i;

  /* Terminate white space on the right side of the name */
  for (i = strlen(tle->sat_name)-1; i >= 0 && isspace((int)tle->sat_name[i]); i--)
	  tle->sat_name[i] = 0;

  /** Decode Card 1 **/
  /* Satellite's catalogue number */
  strncpy( buff, &tle_set[2],5 );
  buff[5] = '\0';
  tle->catnr = atoi(buff);

  /* International Designator for satellite */
  strncpy( tle->idesg, &tle_set[9],8 );
  tle->idesg[8] = '\0';
  for (i = 7; i && isspace((int)tle->idesg[i]); i--)
	  tle->idesg[i] = '\0';

  /* Satellite's epoch */
  strncpy( buff, &tle_set[18],14 );
  buff[14] = '\0';
  tle->epoch = atof(buff);

  /* Satellite's First Time Derivative */
  strncpy( buff, &tle_set[33],10 );
  buff[10]='\0';
  tle->xndt2o = atof(buff);

  /* Satellite's Second Time Derivative */
  strncpy( buff, &tle_set[44],1 );
  buff[1] = '.';
  strncpy( &buff[2], &tle_set[45],5 );
  buff[7] = 'E';
  strncpy( &buff[8], &tle_set[50],2 );
  buff[10]='\0';
  tle->xndd6o = atof(buff);

  /* Satellite's bstar drag term */
  strncpy( buff, &tle_set[53],1 );
  buff[1] = '.';
  strncpy( &buff[2], &tle_set[54],5 );
  buff[7] = 'E';
  strncpy( &buff[8], &tle_set[59],2 );
  buff[10]='\0';
  tle->bstar = atof(buff);

  /* Element Number */
  strncpy( buff, &tle_set[64],4 );
  buff[4]='\0';
  tle->elset = atoi(buff);

  /** Decode Card 2 **/
  /* Satellite's Orbital Inclination (degrees) */
  strncpy( buff, &tle_set[77], 8 );
  buff[8]='\0';
  tle->xincl = atof(buff);

  /* Satellite's RAAN (degrees) */
  strncpy( buff, &tle_set[86], 8 );
  buff[8]='\0';
  tle->xnodeo = atof(buff);

  /* Satellite's Orbital Eccentricity */
  buff[0] = '.';
  strncpy( &buff[1], &tle_set[95], 7 );
  buff[8]='\0';
  tle->eo = atof(buff);

  /* Satellite's Argument of Perigee (degrees) */
  strncpy( buff, &tle_set[103], 8 );
  buff[8]='\0';
  tle->omegao = atof(buff);

  /* Satellite's Mean Anomaly of Orbit (degrees) */
  strncpy( buff, &tle_set[112], 8 );
  buff[8]='\0';
  tle->xmo = atof(buff);

  /* Satellite's Mean Motion (rev/day) */
  strncpy( buff, &tle_set[121], 10 );
  buff[10]='\0';
  tle->xno = atof(buff);

  /* Satellite's Revolution number at epoch */
  strncpy( buff, &tle_set[132], 5 );
  buff[5]='\0';
  double rvn = atof(buff);
  tle->revnum = (int)rvn;

} /* Procedure Convert_Satellite_Data */

/*------------------------------------------------------------------*/

/* Opens TLE source file and reads one (first!) TLE set  */
/* including the name line 0. If only the two Keplerian  */
/* element lines exist, this function will fail. Returns */
/* -1 if file open fails and -2 if TLE set is not valid  */
  int
Input_Tle_Set( char *tleStr, tle_t *tle)
{
  int chr;  /* Used for inputting characters */

  /*char sat_name[sizeof(tle->sat_name)]; /* Temp var for satellite name 
  char tle_set[139]; /* Two lines of a TLE set 

  // File pointer for opening TLE source file
  FILE *fp;

  // Open TLE file, abort on failure
  if( (fp = fopen( tle_file, "r")) == NULL )
	return(-1); */

  /* Read the satellite's name 
  fgets(sat_name, sizeof(sat_name), fp);
  strncpy(tle->sat_name, sat_name, sizeof(tle->sat_name));
  tle->sat_name[sizeof(tle->sat_name)-1] = 0;

  /* Read in first line of TLE set 
  fgets( tle_set, 70, fp );

  /* Dump CR/LF, put back last character 
  while( ((chr = fgetc(fp)) == CR) || (chr == LF) );
  ungetc(chr,fp);

  /* Read in second line of TLE set and terminate string 
  fgets( &tle_set[69], 70, fp );
  tle_set[138]='\0';

  fclose(fp); */

  /* Check TLE set and abort if not valid */
  if( !Good_Elements(tleStr) )
	return(-2);

  /* Convert the TLE set to orbital elements */
  Convert_Satellite_Data( tleStr, tle );

  return(0);
} /* End of Input_Tle_Set() */

/*------------------------------------------------------------------*/

/* Selects the apropriate ephemeris type to be used */
/* for predictions according to the data in the TLE */
/* It also processes values in the tle set so that  */
/* they are apropriate for the sgp4/sdp4 routines   */
  void
select_ephemeris(tle_t *tle)
{
  double ao,xnodp,dd1,dd2,delo,temp,a1,del1,r1;

  /* Preprocess tle set */
  tle->xnodeo *= de2ra;
  tle-> omegao *= de2ra;
  tle->xmo *= de2ra;
  tle->xincl *= de2ra;
  temp = twopi/xmnpda/xmnpda;
  tle->xno = tle->xno*temp*xmnpda;
  tle->xndt2o *= temp;
  tle->xndd6o = tle->xndd6o*temp/xmnpda;
  tle->bstar /= ae;

  /* Period > 225 minutes is deep space */
  dd1 = (xke/tle->xno);
  dd2 = tothrd;
  a1 = pow(dd1, dd2);
  r1 = cos(tle->xincl);
  dd1 = (1.0-tle->eo*tle->eo);
  temp = ck2*1.5f*(r1*r1*3.0-1.0)/pow(dd1, 1.5);
  del1 = temp/(a1*a1);
  ao = a1*(1.0-del1*(tothrd*.5+del1*
		(del1*1.654320987654321+1.0)));
  delo = temp/(ao*ao);
  xnodp = tle->xno/(delo+1.0);

  /* Select a deep-space/near-earth ephemeris */
  if (twopi/xnodp/xmnpda >= .15625)
	SetFlag(DEEP_SPACE_EPHEM_FLAG);
  else
	ClearFlag(DEEP_SPACE_EPHEM_FLAG);

  return;
} /* End of select_ephemeris() */

/*------------------------------------------------------------------*/

/* Returns sign of a double */
  int
Sign(double arg)
{
  if( arg > 0 )
	return( 1 );
  else if( arg < 0 )
	return( -1 );
  else
	return( 0 );
} /* Function Sign*/

/*------------------------------------------------------------------*/

/* Returns square of a double */
  double
Sqr(double arg)
{
  return( arg*arg );
} /* Function Sqr */

/*------------------------------------------------------------------*/

/* Returns cube of a double */
  double
Cube(double arg)
{
  return( arg*arg*arg );
} /*Function Cube*/

/*------------------------------------------------------------------*/

/* Returns angle in radians from arg id degrees */
  double
Radians(double arg)
{
  return( arg*de2ra );
} /*Function Radians*/

/*------------------------------------------------------------------*/

/* Returns angle in degrees from arg in rads */
  double
Degrees(double arg)
{
  return( arg/de2ra );
} /*Function Degrees*/

/*------------------------------------------------------------------*/

/* Returns the arcsine of the argument */
  double
ArcSin(double arg)
{
  if( fabs(arg) >= 1 )
	return( Sign(arg)*pio2 );
  else
	return( atan(arg/sqrt(1-arg*arg)) );
} /*Function ArcSin*/

/*------------------------------------------------------------------*/

/* Returns orccosine of rgument */
  double
ArcCos(double arg)
{
  return( pio2 - ArcSin(arg) );
} /*Function ArcCos*/

/*------------------------------------------------------------------*/

/* Calculates scalar magnitude of a vector_t argument */
  void
Magnitude(vector_t *v)
{
  v->w = sqrt(Sqr(v->x) + Sqr(v->y) + Sqr(v->z));
} /*Procedure Magnitude*/

/*------------------------------------------------------------------*/

/* Adds vectors v1 and v2 together to produce v3 */
  void
Vec_Add(vector_t *v1, vector_t *v2, vector_t *v3)
{
  v3->x = v1->x + v2->x;
  v3->y = v1->y + v2->y;
  v3->z = v1->z + v2->z;

  Magnitude(v3);
} /*Procedure Vec_Add*/

/*------------------------------------------------------------------*/

/* Subtracts vector v2 from v1 to produce v3 */
  void
Vec_Sub(vector_t *v1, vector_t *v2, vector_t *v3)
{
  v3->x = v1->x - v2->x;
  v3->y = v1->y - v2->y;
  v3->z = v1->z - v2->z;

  Magnitude(v3);
} /*Procedure Vec_Sub*/

/*------------------------------------------------------------------*/

/* Multiplies the vector v1 by the scalar k to produce the vector v2 */
  void
Scalar_Multiply(double k, vector_t *v1, vector_t *v2)
{
  v2->x = k * v1->x;
  v2->y = k * v1->y;
  v2->z = k * v1->z;
  v2->w = fabs(k) * v1->w;
} /*Procedure Scalar_Multiply*/

/*------------------------------------------------------------------*/

/* Multiplies the vector v1 by the scalar k */
  void
Scale_Vector(double k, vector_t *v)
{
  v->x *= k;
  v->y *= k;
  v->z *= k;
  Magnitude(v);
} /* Procedure Scale_Vector */

/*------------------------------------------------------------------*/

/* Returns the dot product of two vectors */
  double
Dot(vector_t *v1, vector_t *v2)
{
  return( v1->x*v2->x + v1->y*v2->y + v1->z*v2->z );
}  /*Function Dot*/

/*------------------------------------------------------------------*/

/* Calculates the angle between vectors v1 and v2 */
  double
Angle(vector_t *v1, vector_t *v2)
{
  Magnitude(v1);
  Magnitude(v2);
  return( ArcCos(Dot(v1,v2)/(v1->w*v2->w)) );
} /*Function Angle*/

/*------------------------------------------------------------------*/

/* Produces cross product of v1 and v2, and returns in v3 */
  void
Cross(vector_t *v1, vector_t *v2 ,vector_t *v3)
{
  v3->x = v1->y*v2->z - v1->z*v2->y;
  v3->y = v1->z*v2->x - v1->x*v2->z;
  v3->z = v1->x*v2->y - v1->y*v2->x;
  Magnitude(v3);
} /*Procedure Cross*/

/*------------------------------------------------------------------*/

/* Normalizes a vector */
  void
Normalize( vector_t *v )
{
  v->x /= v->w;
  v->y /= v->w;
  v->z /= v->w;
} /*Procedure Normalize*/

/*------------------------------------------------------------------*/

/* Four-quadrant arctan function */
  double
AcTan(double sinx, double cosx)
{
  if(cosx == 0)
  {
	if(sinx > 0)
	  return (pio2);
	else
	  return (x3pio2);
  }
  else
  {
	if(cosx > 0)
	{
	  if(sinx > 0)
		return ( atan(sinx/cosx) );
	  else
		return ( twopi + atan(sinx/cosx) );
	}
	else
	  return ( pi + atan(sinx/cosx) );
  }

} /* Function AcTan */

/*------------------------------------------------------------------*/

/* Returns mod 2pi of argument */
  double
FMod2p(double x)
{
  int i;
  double ret_val;

  ret_val = x;
  i = (int)(ret_val / twopi);
  ret_val -= (double)i * twopi;
  if (ret_val < 0) ret_val += twopi;

  return (ret_val);
} /* fmod2p */

/*------------------------------------------------------------------*/

/* Returns arg1 mod arg2 */
  double
Modulus(double arg1, double arg2)
{
  int i;
  double ret_val;

  ret_val = arg1;
  i = (int)(ret_val / arg2);
  ret_val -= (double)i * arg2;
  if (ret_val < 0) ret_val += arg2;

  return (ret_val);
} /* modulus */

/*------------------------------------------------------------------*/

/* Returns fractional part of double argument */
  double
Frac( double arg )
{
  return( arg - floor(arg) );
} /* Frac */

/*------------------------------------------------------------------*/

/* Returns argument rounded up to nearest integer */
  int
Round( double arg )
{
  double fl = floor(arg + 0.5);
  return( (int)fl );
} /* Round */

/*------------------------------------------------------------------*/

/* Returns the floor integer of a double arguement, as double */
  double
Int( double arg )
{
  return( floor(arg) );
} /* Int */

/*------------------------------------------------------------------*/

/* Converts the satellite's position and velocity  */
/* vectors from normalised values to km and km/sec */
  void
Convert_Sat_State( vector_t *pos, vector_t *vel )
{
  Scale_Vector( xkmper, pos );
  Scale_Vector( xkmper*xmnpda/secday, vel );

} /* Procedure Convert_Sat_State */

/*------------------------------------------------------------------*/

/* Procedure Calculate_User_PosVel passes the user's geodetic position */
/* and the time of interest and returns the ECI position and velocity  */
/* of the observer. The velocity calculation assumes the geodetic      */
/* position is stationary relative to the earth's surface.             */
void
Calculate_User_PosVel(
	double time,
	geodetic_t *geodetic,
	vector_t *obs_pos,
	vector_t *obs_vel)
{
  /* Reference:  The 1992 Astronomical Almanac, page K11. */

  double c,sq,achcp;

  geodetic->theta = FMod2p(ThetaG_JD(time) + geodetic->lon);/*LMST*/
  c = 1/sqrt(1 + f*(f - 2)*Sqr(sin(geodetic->lat)));
  sq = Sqr(1 - f)*c;
  achcp = (xkmper*c + geodetic->alt)*cos(geodetic->lat);
  obs_pos->x = achcp*cos(geodetic->theta);/*kilometers*/
  obs_pos->y = achcp*sin(geodetic->theta);
  obs_pos->z = (xkmper*sq + geodetic->alt)*sin(geodetic->lat);
  obs_vel->x = -mfactor*obs_pos->y;/*kilometers/second*/
  obs_vel->y =  mfactor*obs_pos->x;
  obs_vel->z =  0;
  Magnitude(obs_pos);
  Magnitude(obs_vel);
} /*Procedure Calculate_User_PosVel*/

/*------------------------------------------------------------------*/

/* Procedure Calculate_LatLonAlt will calculate the geodetic  */
/* position of an object given its ECI position pos and time. */
/* It is intended to be used to determine the ground track of */
/* a satellite.  The calculations  assume the earth to be an  */
/* oblate spheroid as defined in WGS '72.                     */
  void
Calculate_LatLonAlt(double time, vector_t *pos,  geodetic_t *geodetic)
{
  /* Reference:  The 1992 Astronomical Almanac, page K12. */

  double r,e2,phi,c;

  geodetic->theta = AcTan(pos->y,pos->x);/*radians*/
  geodetic->lon = FMod2p(geodetic->theta - ThetaG_JD(time));/*radians*/
  r = sqrt(Sqr(pos->x) + Sqr(pos->y));
  e2 = f*(2 - f);
  geodetic->lat = AcTan(pos->z,r);/*radians*/

  do
  {
	phi = geodetic->lat;
	c = 1/sqrt(1 - e2*Sqr(sin(phi)));
	geodetic->lat = AcTan(pos->z + xkmper*c*e2*sin(phi),r);
  }
  while(fabs(geodetic->lat - phi) >= 1E-10);

  geodetic->alt = r/cos(geodetic->lat) - xkmper*c;/*kilometers*/

  if( geodetic->lat > pio2 )
	geodetic->lat -= twopi;

} /*Procedure Calculate_LatLonAlt*/

/*------------------------------------------------------------------*/

/* The procedures Calculate_Obs and Calculate_RADec calculate         */
/* the *topocentric* coordinates of the object with ECI position,     */
/* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    */
/* The {obs_set} returned for Calculate_Obs consists of azimuth,      */
/* elevation, range, and range rate (in that order) with units of     */
/* radians, radians, kilometers, and kilometers/second, respectively. */
/* The WGS '72 geoid is used and the effect of atmospheric refraction */
/* (under standard temperature and pressure) is incorporated into the */
/* elevation calculation; the effect of atmospheric refraction on     */
/* range and range rate has not yet been quantified.                  */

/* The {obs_set} for Calculate_RADec consists of right ascension and  */
/* declination (in that order) in radians.  Again, calculations are   */
/* based on *topocentric* position using the WGS '72 geoid and        */
/* incorporating atmospheric refraction.                              */

void
Calculate_Obs(
	double time,
	vector_t *pos,
	vector_t *vel,
	geodetic_t *geodetic,
	vector_t *obs_set)
{
  double
	sin_lat,cos_lat,
	sin_theta,cos_theta,
    el,azim,
    top_s,top_e,top_z;

  vector_t
	obs_pos,obs_vel,range,rgvel;

  Calculate_User_PosVel(time, geodetic, &obs_pos, &obs_vel);

  range.x = pos->x - obs_pos.x;
  range.y = pos->y - obs_pos.y;
  range.z = pos->z - obs_pos.z;

  rgvel.x = vel->x - obs_vel.x;
  rgvel.y = vel->y - obs_vel.y;
  rgvel.z = vel->z - obs_vel.z;

  Magnitude(&range);

  sin_lat = sin(geodetic->lat);
  cos_lat = cos(geodetic->lat);
  sin_theta = sin(geodetic->theta);
  cos_theta = cos(geodetic->theta);
  top_s = sin_lat*cos_theta*range.x
	+ sin_lat*sin_theta*range.y
	- cos_lat*range.z;
  top_e = -sin_theta*range.x
	+ cos_theta*range.y;
  top_z = cos_lat*cos_theta*range.x
	+ cos_lat*sin_theta*range.y
	+ sin_lat*range.z;
  azim = atan(-top_e/top_s); /*Azimuth*/
  if( top_s > 0 )
	azim = azim + pi;
  if( azim < 0 )
	azim = azim + twopi;
  el = ArcSin(top_z/range.w);
  obs_set->x = azim;      /* Azimuth (radians)  */
  obs_set->y = el;        /* Elevation (radians)*/
  obs_set->z = range.w; /* Range (kilometers) */

  /*Range Rate (kilometers/second)*/
  obs_set->w = Dot(&range, &rgvel)/range.w;

  /* Corrections for atmospheric refraction */
  /* Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104    */
  /* Correction is meaningless when apparent elevation is below horizon */
  obs_set->y = obs_set->y + Radians((1.02/tan(Radians(Degrees(el)+
			10.3/(Degrees(el)+5.11))))/60);
  if( obs_set->y >= 0 )
	SetFlag(VISIBLE_FLAG);
  else
  {
	obs_set->y = el;  /*Reset to true elevation*/
	ClearFlag(VISIBLE_FLAG);
  } /*else*/
} /*Procedure Calculate_Obs*/

/*------------------------------------------------------------------*/

void
Calculate_RADec(
	double time,
	vector_t *pos,
	vector_t *vel,
	geodetic_t *geodetic,
	vector_t *obs_set)
{
  /* Reference:  Methods of Orbit Determination by  */
  /*                Pedro Ramon Escobal, pp. 401-402 */

  Calculate_Obs(time,pos,vel,geodetic,obs_set);
  if( isFlagSet(VISIBLE_FLAG) )
  {
	double
	  phi,theta,sin_theta,cos_theta,sin_phi,cos_phi,
	  az,el,Lxh,Lyh,Lzh,Sx,Ex,Zx,Sy,Ey,Zy,Sz,Ez,Zz,
	  Lx,Ly,Lz,cos_delta,sin_alpha,cos_alpha;

	az = obs_set->x;
	el = obs_set->y;
	phi   = geodetic->lat;
	theta = FMod2p(ThetaG_JD(time) + geodetic->lon);
	sin_theta = sin(theta);
	cos_theta = cos(theta);
	sin_phi = sin(phi);
	cos_phi = cos(phi);
	Lxh = -cos(az)*cos(el);
	Lyh =  sin(az)*cos(el);
	Lzh =  sin(el);
	Sx = sin_phi*cos_theta;
	Ex = -sin_theta;
	Zx = cos_theta*cos_phi;
	Sy = sin_phi*sin_theta;
	Ey = cos_theta;
	Zy = sin_theta*cos_phi;
	Sz = -cos_phi;
	Ez = 0;
	Zz = sin_phi;
	Lx = Sx*Lxh + Ex*Lyh + Zx*Lzh;
	Ly = Sy*Lxh + Ey*Lyh + Zy*Lzh;
	Lz = Sz*Lxh + Ez*Lyh + Zz*Lzh;
	obs_set->y = ArcSin(Lz);  /*Declination (radians)*/
	cos_delta = sqrt(1 - Sqr(Lz));
	sin_alpha = Ly/cos_delta;
	cos_alpha = Lx/cos_delta;
	obs_set->x = AcTan(sin_alpha,cos_alpha); /*Right Ascension (radians)*/
	obs_set->x = FMod2p(obs_set->x);
  }  /*if*/
} /* Procedure Calculate_RADec */

/*------------------------------------------------------------------*/

/* The function Julian_Date_of_Epoch returns the Julian Date of     */
/* an epoch specified in the format used in the NORAD two-line      */
/* element sets. It has been modified to support dates beyond       */
/* the year 1999 assuming that two-digit years in the range 00-56   */
/* correspond to 2000-2056. Until the two-line element set format   */
/* is changed, it is only valid for dates through 2056 December 31. */

  double
Julian_Date_of_Epoch(double epoch)
{
  double year,day;

  /* Modification to support Y2K */
  /* Valid 1957 through 2056     */
  day = modf(epoch*1E-3, &year)*1E3;
  if( year < 57 )
	year = year + 2000;
  else
	year = year + 1900;
  /* End modification */

  return( Julian_Date_of_Year(year) + day );
} /*Function Julian_Date_of_Epoch*/

/*------------------------------------------------------------------*/

/* Converts a Julian epoch to NORAD TLE epoch format */
  double
Epoch_Time(double jd)
{
  double yr,time,epoch_time;
  struct tm edate;

  Calendar_Date(jd, &edate);
  yr = edate.tm_year - 100*(edate.tm_year/100) ;
  time = Frac(jd + 0.5);
  epoch_time =  yr*1000
	+ DOY(edate.tm_year, edate.tm_mon, edate.tm_mday)
	+ time;

  return( epoch_time );
} /*Function Epoch_Time*/

/*------------------------------------------------------------------*/

/* The function DOY calculates the day of the year for the specified */
/* date. The calculation uses the rules for the Gregorian calendar   */
/* and is valid from the inception of that calendar system.          */
  int
DOY(int yr, int mo, int dy)
{
  const int days[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int i,day;

  day = 0;
  for( i = 0; i < mo-1; i++ )
	day += days[i];
  day = day + dy;

  /* Leap year correction */
  if((yr%4 == 0) && ((yr%100 != 0) || (yr%400 == 0)) && (mo>2) )
	day++;

  return( day );
} /*Function DOY*/

/*------------------------------------------------------------------*/

/* Fraction_of_Day calculates the fraction of */
/* a day passed at the specified input time.  */
  double
Fraction_of_Day(int hr,int mi,int se, suseconds_t usec)
{
  return( (hr + (mi + (se+usec*1e-6)/60.0)/60.0)/24.0 );
} /*Function Fraction_of_Day*/

/*------------------------------------------------------------------*/

/* The function Calendar_Date converts a Julian Date to a struct tm.   */
/* Only the members tm_year, tm_mon and tm_mday are calculated and set */
  void
Calendar_Date(double jd, struct tm *cdate)
{
  /* Astronomical Formulae for Calculators, Jean Meeus, pages 26-27 */
  int Z,month;
  double A,B,C,D,E,F,day,year,factor;

  factor = 0.5/secday/1000;
  F = Frac(jd + 0.5);
  if (F + factor >= 1.0)
  {
	jd = jd + factor;
	F  = 0.0;
  } /*if*/
  Z = Round(jd);
  if( Z < 2299161 )
	A = Z;
  else
  {
	double alpha = Int((Z - 1867216.25)/36524.25);
	A = Z + 1 + alpha - Int(alpha/4);
  } /*else*/
  B = A + 1524;
  C = Int((B - 122.1)/365.25);
  D = Int(365.25 * C);
  E = Int((B - D)/30.6001);
  day = B - D - Int(30.6001 * E) + F;

  if( E < 13.5 )
	month = Round(E - 1);
  else
	month = Round(E - 13);
  if( month > 2.5 )
	year = C - 4716;
  else
	year = C - 4715;

  cdate->tm_year = (int)year;
  cdate->tm_mon = month;
  double fl = floor(day);
  cdate->tm_mday = (int)fl;

} /*Function Calendar_Date*/

/*------------------------------------------------------------------*/

/* Time_of_Day takes a Julian Date and calculates the clock time */
/* portion of that date. Only tm_hour, tm_min and tm_sec are set */
  void
Time_of_Day(double jd, struct tm *cdate)
{
  int hr,mn,sc;
  double time;

  time = Frac(jd - 0.5)*secday;
  time = Round(time);
  double h = floor(time/3600.0);
  hr = (int)h;
  time = time - 3600.0*hr;
  if( hr == 24 )
	hr = 0;
  double m = floor(time/60.0);
  mn = (int)m;
  sc = (int)(time - 60.0 * (double)mn);
  cdate->tm_hour = hr;
  cdate->tm_min = mn;
  cdate->tm_sec = sc;

} /*Function Time_of_Day*/

/*------------------------------------------------------------------*/

/* The function Julian_Date converts a standard calendar   */
/* date and time to a Julian Date. The procedure Date_Time */
/* performs the inverse of this function. */
  double
Julian_Date(struct tm *cdate, struct timeval *tv)
{
  double julian_date;

  suseconds_t usec = 0;
  if (tv != NULL)
	  usec = tv->tv_usec;

  julian_date = Julian_Date_of_Year(cdate->tm_year) +
	DOY(cdate->tm_year,cdate->tm_mon,cdate->tm_mday) +
	Fraction_of_Day(cdate->tm_hour,cdate->tm_min,cdate->tm_sec,usec)
	+ 5.787037e-06; /* Round up to nearest 1 sec */

  return( julian_date );
} /*Function Julian_Date */

/*------------------------------------------------------------------*/


/*  Date_Time()
 *
 *  The function Date_Time() converts a Julian Date to
 *  standard calendar date and time. The function
 *  Julian_Date() performs the inverse of this function.
 */

  void
Date_Time(double julian_date, struct tm *cdate)
{
  time_t jtime;

  jtime = (time_t)((julian_date - 2440587.5) * 86400.0);
  *cdate = *gmtime( &jtime );

} /* End of Date_Time() */


/*------------------------------------------------------------------*/

/* The procedure Check_Date can be used as a check to see if a calendar    */
/* date and time are valid. It works by first converting the calendar      */
/* date and time to a Julian Date (which allows for irregularities, such   */
/* as a time greater than 24 hours) and then converting back and comparing.*/
  int
Check_Date(struct tm *cdate)
{
  double jt;
  struct tm chkdate;

  jt = Julian_Date(cdate, NULL);
  Date_Time(jt, &chkdate);

  if( (cdate->tm_year == chkdate.tm_year) &&
	  (cdate->tm_mon  == chkdate.tm_mon ) &&
	  (cdate->tm_mday == chkdate.tm_mday) &&
	  (cdate->tm_hour == chkdate.tm_hour) &&
	  (cdate->tm_min  == chkdate.tm_min ) &&
	  (cdate->tm_sec  == chkdate.tm_sec ) )
	return ( 1 );
  else
	return( 0 );

} /*Procedure Check_Date*/

/*------------------------------------------------------------------*/

/* Procedures Time_to_UTC and Time_from_UTC are used to  */
/* convert 'struct tm' dates between UTC and local time. */
/* The procedures JD_to_UTC and JD_from_UTC are used to  */
/* do the same thing working directly with Julian dates. */

  struct tm
Time_to_UTC(struct tm *cdate)
{
  time_t tdate;

  tdate = mktime(cdate);
  return( *gmtime(&tdate) );
} /*Procedure Time_to_UTC*/

/*------------------------------------------------------------------*/

  struct tm
Time_from_UTC(struct tm *cdate)
{
  time_t tdate;

  tdate = mktime(cdate);
  return( *localtime(&tdate) );
} /*Procedure Time_from_UTC*/

/*------------------------------------------------------------------*/

  double
JD_to_UTC(double jt)
{
  extern long timezone;
  struct tm cdate;

  time_t t = 0;

  cdate = *localtime( &t );
  jt = jt - timezone/secday;
  if( cdate.tm_isdst )
	jt= jt - 1.0/24.0;

  return( jt );
} /*Procedure JD_to_UTC*/

/*------------------------------------------------------------------*/

  double
JD_from_UTC(double jt)
{
  extern long timezone;
  struct tm cdate;
  time_t t = 0;

  cdate = *localtime( &t );
  jt = jt + timezone/secday;
  if( cdate.tm_isdst )
	jt= jt + 1.0/24.0;

  return( jt );
} /*Procedure JD_from_UTC*/

/*------------------------------------------------------------------*/

/* The function Delta_ET has been added to allow calculations on   */
/* the position of the sun.  It provides the difference between UT */
/* (approximately the same as UTC) and ET (now referred to as TDT).*/
/* This function is based on a least squares fit of data from 1950 */
/* to 1991 and will need to be updated periodically. */

  double
Delta_ET(double year)
{
  /* Values determined using data from 1950-1991 in the 1990
	 Astronomical Almanac.  See DELTA_ET.WQ1 for details. */

  double delta_et;

  delta_et = 26.465 + 0.747622*(year - 1950) +
	1.886913*sin(twopi*(year - 1975)/33);

  return( delta_et );
} /*Function Delta_ET*/

/*------------------------------------------------------------------*/

/* The function Julian_Date_of_Year calculates the Julian Date  */
/* of Day 0.0 of {year}. This function is used to calculate the */
/* Julian Date of any date by using Julian_Date_of_Year, DOY,   */
/* and Fraction_of_Day. */

  double
Julian_Date_of_Year(double year)
{
  /* Astronomical Formulae for Calculators, Jean Meeus, */
  /* pages 23-25. Calculate Julian Date of 0.0 Jan year */

  long A,B,i;
  double jdoy;

  year = year - 1.0;
  i = (long)(year / 100.0);
  A = i;
  i = A / 4;
  B = 2 - A + i;
  i = (long)(365.25 * year);
  i += (long)(30.6001 * 14.0);
  jdoy = (double)(i + 1720994.5 + B);

  return (jdoy);
}  /*Function Julian_Date_of_Year*/

/*------------------------------------------------------------------*/

/* The function ThetaG calculates the Greenwich Mean Sidereal Time */
/* for an epoch specified in the format used in the NORAD two-line */
/* element sets. It has now been adapted for dates beyond the year */
/* 1999, as described above. The function ThetaG_JD provides the   */
/* same calculation except that it is based on an input in the     */
/* form of a Julian Date. */

  double
ThetaG(double epoch, deep_arg_t *deep_arg)
{
  /* Reference:  The 1992 Astronomical Almanac, page B6. */

  double year,day,UT,jd,ThetaG;

  /* Modification to support Y2K */
  /* Valid 1957 through 2056     */
  day = modf(epoch*1E-3,&year)*1E3;
  if(year < 57)
	year += 2000;
  else
	year += 1900;
  /* End modification */

  UT   = modf(day,&day);
  jd   = Julian_Date_of_Year(year)+day;
  deep_arg->ds50 = jd-2433281.5+UT;
  ThetaG = FMod2p(6.3003880987*deep_arg->ds50+1.72944494);

  return (ThetaG);
} /* Function ThetaG */

/*------------------------------------------------------------------*/

  double
ThetaG_JD(double jd)
{
  /* Reference:  The 1992 Astronomical Almanac, page B6. */

  double UT,TU,GMST, ThetaG_JD;

  UT   = Frac(jd + 0.5);
  jd   = jd - UT;
  TU   = (jd - 2451545.0)/36525;
  GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6));
  GMST = Modulus(GMST + secday*omega_E*UT,secday);
  ThetaG_JD = twopi * GMST/secday;

  return( ThetaG_JD );
} /*Function ThetaG_JD*/

/*------------------------------------------------------------------*/

/* Gets calendar time from time() and produces a UTC calendar date */
  void
UTC_Calendar_Now( struct tm *cdate, struct timeval *tv )
{
  time_t t;

  t = time(0);
  *cdate = *gmtime(&t);
  cdate->tm_year += 1900;
  cdate->tm_mon += 1;

  if (tv != NULL)
	  gettimeofday(tv, NULL);

} /* End UTC_Calendar_Now */
/*------------------------------------------------------------------*/

/* SGP4 */
/* This function is used to calculate the position and velocity */
/* of near-earth (period < 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s. */
  void
SGP4(double tsince, tle_t *tle, vector_t *pos, vector_t *vel)
{
  static double
	aodp,aycof,c1,c4,c5,cosio,d2,d3,d4,delmo,omgcof,
	eta,omgdot,sinio,xnodp,sinmo,t2cof,t3cof,t4cof,t5cof,
	x1mth2,x3thm1,x7thm1,xmcof,xmdot,xnodcf,xnodot,xlcof;

  double
	cosuk,sinuk,rfdotk,vx,vy,vz,ux,uy,uz,xmy,xmx,
	cosnok,sinnok,cosik,sinik,rdotk,xinck,xnodek,uk,
	rk,cos2u,sin2u,u,sinu,cosu,betal,rfdot,rdot,r,pl,
	elsq,esine,ecose,epw,cosepw,x1m5th,xhdot1,tfour,
	sinepw,capu,ayn,xlt,aynl,xll,axn,xn,beta,xl,e,a,
	tcube,delm,delomg,templ,tempe,tempa,xnode,tsq,xmp,
	omega,xnoddf,omgadf,xmdf,a1,a3ovk2,ao,betao,betao2,
	c1sq,c2,c3,coef,coef1,del1,delo,eeta,eosq,etasq,
	perige,pinvsq,psisq,qoms24,s4,temp,temp1,temp2,
	temp3,temp4,temp5,temp6,theta2,theta4,tsi;

  int i;

  /* Initialization */
  if (isFlagClear(SGP4_INITIALIZED_FLAG))
  {
	SetFlag(SGP4_INITIALIZED_FLAG);

	/* Recover original mean motion (xnodp) and   */
	/* semimajor axis (aodp) from input elements. */
	a1 = pow(xke/tle->xno,tothrd);
	cosio = cos(tle->xincl);
	theta2 = cosio*cosio;
	x3thm1 = 3*theta2-1.0;
	eosq = tle->eo*tle->eo;
	betao2 = 1-eosq;
	betao = sqrt(betao2);
	del1 = 1.5*ck2*x3thm1/(a1*a1*betao*betao2);
	ao = a1*(1-del1*(0.5*tothrd+del1*(1+134/81*del1)));
	delo = 1.5*ck2*x3thm1/(ao*ao*betao*betao2);
	xnodp = tle->xno/(1+delo);
	aodp = ao/(1-delo);

	/* For perigee less than 220 kilometers, the "simple" flag is set */
	/* and the equations are truncated to linear variation in sqrt a  */
	/* and quadratic variation in mean anomaly.  Also, the c3 term,   */
	/* the delta omega term, and the delta m term are dropped.        */
	if((aodp*(1-tle->eo)/ae) < (220/xkmper+ae))
	  SetFlag(SIMPLE_FLAG);
	else
	  ClearFlag(SIMPLE_FLAG);

	/* For perigee below 156 km, the       */
	/* values of s and qoms2t are altered. */
	s4 = s;
	qoms24 = qoms2t;
	perige = (aodp*(1-tle->eo)-ae)*xkmper;
	if(perige < 156)
	{
	  if(perige <= 98)
		s4 = 20;
	  else
		s4 = perige-78;
	  qoms24 = pow((120-s4)*ae/xkmper,4);
	  s4 = s4/xkmper+ae;
	}; /* End of if(perige <= 98) */

	pinvsq = 1/(aodp*aodp*betao2*betao2);
	tsi = 1/(aodp-s4);
	eta = aodp*tle->eo*tsi;
	etasq = eta*eta;
	eeta = tle->eo*eta;
	psisq = fabs(1-etasq);
	coef = qoms24*pow(tsi,4);
	coef1 = coef/pow(psisq,3.5);
	c2 = coef1*xnodp*(aodp*(1+1.5*etasq+eeta*(4+etasq))+
		0.75*ck2*tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
	c1 = tle->bstar*c2;
	sinio = sin(tle->xincl);
	a3ovk2 = -xj3/ck2*pow(ae,3);
	c3 = coef*tsi*a3ovk2*xnodp*ae*sinio/tle->eo;
	x1mth2 = 1-theta2;
	c4 = 2*xnodp*coef1*aodp*betao2*(eta*(2+0.5*etasq)+
		tle->eo*(0.5+2*etasq)-2*ck2*tsi/(aodp*psisq)*
		(-3*x3thm1*(1-2*eeta+etasq*(1.5-0.5*eeta))+0.75*
		 x1mth2*(2*etasq-eeta*(1+etasq))*cos(2*tle->omegao)));
	c5 = 2*coef1*aodp*betao2*(1+2.75*(etasq+eeta)+eeta*etasq);
	theta4 = theta2*theta2;
	temp1 = 3*ck2*pinvsq*xnodp;
	temp2 = temp1*ck2*pinvsq;
	temp3 = 1.25*ck4*pinvsq*pinvsq*xnodp;
	xmdot = xnodp+0.5*temp1*betao*x3thm1+
	  0.0625*temp2*betao*(13-78*theta2+137*theta4);
	x1m5th = 1-5*theta2;
	omgdot = -0.5*temp1*x1m5th+0.0625*temp2*(7-114*theta2+
		395*theta4)+temp3*(3-36*theta2+49*theta4);
	xhdot1 = -temp1*cosio;
	xnodot = xhdot1+(0.5*temp2*(4-19*theta2)+
		2*temp3*(3-7*theta2))*cosio;
	omgcof = tle->bstar*c3*cos(tle->omegao);
	xmcof = -tothrd*coef*tle->bstar*ae/eeta;
	xnodcf = 3.5*betao2*xhdot1*c1;
	t2cof = 1.5*c1;
	xlcof = 0.125*a3ovk2*sinio*(3+5*cosio)/(1+cosio);
	aycof = 0.25*a3ovk2*sinio;
	delmo = pow(1+eta*cos(tle->xmo),3);
	sinmo = sin(tle->xmo);
	x7thm1 = 7*theta2-1;
	if (isFlagClear(SIMPLE_FLAG))
	{
	  c1sq = c1*c1;
	  d2 = 4*aodp*tsi*c1sq;
	  temp = d2*tsi*c1/3;
	  d3 = (17*aodp+s4)*temp;
	  d4 = 0.5*temp*aodp*tsi*(221*aodp+31*s4)*c1;
	  t3cof = d2+2*c1sq;
	  t4cof = 0.25*(3*d3+c1*(12*d2+10*c1sq));
	  t5cof = 0.2*(3*d4+12*c1*d3+6*d2*d2+15*c1sq*(2*d2+c1sq));
	}; /* End of if (isFlagClear(SIMPLE_FLAG)) */
  }; /* End of SGP4() initialization */

  /* Update for secular gravity and atmospheric drag. */
  xmdf = tle->xmo+xmdot*tsince;
  omgadf = tle->omegao+omgdot*tsince;
  xnoddf = tle->xnodeo+xnodot*tsince;
  omega = omgadf;
  xmp = xmdf;
  tsq = tsince*tsince;
  xnode = xnoddf+xnodcf*tsq;
  tempa = 1-c1*tsince;
  tempe = tle->bstar*c4*tsince;
  templ = t2cof*tsq;
  if (isFlagClear(SIMPLE_FLAG))
  {
	delomg = omgcof*tsince;
	delm = xmcof*(pow(1+eta*cos(xmdf),3)-delmo);
	temp = delomg+delm;
	xmp = xmdf+temp;
	omega = omgadf-temp;
	tcube = tsq*tsince;
	tfour = tsince*tcube;
	tempa = tempa-d2*tsq-d3*tcube-d4*tfour;
	tempe = tempe+tle->bstar*c5*(sin(xmp)-sinmo);
	templ = templ+t3cof*tcube+tfour*(t4cof+tsince*t5cof);
  }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

  a = aodp*pow(tempa,2);
  e = tle->eo-tempe;
  xl = xmp+omega+xnode+xnodp*templ;
  beta = sqrt(1-e*e);
  xn = xke/pow(a,1.5);

  /* Long period periodics */
  axn = e*cos(omega);
  temp = 1/(a*beta*beta);
  xll = temp*xlcof*axn;
  aynl = temp*aycof;
  xlt = xl+xll;
  ayn = e*sin(omega)+aynl;

  /* Solve Kepler's' Equation */
  capu = FMod2p(xlt-xnode);
  temp2 = capu;

  i = 0;
  do
  {
	sinepw = sin(temp2);
	cosepw = cos(temp2);
	temp3 = axn*sinepw;
	temp4 = ayn*cosepw;
	temp5 = axn*cosepw;
	temp6 = ayn*sinepw;
	epw = (capu-temp4+temp3-temp2)/(1-temp5-temp6)+temp2;
	if(fabs(epw-temp2) <= e6a)
	  break;
	temp2 = epw;
  }
  while( i++ < 10 );

  /* Short period preliminary quantities */
  ecose = temp5+temp6;
  esine = temp3-temp4;
  elsq = axn*axn+ayn*ayn;
  temp = 1-elsq;
  pl = a*temp;
  r = a*(1-ecose);
  temp1 = 1/r;
  rdot = xke*sqrt(a)*esine*temp1;
  rfdot = xke*sqrt(pl)*temp1;
  temp2 = a*temp1;
  betal = sqrt(temp);
  temp3 = 1/(1+betal);
  cosu = temp2*(cosepw-axn+ayn*esine*temp3);
  sinu = temp2*(sinepw-ayn-axn*esine*temp3);
  u = AcTan(sinu, cosu);
  sin2u = 2*sinu*cosu;
  cos2u = 2*cosu*cosu-1;
  temp = 1/pl;
  temp1 = ck2*temp;
  temp2 = temp1*temp;

  /* Update for short periodics */
  rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
  uk = u-0.25*temp2*x7thm1*sin2u;
  xnodek = xnode+1.5*temp2*cosio*sin2u;
  xinck = tle->xincl+1.5*temp2*cosio*sinio*cos2u;
  rdotk = rdot-xn*temp1*x1mth2*sin2u;
  rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);

  /* Orientation vectors */
  sinuk = sin(uk);
  cosuk = cos(uk);
  sinik = sin(xinck);
  cosik = cos(xinck);
  sinnok = sin(xnodek);
  cosnok = cos(xnodek);
  xmx = -sinnok*cosik;
  xmy = cosnok*cosik;
  ux = xmx*sinuk+cosnok*cosuk;
  uy = xmy*sinuk+sinnok*cosuk;
  uz = sinik*sinuk;
  vx = xmx*cosuk-cosnok*sinuk;
  vy = xmy*cosuk-sinnok*sinuk;
  vz = sinik*cosuk;

  /* Position and velocity */
  pos->x = rk*ux;
  pos->y = rk*uy;
  pos->z = rk*uz;
  vel->x = rdotk*ux+rfdotk*vx;
  vel->y = rdotk*uy+rfdotk*vy;
  vel->z = rdotk*uz+rfdotk*vz;

} /*SGP4*/

/*------------------------------------------------------------------*/

/* SDP4 */
/* This function is used to calculate the position and velocity */
/* of deep-space (period > 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s. */
  void
SDP4(double tsince, tle_t *tle, vector_t *pos, vector_t *vel)
{
  int i;

  static double
	x3thm1,c1,x1mth2,c4,xnodcf,t2cof,xlcof,aycof,x7thm1;

  double
	a,axn,ayn,aynl,beta,betal,capu,cos2u,cosepw,cosik,
	cosnok,cosu,cosuk,ecose,elsq,epw,esine,pl,theta4,
	rdot,rdotk,rfdot,rfdotk,rk,sin2u,sinepw,sinik,
	sinnok,sinu,sinuk,tempe,templ,tsq,u,uk,ux,uy,uz,
	vx,vy,vz,xinck,xl,xlt,xmam,xmdf,xmx,xmy,xnoddf,
	xnodek,xll,a1,a3ovk2,ao,c2,coef,coef1,x1m5th,
	xhdot1,del1,r,delo,eeta,eta,etasq,perige,
	psisq,tsi,qoms24,s4,pinvsq,temp,tempa,temp1,
	temp2,temp3,temp4,temp5,temp6;

  static deep_arg_t deep_arg;

  /* Initialization */
  if (isFlagClear(SDP4_INITIALIZED_FLAG))
  {
	SetFlag(SDP4_INITIALIZED_FLAG);

	/* Recover original mean motion (xnodp) and   */
	/* semimajor axis (aodp) from input elements. */
	a1 = pow(xke/tle->xno,tothrd);
	deep_arg.cosio = cos(tle->xincl);
	deep_arg.theta2 = deep_arg.cosio*deep_arg.cosio;
	x3thm1 = 3*deep_arg.theta2-1;
	deep_arg.eosq = tle->eo*tle->eo;
	deep_arg.betao2 = 1-deep_arg.eosq;
	deep_arg.betao = sqrt(deep_arg.betao2);
	del1 = 1.5*ck2*x3thm1/(a1*a1*deep_arg.betao*deep_arg.betao2);
	ao = a1*(1-del1*(0.5*tothrd+del1*(1+134/81*del1)));
	delo = 1.5*ck2*x3thm1/(ao*ao*deep_arg.betao*deep_arg.betao2);
	deep_arg.xnodp = tle->xno/(1+delo);
	deep_arg.aodp = ao/(1-delo);

	/* For perigee below 156 km, the values */
	/* of s and qoms2t are altered.         */
	s4 = s;
	qoms24 = qoms2t;
	perige = (deep_arg.aodp*(1-tle->eo)-ae)*xkmper;
	if(perige < 156)
	{
	  if(perige <= 98)
		s4 = 20;
	  else
		s4 = perige-78;
	  qoms24 = pow((120-s4)*ae/xkmper,4);
	  s4 = s4/xkmper+ae;
	}
	pinvsq = 1/(deep_arg.aodp*deep_arg.aodp*
		deep_arg.betao2*deep_arg.betao2);
	deep_arg.sing = sin(tle->omegao);
	deep_arg.cosg = cos(tle->omegao);
	tsi = 1/(deep_arg.aodp-s4);
	eta = deep_arg.aodp*tle->eo*tsi;
	etasq = eta*eta;
	eeta = tle->eo*eta;
	psisq = fabs(1-etasq);
	coef = qoms24*pow(tsi,4);
	coef1 = coef/pow(psisq,3.5);
	c2 = coef1*deep_arg.xnodp*(deep_arg.aodp*(1+1.5*etasq+eeta*
		  (4+etasq))+0.75*ck2*tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
	c1 = tle->bstar*c2;
	deep_arg.sinio = sin(tle->xincl);
	a3ovk2 = -xj3/ck2*pow(ae,3);
	x1mth2 = 1-deep_arg.theta2;
	c4 = 2*deep_arg.xnodp*coef1*deep_arg.aodp*deep_arg.betao2*
	  (eta*(2+0.5*etasq)+tle->eo*(0.5+2*etasq)-2*ck2*tsi/
	  (deep_arg.aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*
	  (1.5-0.5*eeta))+0.75*x1mth2*(2*etasq-eeta*(1+etasq))*
	  cos(2*tle->omegao)));
	theta4 = deep_arg.theta2*deep_arg.theta2;
	temp1 = 3*ck2*pinvsq*deep_arg.xnodp;
	temp2 = temp1*ck2*pinvsq;
	temp3 = 1.25*ck4*pinvsq*pinvsq*deep_arg.xnodp;
	deep_arg.xmdot = deep_arg.xnodp+0.5*temp1*deep_arg.betao*
	  x3thm1+0.0625*temp2*deep_arg.betao*
	  (13-78*deep_arg.theta2+137*theta4);
	x1m5th = 1-5*deep_arg.theta2;
	deep_arg.omgdot = -0.5*temp1*x1m5th+0.0625*temp2*
	  (7-114*deep_arg.theta2+395*theta4)+
	  temp3*(3-36*deep_arg.theta2+49*theta4);
	xhdot1 = -temp1*deep_arg.cosio;
	deep_arg.xnodot = xhdot1+(0.5*temp2*(4-19*deep_arg.theta2)+
		2*temp3*(3-7*deep_arg.theta2))*deep_arg.cosio;
	xnodcf = 3.5*deep_arg.betao2*xhdot1*c1;
	t2cof = 1.5*c1;
	xlcof = 0.125*a3ovk2*deep_arg.sinio*(3+5*deep_arg.cosio)/
	  (1+deep_arg.cosio);
	aycof = 0.25*a3ovk2*deep_arg.sinio;
	x7thm1 = 7*deep_arg.theta2-1;

	/* initialize Deep() */
	Deep(dpinit, tle, &deep_arg);
  }; /*End of SDP4() initialization */

  /* Update for secular gravity and atmospheric drag */
  xmdf = tle->xmo+deep_arg.xmdot*tsince;
  deep_arg.omgadf = tle->omegao+deep_arg.omgdot*tsince;
  xnoddf = tle->xnodeo+deep_arg.xnodot*tsince;
  tsq = tsince*tsince;
  deep_arg.xnode = xnoddf+xnodcf*tsq;
  tempa = 1-c1*tsince;
  tempe = tle->bstar*c4*tsince;
  templ = t2cof*tsq;
  deep_arg.xn = deep_arg.xnodp;

  /* Update for deep-space secular effects */
  deep_arg.xll = xmdf;
  deep_arg.t = tsince;

  Deep(dpsec, tle, &deep_arg);

  xmdf = deep_arg.xll;
  a = pow(xke/deep_arg.xn,tothrd)*tempa*tempa;
  deep_arg.em = deep_arg.em-tempe;
  xmam = xmdf+deep_arg.xnodp*templ;

  /* Update for deep-space periodic effects */
  deep_arg.xll = xmam;

  Deep(dpper, tle, &deep_arg);

  xmam = deep_arg.xll;
  xl = xmam+deep_arg.omgadf+deep_arg.xnode;
  beta = sqrt(1-deep_arg.em*deep_arg.em);
  deep_arg.xn = xke/pow(a,1.5);

  /* Long period periodics */
  axn = deep_arg.em*cos(deep_arg.omgadf);
  temp = 1/(a*beta*beta);
  xll = temp*xlcof*axn;
  aynl = temp*aycof;
  xlt = xl+xll;
  ayn = deep_arg.em*sin(deep_arg.omgadf)+aynl;

  /* Solve Kepler's Equation */
  capu = FMod2p(xlt-deep_arg.xnode);
  temp2 = capu;

  i = 0;
  do
  {
	sinepw = sin(temp2);
	cosepw = cos(temp2);
	temp3 = axn*sinepw;
	temp4 = ayn*cosepw;
	temp5 = axn*cosepw;
	temp6 = ayn*sinepw;
	epw = (capu-temp4+temp3-temp2)/(1-temp5-temp6)+temp2;
	if(fabs(epw-temp2) <= e6a)
	  break;
	temp2 = epw;
  }
  while( i++ < 10 );

  /* Short period preliminary quantities */
  ecose = temp5+temp6;
  esine = temp3-temp4;
  elsq = axn*axn+ayn*ayn;
  temp = 1-elsq;
  pl = a*temp;
  r = a*(1-ecose);
  temp1 = 1/r;
  rdot = xke*sqrt(a)*esine*temp1;
  rfdot = xke*sqrt(pl)*temp1;
  temp2 = a*temp1;
  betal = sqrt(temp);
  temp3 = 1/(1+betal);
  cosu = temp2*(cosepw-axn+ayn*esine*temp3);
  sinu = temp2*(sinepw-ayn-axn*esine*temp3);
  u = AcTan(sinu,cosu);
  sin2u = 2*sinu*cosu;
  cos2u = 2*cosu*cosu-1;
  temp = 1/pl;
  temp1 = ck2*temp;
  temp2 = temp1*temp;

  /* Update for short periodics */
  rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
  uk = u-0.25*temp2*x7thm1*sin2u;
  xnodek = deep_arg.xnode+1.5*temp2*deep_arg.cosio*sin2u;
  xinck = deep_arg.xinc+1.5*temp2*deep_arg.cosio*deep_arg.sinio*cos2u;
  rdotk = rdot-deep_arg.xn*temp1*x1mth2*sin2u;
  rfdotk = rfdot+deep_arg.xn*temp1*(x1mth2*cos2u+1.5*x3thm1);

  /* Orientation vectors */
  sinuk = sin(uk);
  cosuk = cos(uk);
  sinik = sin(xinck);
  cosik = cos(xinck);
  sinnok = sin(xnodek);
  cosnok = cos(xnodek);
  xmx = -sinnok*cosik;
  xmy = cosnok*cosik;
  ux = xmx*sinuk+cosnok*cosuk;
  uy = xmy*sinuk+sinnok*cosuk;
  uz = sinik*sinuk;
  vx = xmx*cosuk-cosnok*sinuk;
  vy = xmy*cosuk-sinnok*sinuk;
  vz = sinik*cosuk;

  /* Position and velocity */
  pos->x = rk*ux;
  pos->y = rk*uy;
  pos->z = rk*uz;
  vel->x = rdotk*ux+rfdotk*vx;
  vel->y = rdotk*uy+rfdotk*vy;
  vel->z = rdotk*uz+rfdotk*vz;

} /* SDP4 */

/*------------------------------------------------------------------*/

/* DEEP */
/* This function is used by SDP4 to add lunar and solar */
/* perturbation effects to deep-space orbit objects.    */
  void
Deep(int ientry, tle_t *tle, deep_arg_t *deep_arg)
{
  static double
	thgr,xnq,xqncl,omegaq,zmol,zmos,savtsn,ee2,e3,xi2,
	xl2,xl3,xl4,xgh2,xgh3,xgh4,xh2,xh3,sse,ssi,ssg,xi3,
	se2,si2,sl2,sgh2,sh2,se3,si3,sl3,sgh3,sh3,sl4,sgh4,
	ssl,ssh,d3210,d3222,d4410,d4422,d5220,d5232,d5421,
	d5433,del1,del2,del3,fasx2,fasx4,fasx6,xlamo,xfact,
	xni,atime,stepp,stepn,step2,preep,pl,sghs,xli,
	d2201,d2211,sghl,sh1,pinc,pe,shs,zsingl,zcosgl,
	zsinhl,zcoshl,zsinil,zcosil;

  double
	a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,ainv2,alfdp,aqnv,
	sgh,sini2,sinis,sinok,sh,si,sil,day,betdp,dalf,
	bfact,c,cc,cosis,cosok,cosq,ctem,f322,zx,zy,
	dbet,dls,eoc,eq,f2,f220,f221,f3,f311,f321,xnoh,
	f330,f441,f442,f522,f523,f542,f543,g200,g201,
	g211,pgh,ph,s1,s2,s3,s4,s5,s6,s7,se,sel,ses,xls,
	g300,g310,g322,g410,g422,g520,g521,g532,g533,gam,
	sinq,sinzf,sis,sl,sll,sls,stem,temp,temp1,x1,x2,
	x2li,x2omi,x3,x4,x5,x6,x7,x8,xl,xldot,xmao,xnddt,
	xndot,xno2,xnodce,xnoi,xomi,xpidot,z1,z11,z12,z13,
	z2,z21,z22,z23,z3,z31,z32,z33,ze,zf,zm,zn,
	zsing,zsinh,zsini,zcosg,zcosh,zcosi,delt=0,ft=0;

  switch(ientry)
  {
	case dpinit : /* Entrance for deep space initialization */
	  thgr = ThetaG(tle->epoch, deep_arg);
	  eq = tle->eo;
	  xnq = deep_arg->xnodp;
	  aqnv = 1/deep_arg->aodp;
	  xqncl = tle->xincl;
	  xmao = tle->xmo;
	  xpidot = deep_arg->omgdot+deep_arg->xnodot;
	  sinq = sin(tle->xnodeo);
	  cosq = cos(tle->xnodeo);
	  omegaq = tle->omegao;

	  /* Initialize lunar solar terms */
	  day = deep_arg->ds50+18261.5;  /*Days since 1900 Jan 0.5*/
	  if (day != preep)
	  {
		preep = day;
		xnodce = 4.5236020-9.2422029E-4*day;
		stem = sin(xnodce);
		ctem = cos(xnodce);
		zcosil = 0.91375164-0.03568096*ctem;
		zsinil = sqrt(1-zcosil*zcosil);
		zsinhl = 0.089683511*stem/zsinil;
		zcoshl = sqrt(1-zsinhl*zsinhl);
		c = 4.7199672+0.22997150*day;
		gam = 5.8351514+0.0019443680*day;
		zmol = FMod2p(c-gam);
		zx = 0.39785416*stem/zsinil;
		zy = zcoshl*ctem+0.91744867*zsinhl*stem;
		zx = AcTan(zx,zy);
		zx = gam+zx-xnodce;
		zcosgl = cos(zx);
		zsingl = sin(zx);
		zmos = 6.2565837+0.017201977*day;
		zmos = FMod2p(zmos);
	  } /* End if(day != preep) */

	  /* Do solar terms */
	  savtsn = 1E20;
	  zcosg = zcosgs;
	  zsing = zsings;
	  zcosi = zcosis;
	  zsini = zsinis;
	  zcosh = cosq;
	  zsinh = sinq;
	  cc = c1ss;
	  zn = zns;
	  ze = zes;
	  xnoi = 1/xnq;

	  /* Loop breaks when Solar terms are done a second */
	  /* time, after Lunar terms are initialized        */
	  for(;;)
	  {
		/* Solar terms done again after Lunar terms are done */
		a1 = zcosg*zcosh+zsing*zcosi*zsinh;
		a3 = -zsing*zcosh+zcosg*zcosi*zsinh;
		a7 = -zcosg*zsinh+zsing*zcosi*zcosh;
		a8 = zsing*zsini;
		a9 = zsing*zsinh+zcosg*zcosi*zcosh;
		a10 = zcosg*zsini;
		a2 = deep_arg->cosio*a7+ deep_arg->sinio*a8;
		a4 = deep_arg->cosio*a9+ deep_arg->sinio*a10;
		a5 = -deep_arg->sinio*a7+ deep_arg->cosio*a8;
		a6 = -deep_arg->sinio*a9+ deep_arg->cosio*a10;
		x1 = a1*deep_arg->cosg+a2*deep_arg->sing;
		x2 = a3*deep_arg->cosg+a4*deep_arg->sing;
		x3 = -a1*deep_arg->sing+a2*deep_arg->cosg;
		x4 = -a3*deep_arg->sing+a4*deep_arg->cosg;
		x5 = a5*deep_arg->sing;
		x6 = a6*deep_arg->sing;
		x7 = a5*deep_arg->cosg;
		x8 = a6*deep_arg->cosg;
		z31 = 12*x1*x1-3*x3*x3;
		z32 = 24*x1*x2-6*x3*x4;
		z33 = 12*x2*x2-3*x4*x4;
		z1 = 3*(a1*a1+a2*a2)+z31*deep_arg->eosq;
		z2 = 6*(a1*a3+a2*a4)+z32*deep_arg->eosq;
		z3 = 3*(a3*a3+a4*a4)+z33*deep_arg->eosq;
		z11 = -6*a1*a5+deep_arg->eosq*(-24*x1*x7-6*x3*x5);
		z12 = -6*(a1*a6+a3*a5)+ deep_arg->eosq*
		  (-24*(x2*x7+x1*x8)-6*(x3*x6+x4*x5));
		z13 = -6*a3*a6+deep_arg->eosq*(-24*x2*x8-6*x4*x6);
		z21 = 6*a2*a5+deep_arg->eosq*(24*x1*x5-6*x3*x7);
		z22 = 6*(a4*a5+a2*a6)+ deep_arg->eosq*
		  (24*(x2*x5+x1*x6)-6*(x4*x7+x3*x8));
		z23 = 6*a4*a6+deep_arg->eosq*(24*x2*x6-6*x4*x8);
		z1 = z1+z1+deep_arg->betao2*z31;
		z2 = z2+z2+deep_arg->betao2*z32;
		z3 = z3+z3+deep_arg->betao2*z33;
		s3 = cc*xnoi;
		s2 = -0.5*s3/deep_arg->betao;
		s4 = s3*deep_arg->betao;
		s1 = -15*eq*s4;
		s5 = x1*x3+x2*x4;
		s6 = x2*x3+x1*x4;
		s7 = x2*x4-x1*x3;
		se = s1*zn*s5;
		si = s2*zn*(z11+z13);
		sl = -zn*s3*(z1+z3-14-6*deep_arg->eosq);
		sgh = s4*zn*(z31+z33-6);
		sh = -zn*s2*(z21+z23);
		if (xqncl < 5.2359877E-2) sh = 0;
		ee2 = 2*s1*s6;
		e3 = 2*s1*s7;
		xi2 = 2*s2*z12;
		xi3 = 2*s2*(z13-z11);
		xl2 = -2*s3*z2;
		xl3 = -2*s3*(z3-z1);
		xl4 = -2*s3*(-21-9*deep_arg->eosq)*ze;
		xgh2 = 2*s4*z32;
		xgh3 = 2*s4*(z33-z31);
		xgh4 = -18*s4*ze;
		xh2 = -2*s2*z22;
		xh3 = -2*s2*(z23-z21);

		if(isFlagSet(LUNAR_TERMS_DONE_FLAG))
		  break;

		/* Do lunar terms */
		sse = se;
		ssi = si;
		ssl = sl;
		ssh = sh/deep_arg->sinio;
		ssg = sgh-deep_arg->cosio*ssh;
		se2 = ee2;
		si2 = xi2;
		sl2 = xl2;
		sgh2 = xgh2;
		sh2 = xh2;
		se3 = e3;
		si3 = xi3;
		sl3 = xl3;
		sgh3 = xgh3;
		sh3 = xh3;
		sl4 = xl4;
		sgh4 = xgh4;
		zcosg = zcosgl;
		zsing = zsingl;
		zcosi = zcosil;
		zsini = zsinil;
		zcosh = zcoshl*cosq+zsinhl*sinq;
		zsinh = sinq*zcoshl-cosq*zsinhl;
		zn = znl;
		cc = c1l;
		ze = zel;
		SetFlag(LUNAR_TERMS_DONE_FLAG);
	  } /* End of for(;;) */

	  sse = sse+se;
	  ssi = ssi+si;
	  ssl = ssl+sl;
	  ssg = ssg+sgh-deep_arg->cosio/deep_arg->sinio*sh;
	  ssh = ssh+sh/deep_arg->sinio;

	  /* Geopotential resonance initialization for 12 hour orbits */
	  ClearFlag(RESONANCE_FLAG);
	  ClearFlag(SYNCHRONOUS_FLAG);

	  if( !((xnq < 0.0052359877) && (xnq > 0.0034906585)) )
	  {
		if( (xnq < 0.00826) || (xnq > 0.00924) )
		  return;
		if (eq < 0.5) return;
		SetFlag(RESONANCE_FLAG);
		eoc = eq*deep_arg->eosq;
		g201 = -0.306-(eq-0.64)*0.440;
		if (eq <= 0.65)
		{
		  g211 = 3.616-13.247*eq+16.290*deep_arg->eosq;
		  g310 = -19.302+117.390*eq-228.419*
			deep_arg->eosq+156.591*eoc;
		  g322 = -18.9068+109.7927*eq-214.6334*
			deep_arg->eosq+146.5816*eoc;
		  g410 = -41.122+242.694*eq-471.094*
			deep_arg->eosq+313.953*eoc;
		  g422 = -146.407+841.880*eq-1629.014*
			deep_arg->eosq+1083.435*eoc;
		  g520 = -532.114+3017.977*eq-5740*
			deep_arg->eosq+3708.276*eoc;
		}
		else
		{
		  g211 = -72.099+331.819*eq-508.738*
			deep_arg->eosq+266.724*eoc;
		  g310 = -346.844+1582.851*eq-2415.925*
			deep_arg->eosq+1246.113*eoc;
		  g322 = -342.585+1554.908*eq-2366.899*
			deep_arg->eosq+1215.972*eoc;
		  g410 = -1052.797+4758.686*eq-7193.992*
			deep_arg->eosq+3651.957*eoc;
		  g422 = -3581.69+16178.11*eq-24462.77*
			deep_arg->eosq+ 12422.52*eoc;
		  if (eq <= 0.715)
			g520 = 1464.74-4664.75*eq+3763.64*deep_arg->eosq;
		  else
			g520 = -5149.66+29936.92*eq-54087.36*
			  deep_arg->eosq+31324.56*eoc;
		} /* End if (eq <= 0.65) */

		if (eq < 0.7)
		{
		  g533 = -919.2277+4988.61*eq-9064.77*
			deep_arg->eosq+5542.21*eoc;
		  g521 = -822.71072+4568.6173*eq-8491.4146*
			deep_arg->eosq+5337.524*eoc;
		  g532 = -853.666+4690.25*eq-8624.77*
			deep_arg->eosq+ 5341.4*eoc;
		}
		else
		{
		  g533 = -37995.78+161616.52*eq-229838.2*
			deep_arg->eosq+109377.94*eoc;
		  g521 = -51752.104+218913.95*eq-309468.16*
			deep_arg->eosq+146349.42*eoc;
		  g532 = -40023.88+170470.89*eq-242699.48*
			deep_arg->eosq+115605.82*eoc;
		} /* End if (eq <= 0.7) */

		sini2 = deep_arg->sinio*deep_arg->sinio;
		f220 = 0.75*(1+2*deep_arg->cosio+deep_arg->theta2);
		f221 = 1.5*sini2;
		f321 = 1.875*deep_arg->sinio*(1-2*\
			deep_arg->cosio-3*deep_arg->theta2);
		f322 = -1.875*deep_arg->sinio*(1+2*
			deep_arg->cosio-3*deep_arg->theta2);
		f441 = 35*sini2*f220;
		f442 = 39.3750*sini2*sini2;
		f522 = 9.84375*deep_arg->sinio*(sini2*(1-2*deep_arg->cosio-5*
			  deep_arg->theta2)+0.33333333*(-2+4*deep_arg->cosio+
			  6*deep_arg->theta2));
		f523 = deep_arg->sinio*(4.92187512*sini2*(-2-4*
			  deep_arg->cosio+10*deep_arg->theta2)+6.56250012
			  *(1+2*deep_arg->cosio-3*deep_arg->theta2));
		f542 = 29.53125*deep_arg->sinio*(2-8*
			deep_arg->cosio+deep_arg->theta2*
			(-12+8*deep_arg->cosio+10*deep_arg->theta2));
		f543 = 29.53125*deep_arg->sinio*(-2-8*deep_arg->cosio+
			deep_arg->theta2*(12+8*deep_arg->cosio-10*
			deep_arg->theta2));
		xno2 = xnq*xnq;
		ainv2 = aqnv*aqnv;
		temp1 = 3*xno2*ainv2;
		temp = temp1*root22;
		d2201 = temp*f220*g201;
		d2211 = temp*f221*g211;
		temp1 = temp1*aqnv;
		temp = temp1*root32;
		d3210 = temp*f321*g310;
		d3222 = temp*f322*g322;
		temp1 = temp1*aqnv;
		temp = 2*temp1*root44;
		d4410 = temp*f441*g410;
		d4422 = temp*f442*g422;
		temp1 = temp1*aqnv;
		temp = temp1*root52;
		d5220 = temp*f522*g520;
		d5232 = temp*f523*g532;
		temp = 2*temp1*root54;
		d5421 = temp*f542*g521;
		d5433 = temp*f543*g533;
		xlamo = xmao+tle->xnodeo+tle->xnodeo-thgr-thgr;
		bfact = deep_arg->xmdot+deep_arg->xnodot+
		  deep_arg->xnodot-thdt-thdt;
		bfact = bfact+ssl+ssh+ssh;
	  } /* if( !(xnq < 0.0052359877) && (xnq > 0.0034906585) ) */
	  else
	  {
		SetFlag(RESONANCE_FLAG);
		SetFlag(SYNCHRONOUS_FLAG);
		/* Synchronous resonance terms initialization */
		g200 = 1+deep_arg->eosq*(-2.5+0.8125*deep_arg->eosq);
		g310 = 1+2*deep_arg->eosq;
		g300 = 1+deep_arg->eosq*(-6+6.60937*deep_arg->eosq);
		f220 = 0.75*(1+deep_arg->cosio)*(1+deep_arg->cosio);
		f311 = 0.9375*deep_arg->sinio*deep_arg->sinio*
		  (1+3*deep_arg->cosio)-0.75*(1+deep_arg->cosio);
		f330 = 1+deep_arg->cosio;
		f330 = 1.875*f330*f330*f330;
		del1 = 3*xnq*xnq*aqnv*aqnv;
		del2 = 2*del1*f220*g200*q22;
		del3 = 3*del1*f330*g300*q33*aqnv;
		del1 = del1*f311*g310*q31*aqnv;
		fasx2 = 0.13130908;
		fasx4 = 2.8843198;
		fasx6 = 0.37448087;
		xlamo = xmao+tle->xnodeo+tle->omegao-thgr;
		bfact = deep_arg->xmdot+xpidot-thdt;
		bfact = bfact+ssl+ssg+ssh;
	  } /* End if( !(xnq < 0.0052359877) && (xnq > 0.0034906585) ) */

	  xfact = bfact-xnq;

	  /* Initialize integrator */
	  xli = xlamo;
	  xni = xnq;
	  atime = 0;
	  stepp = 720;
	  stepn = -720;
	  step2 = 259200;
	  /* End case dpinit: */
	  return;

	case dpsec: /* Entrance for deep space secular effects */
	  deep_arg->xll = deep_arg->xll+ssl*deep_arg->t;
	  deep_arg->omgadf = deep_arg->omgadf+ssg*deep_arg->t;
	  deep_arg->xnode = deep_arg->xnode+ssh*deep_arg->t;
	  deep_arg->em = tle->eo+sse*deep_arg->t;
	  deep_arg->xinc = tle->xincl+ssi*deep_arg->t;
	  if (deep_arg->xinc < 0)
	  {
		deep_arg->xinc = -deep_arg->xinc;
		deep_arg->xnode = deep_arg->xnode + pi;
		deep_arg->omgadf = deep_arg->omgadf-pi;
	  }
	  if( isFlagClear(RESONANCE_FLAG) )
		return;

	  do
	  {
		if( (atime == 0.0) ||
			((deep_arg->t >= 0) && (atime < 0 )) ||
			((deep_arg->t <  0) && (atime >= 0)) )
		{
		  /* Epoch restart */
		  if( deep_arg->t >= 0 )
			delt = stepp;
		  else
			delt = stepn;

		  atime = 0;
		  xni = xnq;
		  xli = xlamo;
		}
		else
		{
		  if( fabs(deep_arg->t) >= fabs(atime) )
		  {
			if ( deep_arg->t > 0 )
			  delt = stepp;
			else
			  delt = stepn;
		  }
		}

		do
		{
		  if ( fabs(deep_arg->t-atime) >= stepp )
		  {
			SetFlag(DO_LOOP_FLAG);
			ClearFlag(EPOCH_RESTART_FLAG);
		  }
		  else
		  {
			ft = deep_arg->t-atime;
			ClearFlag(DO_LOOP_FLAG);
		  }

		  if( fabs(deep_arg->t) < fabs(atime) )
		  {
			if (deep_arg->t >= 0)
			  delt = stepn;
			else
			  delt = stepp;
			SetFlag(DO_LOOP_FLAG | EPOCH_RESTART_FLAG);
		  }

		  /* Dot terms calculated */
		  if( isFlagSet(SYNCHRONOUS_FLAG) )
		  {
			xndot = del1*sin(xli-fasx2)+del2*sin(2*(xli-fasx4))
			  +del3*sin(3*(xli-fasx6));
			xnddt = del1*cos(xli-fasx2)+2*del2*cos(2*(xli-fasx4))
			  +3*del3*cos(3*(xli-fasx6));
		  }
		  else
		  {
			xomi = omegaq+deep_arg->omgdot*atime;
			x2omi = xomi+xomi;
			x2li = xli+xli;
			xndot = d2201*sin(x2omi+xli-g22)
			  +d2211*sin(xli-g22)
			  +d3210*sin(xomi+xli-g32)
			  +d3222*sin(-xomi+xli-g32)
			  +d4410*sin(x2omi+x2li-g44)
			  +d4422*sin(x2li-g44)
			  +d5220*sin(xomi+xli-g52)
			  +d5232*sin(-xomi+xli-g52)
			  +d5421*sin(xomi+x2li-g54)
			  +d5433*sin(-xomi+x2li-g54);
			xnddt = d2201*cos(x2omi+xli-g22)
			  +d2211*cos(xli-g22)
			  +d3210*cos(xomi+xli-g32)
			  +d3222*cos(-xomi+xli-g32)
			  +d5220*cos(xomi+xli-g52)
			  +d5232*cos(-xomi+xli-g52)
			  +2*(d4410*cos(x2omi+x2li-g44)
				  +d4422*cos(x2li-g44)
				  +d5421*cos(xomi+x2li-g54)
				  +d5433*cos(-xomi+x2li-g54));
		  } /* End of if (isFlagSet(SYNCHRONOUS_FLAG)) */

		  xldot = xni+xfact;
		  xnddt = xnddt*xldot;

		  if(isFlagSet(DO_LOOP_FLAG))
		  {
			xli = xli+xldot*delt+xndot*step2;
			xni = xni+xndot*delt+xnddt*step2;
			atime = atime+delt;
		  }
		}
		while(isFlagSet(DO_LOOP_FLAG) && isFlagClear(EPOCH_RESTART_FLAG));
	  }
	  while(isFlagSet(DO_LOOP_FLAG) && isFlagSet(EPOCH_RESTART_FLAG));

	  deep_arg->xn = xni+xndot*ft+xnddt*ft*ft*0.5;
	  xl = xli+xldot*ft+xndot*ft*ft*0.5;
	  temp = -deep_arg->xnode+thgr+deep_arg->t*thdt;

	  if (isFlagClear(SYNCHRONOUS_FLAG))
		deep_arg->xll = xl+temp+temp;
	  else
		deep_arg->xll = xl-deep_arg->omgadf+temp;

	  return;
	  /*End case dpsec: */

	case dpper: /* Entrance for lunar-solar periodics */
	  sinis = sin(deep_arg->xinc);
	  cosis = cos(deep_arg->xinc);
	  if (fabs(savtsn-deep_arg->t) >= 30)
	  {
		savtsn = deep_arg->t;
		zm = zmos+zns*deep_arg->t;
		zf = zm+2*zes*sin(zm);
		sinzf = sin(zf);
		f2 = 0.5*sinzf*sinzf-0.25;
		f3 = -0.5*sinzf*cos(zf);
		ses = se2*f2+se3*f3;
		sis = si2*f2+si3*f3;
		sls = sl2*f2+sl3*f3+sl4*sinzf;
		sghs = sgh2*f2+sgh3*f3+sgh4*sinzf;
		shs = sh2*f2+sh3*f3;
		zm = zmol+znl*deep_arg->t;
		zf = zm+2*zel*sin(zm);
		sinzf = sin(zf);
		f2 = 0.5*sinzf*sinzf-0.25;
		f3 = -0.5*sinzf*cos(zf);
		sel = ee2*f2+e3*f3;
		sil = xi2*f2+xi3*f3;
		sll = xl2*f2+xl3*f3+xl4*sinzf;
		sghl = xgh2*f2+xgh3*f3+xgh4*sinzf;
		sh1 = xh2*f2+xh3*f3;
		pe = ses+sel;
		pinc = sis+sil;
		pl = sls+sll;
	  }

	  pgh = sghs+sghl;
	  ph = shs+sh1;
	  deep_arg->xinc = deep_arg->xinc+pinc;
	  deep_arg->em = deep_arg->em+pe;

	  if (xqncl >= 0.2)
	  {
		/* Apply periodics directly */
		ph = ph/deep_arg->sinio;
		pgh = pgh-deep_arg->cosio*ph;
		deep_arg->omgadf = deep_arg->omgadf+pgh;
		deep_arg->xnode = deep_arg->xnode+ph;
		deep_arg->xll = deep_arg->xll+pl;
	  }
	  else
	  {
		/* Apply periodics with Lyddane modification */
		sinok = sin(deep_arg->xnode);
		cosok = cos(deep_arg->xnode);
		alfdp = sinis*sinok;
		betdp = sinis*cosok;
		dalf = ph*cosok+pinc*cosis*sinok;
		dbet = -ph*sinok+pinc*cosis*cosok;
		alfdp = alfdp+dalf;
		betdp = betdp+dbet;
		deep_arg->xnode = FMod2p(deep_arg->xnode);
		xls = deep_arg->xll+deep_arg->omgadf+cosis*deep_arg->xnode;
		dls = pl+pgh-pinc*deep_arg->xnode*sinis;
		xls = xls+dls;
		xnoh = deep_arg->xnode;
		deep_arg->xnode = AcTan(alfdp,betdp);

		/* This is a patch to Lyddane modification */
		/* suggested by Rob Matson. */
		if(fabs(xnoh-deep_arg->xnode) > pi)
		{
		  if(deep_arg->xnode < xnoh)
			deep_arg->xnode +=twopi;
		  else
			deep_arg->xnode -=twopi;
		}

		deep_arg->xll = deep_arg->xll+pl;
		deep_arg->omgadf = xls-deep_arg->xll-cos(deep_arg->xinc)*
		  deep_arg->xnode;
	  } /* End case dpper: */
	  return;

  } /* End switch(ientry) */

} /* End of Deep() */

/*------------------------------------------------------------------*/

/* Functions for testing and setting/clearing flags */

/* An int variable holding the single-bit flags */
static int Flags = 0;

  int
isFlagSet(int flag)
{
  return (Flags & flag);
}

  int
isFlagClear(int flag)
{
  return (~Flags & flag);
}

  void
SetFlag(int flag)
{
  Flags |= flag;
}

  void
ClearFlag(int flag)
{
  Flags &= ~flag;
}

/*------------------------------------------------------------------*/

/* Calculates solar position vector */
  void
Calculate_Solar_Position(double time, vector_t *solar_vector)
{
  double mjd,year,T,M,L,e,C,O,Lsa,nu,R,eps;

  mjd = time - 2415020.0;
  year = 1900 + mjd/365.25;
  T = (mjd + Delta_ET(year)/secday)/36525.0;
  M = Radians(Modulus(358.47583 + Modulus(35999.04975*T,360.0)
		- (0.000150 + 0.0000033*T)*Sqr(T),360.0));
  L = Radians(Modulus(279.69668 + Modulus(36000.76892*T,360.0)
		+ 0.0003025*Sqr(T),360.0));
  e = 0.01675104 - (0.0000418 + 0.000000126*T)*T;
  C = Radians((1.919460 - (0.004789 + 0.000014*T)*T)*sin(M)
	  + (0.020094 - 0.000100*T)*sin(2*M) + 0.000293*sin(3*M));
  O = Radians(Modulus(259.18 - 1934.142*T,360.0));
  Lsa = Modulus(L + C - Radians(0.00569 - 0.00479*sin(O)),twopi);
  nu = Modulus(M + C,twopi);
  R = 1.0000002*(1 - Sqr(e))/(1 + e*cos(nu));
  eps = Radians(23.452294 - (0.0130125 + (0.00000164 -
		  0.000000503*T)*T)*T + 0.00256*cos(O));
  R = AU*R;
  solar_vector->x = R*cos(Lsa);
  solar_vector->y = R*sin(Lsa)*cos(eps);
  solar_vector->z = R*sin(Lsa)*sin(eps);
  solar_vector->w = R;
} /*Procedure Calculate_Solar_Position*/

/*------------------------------------------------------------------*/

/* Calculates stellite's eclipse status and depth */
  int
Sat_Eclipsed(vector_t *pos, vector_t *sol, double *depth)
{
  double sd_sun, sd_earth, delta;
  vector_t Rho, earth;

  /* Determine partial eclipse */
  sd_earth = ArcSin(xkmper/pos->w);
  Vec_Sub(sol,pos,&Rho);
  sd_sun = ArcSin(sr/Rho.w);
  Scalar_Multiply(-1,pos,&earth);
  delta = Angle(sol,&earth);
  *depth = sd_earth - sd_sun - delta;
  if( sd_earth < sd_sun )
	return( 0 );
  else
	if( *depth >= 0 )
	  return( 1 );
	else
	  return( 0 );

} /*Function Sat_Eclipsed*/

/*------------------------------------------------------------------*/










#endif


