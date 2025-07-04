#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include "sgp4sdp4.h"

int main() {

	int server = socket(AF_INET, SOCK_STREAM, 0);
	struct sockaddr_in server_addr;
	server_addr.sin_family = AF_INET;
	server_addr.sin_port = htons(80);
	server_addr.sin_addr.s_addr = INADDR_ANY;
	bind(server, (struct sockaddr *)&server_addr, sizeof(server_addr));
	listen(server, 5);
	int client = accept(server, NULL, NULL);

    // tle input
    char tleStr[139] = "\0";
    char temp1[70] = "\0";
    char temp2[70] = "\0";
    printf("Enter Line 1:\n");
    scanf("%[^\n]%*c", temp1);
    printf("Enter Line 2:\n");
    scanf("%[^\n]%*c", temp2);
    strncat(tleStr, temp1, 69);
    strncat(tleStr, temp2, 69);
    tleStr[138] = '\0';
    
    // coords input
    double lat, lon, height;
    printf("\nEnter Latitude, Longitude, and Height above Sea Level (in fractional degrees and kilometers):\n");
    scanf("%lf %lf %lf", &lat, &lon, &height);

    // user geodetic
    geodetic_t obs_geodetic = {lat, lon, height, 0.0};
    // zero vector
    vector_t zero_vector = {0,0,0,0};
    // sat vel and pos (initialization)
    vector_t vel = zero_vector;
    vector_t pos = zero_vector;
    // azimuth, elevation, range, range rate
    vector_t obs_set;
    // solar eci
    vector_t solar_vector = zero_vector;
    // solar observed azimuth and elevation 
    vector_t solar_set;
    // date and time
    struct tm utc;
    struct timeval tv;
    // sat geodetic
    geodetic_t sat_geodetic;
    // tle struct
    tle_t tle;

    // from sgp4sdp4 main.c
    double
	tsince,            /* Time since epoch (in minutes) */
	jul_epoch,         /* Julian date of epoch          */
	jul_utc,           /* Julian UTC date               */
	eclipse_depth = 0, /* Depth of satellite eclipse    */
	/* Satellite's observed position, range, range rate */
	sat_azi, sat_ele, sat_range, sat_range_rate,
	/* Satellites geodetic position and velocity */
	sat_lat, sat_lon, sat_alt, sat_vel,
	/* Solar azimuth and elevation */
	sun_azi, sun_ele;
    /* Used for storing function return codes */
    int flg;
    char
	ephem[5],       /* Ephemeris in use string  */
	sat_status[12]; /* Satellite eclipse status */

    flg = Input_Tle_Set(tleStr, &tle);

    /* Print satellite name and TLE read status */
    printf(" %s: ", tle.sat_name);
    if( flg == -2 )
    {
	    printf("TLE set bad - Exiting!\n");
	    exit(-2);
    }
    else {
	    printf("TLE set good - Happy Tracking!\n");
    }
    ClearFlag(ALL_FLAGS);
    select_ephemeris(&tle);

    do  /* Loop */
  {
	/* Get UTC calendar and convert to Julian */
	UTC_Calendar_Now(&utc, &tv);
	jul_utc = Julian_Date(&utc, &tv);

	/* Convert satellite's epoch time to Julian  */
	/* and calculate time since epoch in minutes */
	jul_epoch = Julian_Date_of_Epoch(tle.epoch);
	tsince = (jul_utc - jul_epoch) * xmnpda;

	/* Copy the ephemeris type in use to ephem string */
	if( isFlagSet(DEEP_SPACE_EPHEM_FLAG) )
	  strcpy(ephem,"SDP4");
	else
	  strcpy(ephem,"SGP4");

	/* Call NORAD routines according to deep-space flag */
	if( isFlagSet(DEEP_SPACE_EPHEM_FLAG) )
	  SDP4(tsince, &tle, &pos, &vel);
	else
	  SGP4(tsince, &tle, &pos, &vel);

	/* Scale position and velocity vectors to km and km/sec */
	Convert_Sat_State( &pos, &vel );

	/* Calculate velocity of satellite */
	Magnitude( &vel );
	sat_vel = vel.w;

	/** All angles in rads. Distance in km. Velocity in km/s **/
	/* Calculate satellite Azi, Ele, Range and Range-rate */
	Calculate_Obs(jul_utc, &pos, &vel, &obs_geodetic, &obs_set);

	/* Calculate satellite Lat North, Lon East and Alt. */
	Calculate_LatLonAlt(jul_utc, &pos, &sat_geodetic);

	/* Calculate solar position and satellite eclipse depth */
	/* Also set or clear the satellite eclipsed flag accordingly */
	Calculate_Solar_Position(jul_utc, &solar_vector);
	Calculate_Obs(jul_utc,&solar_vector,&zero_vector,&obs_geodetic,&solar_set);

	if( Sat_Eclipsed(&pos, &solar_vector, &eclipse_depth) )
	  SetFlag( SAT_ECLIPSED_FLAG );
	else
	  ClearFlag( SAT_ECLIPSED_FLAG );

	/* Copy a satellite eclipse status string in sat_status */
	if( isFlagSet( SAT_ECLIPSED_FLAG ) )
	  strcpy( sat_status, "Eclipsed" );
	else
	  strcpy( sat_status, "In Sunlight" );

	/* Convert and print satellite and solar data */
	sat_azi = Degrees(obs_set.x);
	sat_ele = Degrees(obs_set.y);
	sat_range = obs_set.z;
	sat_range_rate = obs_set.w;
	sat_lat = Degrees(sat_geodetic.lat);
	sat_lon = Degrees(sat_geodetic.lon);
	sat_alt = sat_geodetic.alt;

	sun_azi = Degrees(solar_set.x);
	sun_ele = Degrees(solar_set.y);

	printf("\n Date: %02d/%02d/%04d UTC: %02d:%02d:%02d  Ephemeris: %s"
		"\n Azi=%6.1f Ele=%6.1f Range=%8.1f Range Rate=%6.2f"
		"\n Lat=%6.1f Lon=%6.1f  Alt=%8.1f  Vel=%8.3f"
		"\n Stellite Status: %s - Depth: %2.3f"
		"\n Sun Azi=%6.1f Sun Ele=%6.1f\n",
		utc.tm_mday, utc.tm_mon, utc.tm_year,
		utc.tm_hour, utc.tm_min, utc.tm_sec, ephem,
		sat_azi, sat_ele, sat_range, sat_range_rate,
		sat_lat, sat_lon, sat_alt, sat_vel,
		sat_status, eclipse_depth,
		sun_azi, sun_ele);

	usleep(1e5);
  }  /* End of do */
  while( 1 );






    printf("\ndone\n");

    return 0;

}
