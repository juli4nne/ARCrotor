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
#include <poll.h>
#include "sgp4sdp4.h"
#include "rotor.h"

struct Server serverConstructor(int domain, int service, int protocol, unsigned long interface, int port, int backlog){
	struct Server server;

	server.domain = domain;
	server.service = service;
	server.protocol = protocol;
	server.interface = interface;
	server.port = port;
	server.backlog = backlog;

	server.address.sin_family = domain;
	server.address.sin_port = htons(port);
	server.address.sin_addr.s_addr = inet_addr("127.0.0.1");

	server.socket = socket(domain, service, protocol);
	if(server.socket == 0){
		perror("Failed to connect socket.\n");
		exit(1);
	}

	if(bind(server.socket, (struct sockaddr*)&server.address, sizeof(server.address)) < 0){
		perror("Failed to bind socket.\n");
		exit(1);
	}

	if(listen(server.socket, server.backlog) < 0){
		perror("Failed to listen to socket.\n");
		exit(1);
	}

	return server;
}

char* readHTML(char* filename, int type){
    //write in HTML file
    FILE *fp = fopen(filename, "r");;
	if(fp == NULL){
		printf("Error opening file.\n");
		exit(1);
	}
	fseek(fp, 0, SEEK_END);
	int fileSize = ftell(fp);
	rewind(fp);
	char* html = (char*)malloc(fileSize + 1);
	fread(html, 1, fileSize, fp);
	html[fileSize] = '\0';
	fclose(fp);

    char header[300];
    switch(type){
        case 1: strcpy(header, "HTTP/1.1 200 OK\nContent-Type: text/html\nContent-Length: ");
        break;

        case 2: 
        sprintf(header, "HTTP/1.1 302 Found\nLocation: %s\nContent-Type: text/html\nContent-Length: ", filename);
        break;

        case 3: strcpy(header, "HTTP/1.1 200 OK\nContent-Type: text/javascript\nContent-Length: ");
        break;
    }

    char contentLength[6];
    sprintf(contentLength, "%d", fileSize);
    char* web = (char*)malloc((strlen(header) + strlen(contentLength) + fileSize + 3)*sizeof(char));
    strcpy(web, header);
    strcat(web, contentLength);
    strcat(web, "\n\n");
    strcat(web, html);
    free(html);

    return web;
}

int main() {
    struct Server server = serverConstructor(AF_INET, SOCK_STREAM, 0, INADDR_ANY, 1234, 10);

    printf("Waiting for connection at 127.0.0.1:1234\n");
    int addrlen = sizeof(server.address);

    char buffer[1000];

    int newSocket = accept(server.socket, (struct sockaddr *)&server.address, (socklen_t*)&addrlen);
    int bytesRead = recv(newSocket, buffer, sizeof(buffer) - 1, 0);
    int load1 = 0;
    int load2 = 0;
    char* web = (char*)malloc(100000*sizeof(char));
    int bytesSent = 0;
    while(load1 == 0 || load2 == 0){
        if(strncmp(buffer, "GET / ", 6) == 0){
            strcpy(web, readHTML("geo_input.html",1));
            bytesSent = send(newSocket, web, strlen(web), 0);
            load1 = 1;
            recv(newSocket, buffer, sizeof(buffer) - 1, 0);
        } else if(strncmp(buffer, "GET /jquery.js", 14) == 0){
            strcpy(web, readHTML("jquery.js",3));
            bytesSent = send(newSocket, web, strlen(web), 0);
            load2 = 1;
            recv(newSocket, buffer, sizeof(buffer) - 1, 0);
        }
    }
    

    struct pollfd serverPoll[1];
    serverPoll[0].fd = newSocket;
    serverPoll[0].events = POLLIN;

    double lat, lon, height;

    while(1){
        long pollStatus = poll(serverPoll, 1, 500);
        if(pollStatus < 0){
            perror("Poll Failed");
            break;
        }
        if(serverPoll[0].revents & POLLIN){
            strcpy(buffer, "\0");
            int newSocket = accept(server.socket, (struct sockaddr *)&server.address, (socklen_t*)&addrlen);
            int bytesRead = recv(newSocket, buffer, sizeof(buffer) - 1, 0);
            if (bytesRead < 0){
                perror("Error reading form.");
                break;
            }
            buffer[bytesRead] = '\0';

            if(strstr(buffer, "POST") != NULL){
                char* temp = (char*)malloc(16*sizeof(char));
                char* latData = strstr(buffer, "lat=");
                if(latData != NULL){
                    latData += 4;
                    for(int i = 0; i < 15; i++){
                        if(isdigit(latData[i]) || latData[i] == '.' || latData[i] == '-'){
                            temp[i] = latData[i];
                        } else{
                            temp[i] = '\0';
                            break;
                        }
                    }
                    temp[15] = '\0';
                    lat = atof(temp);
                }
                char* lonData = strstr(buffer, "lon=");
                if(lonData != NULL){
                    lonData += 4;
                    for(int i = 0; i < 15; i++){
                        if(isdigit(lonData[i]) || lonData[i] == '.'|| lonData[i] == '-'){
                            temp[i] = lonData[i];
                        } else{
                            temp[i] = '\0';
                            break;
                        }
                    }
                    temp[15] = '\0';
                    lon = atof(temp);
                }
                char* heiData = strstr(buffer, "hei=");
                if(heiData != NULL){
                    heiData += 4;
                    for(int i = 0; i < 15; i++){
                        if(isdigit(heiData[i]) || heiData[i] == '.' || heiData[i] == '-'){
                            temp[i] = heiData[i];
                        } else{
                            temp[i] = '\0';
                            break;
                        }
                    }
                    temp[15] = '\0';
                    height = atof(temp);
                }


                printf("\nSubmitted!\n");
                break;
            }
        }
    }

    // tle input buffers
    char tleStr[139] = "\0";
    char temp1[70] = "\0";
    char temp2[70] = "\0";
    
    // newSocket = accept(server.socket, (struct sockaddr *)&server.address, (socklen_t*)&addrlen);
    strcpy(web, readHTML("tle_input.html",2));
	bytesSent = send(newSocket, web, strlen(web), 0);


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
	tsince,            // Time since epoch (in minutes) 
	jul_epoch,         // Julian date of epoch          
	jul_utc,           // Julian UTC date               
	eclipse_depth = 0, // Depth of satellite eclipse    
	// Satellite's observed position, range, range rate 
	sat_azi, sat_ele, sat_range, sat_range_rate,
	// Satellites geodetic position and velocity 
	sat_lat, sat_lon, sat_alt, sat_vel,
	// Solar azimuth and elevation 
	sun_azi, sun_ele;
    // Used for storing function return codes 
    int flg;
    char
	ephem[5],       // Ephemeris in use string  
	sat_status[12]; // Satellite eclipse status 
	
    while(1){
        long pollStatus = poll(serverPoll, 1, 500);
        if(pollStatus < 0){
            perror("Poll Failed");
            break;
        }
        if(serverPoll[0].revents & POLLIN){
            char buffer[1000];
            int newSocket = accept(server.socket, (struct sockaddr *)&server.address, (socklen_t*)&addrlen);
            int bytesRead = recv(newSocket, buffer, sizeof(buffer) - 1, 0);
            if (bytesRead < 0){
                perror("Error reading form.");
                break;
            }
            buffer[bytesRead] = '\0';

            if(strstr(buffer, "POST") != NULL){
                char* line1Data = strstr(buffer, "line1=");
                if(line1Data != NULL){
                    line1Data += 6;
                    for(int i = 0; i < 69; i++){
                        if(line1Data[i] == '+'){
                            temp1[i] = ' ';
                        } else {
                            temp1[i] = line1Data[i];
                        }
                    }
                }
                char* line2Data = strstr(buffer, "line2=");
                if(line2Data != NULL){
                    line2Data += 6;
                    for(int i = 0; i < 69; i++){
                        if(line2Data[i] == '+'){
                            temp2[i] = ' ';
                        } else {
                            temp2[i] = line2Data[i];
                        }
                    }
                }
                strncat(tleStr, temp1, 69);
                strncat(tleStr, temp2, 69);
                tleStr[138] = '\0';

                printf("\nSubmitted!\n");
                break;
            } else {
                printf("\nNo Submission."); 
            }
            
        }
    }

    flg = Input_Tle_Set(tleStr, &tle);

    // Print satellite name and TLE read status 
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

    
    do  // Loop 
  {
	// Get UTC calendar and convert to Julian 
	UTC_Calendar_Now(&utc, &tv);
	jul_utc = Julian_Date(&utc, &tv);

	// Convert satellite's epoch time to Julian  
	// and calculate time since epoch in minutes 
	jul_epoch = Julian_Date_of_Epoch(tle.epoch);
	tsince = (jul_utc - jul_epoch) * xmnpda;

	// Copy the ephemeris type in use to ephem string 
	if( isFlagSet(DEEP_SPACE_EPHEM_FLAG) )
	  strcpy(ephem,"SDP4");
	else
	  strcpy(ephem,"SGP4");

	// Call NORAD routines according to deep-space flag 
	if( isFlagSet(DEEP_SPACE_EPHEM_FLAG) )
	  SDP4(tsince, &tle, &pos, &vel);
	else
	  SGP4(tsince, &tle, &pos, &vel);

	// Scale position and velocity vectors to km and km/sec 
	Convert_Sat_State( &pos, &vel );

	// Calculate velocity of satellite 
	Magnitude( &vel );
	sat_vel = vel.w;

	// All angles in rads. Distance in km. Velocity in km/s 
	// Calculate satellite Azi, Ele, Range and Range-rate 
	Calculate_Obs(jul_utc, &pos, &vel, &obs_geodetic, &obs_set);

	// Calculate satellite Lat North, Lon East and Alt. 
	Calculate_LatLonAlt(jul_utc, &pos, &sat_geodetic);

	// Calculate solar position and satellite eclipse depth 
	// Also set or clear the satellite eclipsed flag accordingly 
	Calculate_Solar_Position(jul_utc, &solar_vector);
	Calculate_Obs(jul_utc,&solar_vector,&zero_vector,&obs_geodetic,&solar_set);

	if( Sat_Eclipsed(&pos, &solar_vector, &eclipse_depth) )
	  SetFlag( SAT_ECLIPSED_FLAG );
	else
	  ClearFlag( SAT_ECLIPSED_FLAG );

	// Copy a satellite eclipse status string in sat_status 
	if( isFlagSet( SAT_ECLIPSED_FLAG ) )
	  strcpy( sat_status, "Eclipsed" );
	else
	  strcpy( sat_status, "In Sunlight" );

	// Convert and print satellite and solar data 
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
  }  // End of do 
  while( 1 );






    printf("\ndone\n");

    return 0;

}
