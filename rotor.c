#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
// change this out for linux
#include <windows.h>


#define M_PI 3.14159265358979323846
double downlink_frequency = 1000000.00;

char tle_line_1[70];
char tle_line_2[70];

char satNum1[15] = "\0";
char intDesignator[15] = "\0";
double epochYear;
double firstDer;
char secondDer[15] = "\0";
char dragTerm[15] = "\0";
int ephemeris;
int elementNum;
int checkSum1;

char satNum2[15] = "\0";
double inclination;
double rightAscension;
int eccentricity;
double argPerigee;
double meanAnomaly;
double meanMotion;
int revNum;
int checkSum2;

int checkSumChecker(char* line){
    int lineProgress = 0;
    int bufferProgress = 0;
    int buffer[68] = {0};
    while(lineProgress != 68){
        if(isdigit(line[lineProgress])){
            buffer[bufferProgress] = line[lineProgress] - '0';
            bufferProgress++;
        } else if(line[lineProgress] == '-'){
            buffer[bufferProgress] = 1;
            bufferProgress++;
        }
        lineProgress++;
    }
    int checkSum = line[lineProgress]- '0';
    int value;
    for(int i = 0; i < bufferProgress; i++){
        value += buffer[i];
    }
    value = value % 10;

    if(checkSum == value){
        return 1;
    }
    return 0;
}

int doThis() {

    printf("Enter Line 1:\n");
    scanf ("%[^\n]%*c", tle_line_1);
    printf("Enter Line 2:\n");
    scanf ("%[^\n]%*c", tle_line_2);

    //reading line 1
    int progress = 2;
    int wordCount = 0;
    char ch = '\0';
    char lookAhead = '\0';
    while(1){
        char buffer[20] = {0};
        int len = 0;
        ch = tle_line_1[progress];
        if(ch == '\0'){
            break;
        }
        if(ch == ' '){
            progress++;
            continue;
        }
        lookAhead = ch;
        while(lookAhead != '\0' && lookAhead != ' '){
            buffer[len] = lookAhead;
            len++;
            progress++;
            lookAhead = tle_line_1[progress];
        }
        progress++;

        switch(wordCount){
            case 0: strcpy(satNum1, buffer); wordCount++; continue;
            case 1: strcpy(intDesignator, buffer); wordCount++; continue;
            case 2: epochYear = atof(buffer); wordCount++; continue;
            case 3: firstDer = atof(buffer); wordCount++; continue;
            case 4: strcpy(secondDer, buffer); wordCount++; continue;
            case 5: strcpy(dragTerm, buffer); wordCount++; continue;
            case 6: ephemeris = atoi(buffer); wordCount++; continue;
            case 7: elementNum = atoi(buffer); wordCount++; continue;
            default: break;
        }
    }

    // reading line 2
    progress = 2;
    wordCount = 0;
    ch = '\0';
    lookAhead = '\0';
    int lineLen = strlen(tle_line_2) - 6;
    int pointPassed = 0;
    while(1){
        char buffer[20] = {0};
        int len = 0;
        if(pointPassed == 0){
            ch = tle_line_2[progress];
        } else {
            progress--;
            ch = tle_line_2[progress];
        }
        if(ch == '\0'){
            break;
        }
        if(ch == ' '){
            progress++;
            continue;
        }
        lookAhead = ch;
        while(lookAhead != '\0' && lookAhead != ' '){
            if((progress != lineLen) || pointPassed){
                buffer[len] = lookAhead;
                len++;
                progress++;
                lookAhead = tle_line_2[progress];
            } else {
                pointPassed = 1;
                break;
            }
            
        }
        progress++;

        switch(wordCount){
            case 0: strcpy(satNum2, buffer); wordCount++; continue;
            case 1: inclination = atof(buffer); wordCount++; continue;
            case 2: rightAscension = atof(buffer); wordCount++; continue;
            case 3: eccentricity = atoi(buffer); wordCount++; continue;
            case 4: argPerigee = atof(buffer); wordCount++; continue;
            case 5: meanAnomaly = atof(buffer); wordCount++; continue;
            case 6: meanMotion = atof(buffer); wordCount++; continue;
            case 7: revNum = atoi(buffer); wordCount++; continue;
            default: break;
        }
    }

    checkSum1 = elementNum % 10;
    elementNum = elementNum / 10;
    checkSum2 = revNum % 10;
    revNum = revNum / 10;

    int test = checkSumChecker(tle_line_1);
    if(test){
        printf("\n\n%s\n%s\n%f\n%f\n%s\n%s\n%d\n%d\n%d\n\n", satNum1, intDesignator, epochYear, firstDer, secondDer, dragTerm, ephemeris, elementNum, checkSum1);
        printf("%s\n%f\n%f\n%d\n%f\n%f\n%f\n%d\n%d\n", satNum2, inclination, rightAscension, eccentricity, argPerigee, meanAnomaly, meanMotion, revNum, checkSum2);
        printf("done\n\n");
    } else {
        printf("Invalid Checksum!");
    }

    
    return 0;

	/* //satellite representation
	predict_orbital_elements_t *orbital_elements;
	orbital_elements = predict_parse_tle(tle_line_1, tle_line_2);

	//containers for time-dependent properties
	struct predict_position orbit; //independent of an observer
	struct predict_observation observation; //relative to an observer

	//prediction time
	predict_julian_date_t pred_time;
	pred_time = predict_to_julian(time(NULL));

	//calculate point along satellite orbit
	predict_orbit(orbital_elements, &orbit, pred_time);
	printf("%f %f\n", orbit.latitude, orbit.longitude);

	//observe the orbit
	double latitude = 63.42;
	double longitude = 10.39;
	predict_observer_t *observer;
	observer = predict_create_observer("LA1K", latitude/180.0*M_PI, longitude/180.0*M_PI, 0);
	predict_observe_orbit(observer, &orbit, &observation);
	printf("%f %f\n", observation.azimuth, observation.elevation);

	//calculate time for the next time satellite passes over the horizon
	struct predict_observation aos = predict_next_aos(observer, orbital_elements, pred_time);
	double time_to_aos = (aos.time - pred_time)*24*60;
	printf("Time to next AOS: %d minutes\n", time_to_aos);

	//calculate doppler shift
	double doppler_shift = predict_doppler_shift(&observation, downlink_frequency); */
}

int main(){
    Sleep(1000);
    for(int i = 0; i < 10; i++){
        doThis();
    }
}