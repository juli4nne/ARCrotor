#include "sgp4sdp4/sgp4sdp4.h"

char tle_line_1[70];
char tle_line_2[70];

int main() {
    printf("Enter Line 1:\n");
    scanf ("%[^\n]%*c", tle_line_1);
    printf("Enter Line 2:\n");
    scanf ("%[^\n]%*c", tle_line_2);
    
    return 0;

}
