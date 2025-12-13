// GeoCoords.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <cctype>
#include <cstring>
#include <cstdlib>
using namespace std;

class GeoCoord {
    double latitude;
    double longitude;
    double latdegrees;
    double londegrees;
    double latminutes;
    double lonminutes;
    double latseconds;
    double lonseconds;
    char latDir;
    char lonDir;
public:
    GeoCoord();
    GeoCoord(double, double);
    GeoCoord(double, double, double, double, double, double, char, char);
    void SetDMS();
    void SetFloat();
    friend ostream& operator<<(ostream& out, const GeoCoord& c);
};

GeoCoord::GeoCoord() {
    latitude = longitude = latdegrees = londegrees = latminutes = lonminutes = latseconds = lonseconds = 0;
    latDir = lonDir = '\0';
}

GeoCoord::GeoCoord(double latitude, double longitude) {
    this->latitude = latitude;
    this->longitude = longitude;
    latdegrees = londegrees = latminutes = lonminutes = latseconds = lonseconds = 0;
    latDir = lonDir = '\0';
}

GeoCoord::GeoCoord(double latdegrees, double londegrees, double latminutes, double lonminutes, double latseconds, double lonseconds, char latDir, char lonDir) {
    latitude = 0;
    longitude = 0;
    this->latdegrees = latdegrees;
    this->londegrees = londegrees;
    this->latminutes = latminutes;
    this->lonminutes = lonminutes;
    this->latseconds = latseconds;
    this->lonseconds = lonseconds;
    this->latDir = latDir;
    this->lonDir = lonDir;
}

void GeoCoord::SetFloat() {
    latitude = latdegrees + (latminutes / 60) + (latseconds / 3600);
    longitude = londegrees + (lonminutes / 60) + (lonseconds / 3600);
    if (toupper(latDir) == 'S') {
        latitude = -latitude;
    }
    if (toupper(lonDir) == 'W') {
        longitude = -longitude;
    }
}

void GeoCoord::SetDMS() {
    double latint, latfloat, lonint, lonfloat;

    double abslat = abs(latitude);
    double abslon = abs(longitude);
    latfloat = modf(abslat, &latint);
    lonfloat = modf(abslon, &lonint);
    latdegrees = latint;
    londegrees = lonint;

    double minlatint, minlatfloat, minlonint, minlonfloat;
    latminutes = latfloat * 60;
    lonminutes = lonfloat * 60;
    minlatfloat = modf(latminutes, &minlatint);
    minlonfloat = modf(lonminutes, &minlonint);

    latseconds = minlatfloat * 60;
    lonseconds = minlonfloat * 60;

    if (latitude < 0) {
        latDir = 'S';
    }
    else {
        latDir = 'N';
    }

    if (longitude < 0) {
        lonDir = 'W';
    }
    else {
        lonDir = 'E';
    }
}

// Overload << for output
ostream& operator<<(ostream& out, const GeoCoord& c)
{
    out << "(" << c.latitude << ", " << c.longitude << ")" << endl
        << c.latdegrees << "° " << c.latminutes << "' " << c.latseconds << "\" " << c.latDir << " by "
        << c.londegrees << "° " << c.lonminutes << "' " << c.lonseconds << "\" " << c.lonDir << endl;
    return out;
}


int main(int argc, char* argv[])
{
    if (argc == 2) {
        if (strcmp(argv[1], "-float") == 0) {
            double lat;
            double lon;
            cout << "Please input the latitude and longitude coordinate pair: ";
            cin >> lat >> lon;
            GeoCoord g = GeoCoord(lat, lon);
            g.SetDMS();
            cout << g;
        }
        else if (strcmp(argv[1], "-dms") == 0) {
            double latdeg;
            double londeg;
            double latmin;
            double lonmin;
            double latsec;
            double lonsec;
            char lator;
            char lonor;
            cout << "Please input the latitude degrees, minutes, seconds, and orientation (N, S): ";
            cin >> latdeg >> latmin >> latsec >> lator;
            cout << "Please input the longitude degrees, minutes, seconds, and orientation (E, W): ";
            cin >> londeg >> lonmin >> lonsec >> lonor;
            GeoCoord g = GeoCoord(latdeg, londeg, latmin, lonmin, latsec, lonsec, lator, lonor);
            g.SetFloat();
            cout << g;
        }
    }
    return EXIT_SUCCESS;
}
