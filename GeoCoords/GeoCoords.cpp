#include <iostream>
#include <cmath>
#include <cctype>
#include <cstring>
#include <cstdlib>
using namespace std;

class GeoCoord
{
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

    static constexpr double pi = 3.1415926535897932384626433832795;
    static double deg2rad(double deg) { return deg * pi / 180.0; }
    static double rad2deg(double rad) { return rad * 180.0 / pi; }

public:
    typedef struct BoundingBox
    {
        double minLat, maxLat;
        double minLon, maxLon;

        // Check if a point is inside this bounding box
        bool contains(const GeoCoord& point) const {
            return (point.latitude >= minLat && point.latitude <= maxLat &&
                point.longitude >= minLon && point.longitude <= maxLon);
        }

        // Check if this bounding box overlaps another
        bool overlaps(const BoundingBox& other) const {
            bool latOverlap = (minLat <= other.maxLat && maxLat >= other.minLat);
            bool lonOverlap = (minLon <= other.maxLon && maxLon >= other.minLon);
            return latOverlap && lonOverlap;
        }


    } BoundingBox;

    GeoCoord();
    GeoCoord(double, double);
    GeoCoord(double, double, double, double, double, double, char, char);
    void SetDMS();
    void SetFloat();
    void normalizeLongitude(double);
    void clampLatitude(double);
    BoundingBox CalcBoundingBox(double);
    GeoCoord interpolateTo(const GeoCoord& other, double f);
    double bearingTo(const GeoCoord& other);
    double distanceTo(const GeoCoord& other);
    void printConversion();
    
    friend ostream& operator<<(ostream& out, const GeoCoord& c);
    friend istream& operator>>(istream& in, GeoCoord& c);
};

// Bounding box around this point
GeoCoord::BoundingBox GeoCoord::CalcBoundingBox(double radius)
{
    //const double degLatKm = 111.0;
    const double degLatMi = 69.0;
    double deltaLat = radius / degLatMi;
    double deltaLon = radius / (degLatMi * cos(latitude * pi / 180.0));

    BoundingBox box;
    box.minLat = latitude - deltaLat;
    box.maxLat = latitude + deltaLat;
    box.minLon = longitude - deltaLon;
    box.maxLon = longitude + deltaLon;
    return box;

}

/****************************************
The Haversine formula itself gives you an angular distance (in radians) between two points on a sphere.
To turn that into a real-world distance, you multiply by the radius of the sphere you’re modeling.
- If you use Earth’s mean radius = 6,371 km, the result will be in kilometers.
- If you use Earth’s radius in miles ≈ 3,959 mi, the result will be in miles.
- If you use meters (6,371,000 m), the result will be in meters.
***************************************/
double GeoCoord::distanceTo(const GeoCoord& other)
{
    constexpr double R = 3959.0; // Earth radius mi
    double phi1 = deg2rad(latitude);
    double phi2 = deg2rad(other.latitude);
    double dPhi = deg2rad(other.latitude - latitude);
    double dLambda = deg2rad(other.longitude - longitude);

    double a = sin(dPhi / 2) * sin(dPhi / 2) +
        cos(phi1) * cos(phi2) *
        sin(dLambda / 2) * sin(dLambda / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;

}

// Bearing to another coordinate
double GeoCoord::bearingTo(const GeoCoord& other)
{
    double phi1 = deg2rad(latitude);
    double phi2 = deg2rad(other.latitude);
    double dLambda = deg2rad(other.longitude - longitude);

    double y = sin(dLambda) * cos(phi2);
    double x = cos(phi1) * sin(phi2) -
        sin(phi1) * cos(phi2) * cos(dLambda);

    double theta = atan2(y, x);
    return fmod((rad2deg(theta) + 360.0), 360.0);

}

// Interpolate toward another coordinate at fraction f
GeoCoord GeoCoord::interpolateTo(const GeoCoord& other, double f)
{
    double phi1 = deg2rad(latitude), lambda1 = deg2rad(longitude);
    double phi2 = deg2rad(other.latitude), lambda2 = deg2rad(other.longitude);

    double d = 2 * asin(sqrt(pow(sin((phi2 - phi1) / 2), 2) +
        cos(phi1) * cos(phi2) * pow(sin((lambda2 - lambda1) / 2), 2)));

    double A = sin((1 - f) * d) / sin(d);
    double B = sin(f * d) / sin(d);

    double x = A * cos(phi1) * cos(lambda1) + B * cos(phi2) * cos(lambda2);
    double y = A * cos(phi1) * sin(lambda1) + B * cos(phi2) * sin(lambda2);
    double z = A * sin(phi1) + B * sin(phi2);

    double phi3 = atan2(z, sqrt(x * x + y * y));
    double lambda3 = atan2(y, x);

    return GeoCoord(rad2deg(phi3), rad2deg(lambda3));

}

// Normalizing and clamping are useful for error correction when crossing the international dateline
void GeoCoord::normalizeLongitude(double lval)
{
    while (lval > 180.0) longitude -= 360.0;
    while (lval < -180.0) longitude += 360.0;
}

void GeoCoord::clampLatitude(double lval)
{
    if (lval > 90.0) latitude = 90.0;
    if (lval < -90.0) latitude = -90.0;
}


GeoCoord::GeoCoord()
{
    latitude = longitude = latdegrees = londegrees = latminutes = lonminutes = latseconds = lonseconds = 0;
    latDir = lonDir = '\0';
}

GeoCoord::GeoCoord(double latitude, double longitude)
{
    this->latitude = latitude;
    this->longitude = longitude;
    latdegrees = londegrees = latminutes = lonminutes = latseconds = lonseconds = 0;
    latDir = lonDir = '\0';
}

GeoCoord::GeoCoord(double latdegrees, double londegrees, double latminutes, double lonminutes, double latseconds, double lonseconds, char latDir, char lonDir)
{
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

void GeoCoord::SetFloat()
{
    latitude = latdegrees + (latminutes / 60) + (latseconds / 3600);
    longitude = londegrees + (lonminutes / 60) + (lonseconds / 3600);
    if (toupper(latDir) == 'S') {
        latitude = -latitude;
    }
    if (toupper(lonDir) == 'W') {
        longitude = -longitude;
    }
}

void GeoCoord::SetDMS()
{
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

// Overload >> for input
istream& operator>>(istream& in, GeoCoord& c)
{
    cout << "Please input the latitude and longitude coordinate pair: ";
    in >> c.latitude >> c.longitude;
    return in;
}

// Overload << for output
ostream& operator<<(ostream& out, const GeoCoord& c)
{
    out << "(" << c.latitude << ", " << c.longitude << ")";
    return out;
}

// Overload << for output
ostream& operator<<(ostream& out, const GeoCoord::BoundingBox& c)
{
    GeoCoord min = GeoCoord(c.minLat, c.minLon);
    GeoCoord max = GeoCoord(c.maxLat, c.maxLon);
    out << min << " to " << max;
    return out;
}

void GeoCoord::printConversion()
{
    cout << this << endl
        << latdegrees << "°" << latminutes << "'" << latseconds << "\"" << latDir << " by "
        << londegrees << "°" << lonminutes << "'" << lonseconds << "\"" << lonDir << endl;
}


int main(int argc, char* argv[])
{
    if (argc == 2)
    {
        if (strcmp(argv[1], "-float") == 0)
        {
            GeoCoord g;
            cin >> g;
            g.SetDMS();
            cout << g;
        }
        else if (strcmp(argv[1], "-dms") == 0)
        {
            double latdeg, londeg, latmin, lonmin, latsec, lonsec;
            char lator, lonor;
            cout << "Please input the latitude degrees, minutes, seconds, and orientation (N, S): ";
            cin >> latdeg >> latmin >> latsec >> lator;
            cout << "Please input the longitude degrees, minutes, seconds, and orientation (E, W): ";
            cin >> londeg >> lonmin >> lonsec >> lonor;
            GeoCoord g = GeoCoord(latdeg, londeg, latmin, lonmin, latsec, lonsec, lator, lonor);
            g.SetFloat();
            cout << g;

        }
        else if (strcmp(argv[1], "-toolkit") == 0)
        {
            GeoCoord one, two;
            double r = 30;
            cin >> one >> two;
            auto box1 = one.CalcBoundingBox(r);
            auto box2 = two.CalcBoundingBox(r);
            cout << "The distance between " << one << " and " << two << " is " << one.distanceTo(two) << " miles." << endl;
            cout << "The midpoint between " << one << " and " << two << " is " << one.interpolateTo(two, .5) << endl;
            cout << "A bounding box around " << one << " having a radius of " << r << " miles is " << box1 << endl;
            cout << "A bounding box around " << two << " having a radius of " << r << " miles is " << box2 << endl;
            if (box1.overlaps(box2))
            {
                cout << "These bounding boxes overlap!" << endl;
            }
            else
            {
                cout << "These bounding boxes DO NOT overlap!" << endl;
            }
            GeoCoord otherPoint = two.interpolateTo(one, .33);
            if (box1.contains(otherPoint))
            {
                cout << "The point one third of the  way out between points one and two is within the first's bounding box!" << endl;
            }
            else
            {
                cout << "The point one third of the  way out between points one and two is NOT within the first's bounding box!" << endl;
            }
            cout << "From " << one << " facing True North, you would rotate " << one.bearingTo(two) << "° clockwise to bear toward " << two << ".\r\n";
        }
    }
    return EXIT_SUCCESS;
}
