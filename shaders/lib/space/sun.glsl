#if !defined LIB_SPACE_SUN
#define LIB_SPACE_SUN
    float PhysicalSun96(in vec3 sceneDirection, in vec3 lightDirection, in float wavelength, in float radiance) {
        float cosTheta = dot(sceneDirection, lightDirection);
        float centerToEdge = clamp(acos(cosTheta) / sunAngularRadius, 0.0, 1.0);
        float cosSunRadius = cos(sunAngularRadius);
        if(cosTheta < cosSunRadius) return 0.0;

        float fit = -0.023 + 0.292 * pow(((wavelength - 390.0) / 441.0), -1.0);
        if(((wavelength - 390.0) / 441.0) < 0.63) {
            fit = 0.4 * pow(0.2 + ((wavelength - 390.0) / 441.0), -0.5);
        }

        float mu = sqrt(1.0 - (centerToEdge*centerToEdge));
        float factor = 1.0 - 1.0 * (1.0 - pow(mu, fit));

        return radiance * factor;
    }

    float PhysicalSun07(in vec3 sceneDirection, in vec3 lightDirection, in float wavelength, in float radiance) {
        /*
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_u.pro
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_v.pro
        */
        const float au = -8.9829751;
        const float bu = 0.0069093916;
        const float cu = -1.8144591e-6;
        const float du = 2.2540875e-10;
        const float eu = -1.3389747e-14;
        const float fu = 3.0453572e-19;

        const float av = 9.2891180;
        const float bv = -0.0062212632;
        const float cv = 1.5788029e-6;
        const float dv = -1.9359644e-10;
        const float ev = 1.1444469e-14;
        const float fv = -2.599494e-19;

        float cosTheta = dot(sceneDirection, lightDirection);
        float cosSunRadius = cos(sunAngularRadius);
        if(cosTheta < cosSunRadius) return 0.0;

        float lambda = wavelength * 10.0; //Angstroms
        float u = au + (bu * lambda) + (cu * pow(lambda, 2.0)) + (du * pow(lambda, 3.0)) + (eu * pow(lambda, 4.0)) + (fu * pow(lambda, 5.0));
        float v = av + (bv * lambda) + (cv * pow(lambda, 2.0)) + (dv * pow(lambda, 3.0)) + (ev * pow(lambda, 4.0)) + (fv * pow(lambda, 5.0));

        float centerToEdge = sin(acos(cosTheta)) / sin(sunAngularRadius);
        float darkening = 1.0 - u - v + u*cos(asin(centerToEdge)) + v*pow(cos(asin(centerToEdge)), 2.0);

        return radiance * darkening;
    }

    float PhysicalSun(in vec3 sceneDirection, in vec3 lightDirection, in float wavelength, in float radiance) {
        float cosTheta = dot(sceneDirection, lightDirection);
        float cosSunRadius = cos(sunAngularRadius);
        if(cosTheta < cosSunRadius) return 0.0;

        return PhysicalSun07(sceneDirection, lightDirection, wavelength, radiance);
    }

    float CorrectDegreeRange(in float d) {
        d = mod(d, 360.0);
        if(d < 0.0) {
            d += 360.0;
        }
        return d;
    }

    vec3 CalculateSunPosition(in float JD, in float longitude, in float latitude) {
        /*
            Equations from Astronomical Algorithms by Jean Meeus, second edition.
        */
        float T = (JD - J2000) / 36525.0; //Eq 12.1
        float siderealTime = 280.46061837 + 360.98564736629 * (JD - J2000) + 0.000387933 * T*T - pow(T, 3.0) / 38710000.0; //Eq 12.4
        float geometricMeanLongitude = CorrectDegreeRange(280.46646 + 36000.76983 * T + 0.0003032 * T*T); //Eq 25.2
        float meanAnomaly = CorrectDegreeRange(357.52911 + 35999.05029 * T + 0.0001537 * T*T); //Eq 25.3

        float eccentricity = 0.016708634 - 0.000042037 * T - 0.0000001267 * T*T; //Eq 25.4

        float center = +(1.914602 - 0.004817 * T - 0.000014 * T*T) * sin(meanAnomaly) + (0.019993 - 0.000101 * T) * sin(2.0 * meanAnomaly) + 0.000289 * sin(3.0 * meanAnomaly);

        float trueLongitude = geometricMeanLongitude + center;
        float trueAnomaly = meanAnomaly + center;

        float R = (1.000001018 * (1.0 - square(eccentricity))) / (1.0 + eccentricity * cos(trueAnomaly)); //Eq 25.5

        float tan_rightAcension = (cos(eccentricity) * sin(trueLongitude)) / cos(trueLongitude); //Eq 25.6
        float sin_declination = sin(eccentricity) * sin(trueLongitude); //Eq 25.7

        float localHourAngle = siderealTime - longitude - atan(tan_rightAcension);

        float tan_azimuth = sin(localHourAngle) / (cos(localHourAngle) * sin(latitude) - tan(asin(sin_declination)) * cos(latitude)); //Eq 13.5
        float sin_altitude = sin(latitude) * sin_declination + cos(latitude) * cos(asin(sin_declination)) * cos(localHourAngle); //Eq 13.6

        return vec3(atan(tan_azimuth), asin(sin_altitude), R);
    }
#endif