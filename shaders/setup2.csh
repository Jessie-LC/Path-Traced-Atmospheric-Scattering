#version 450

/*
    Copyright (C) 2023-2025  Jessie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

layout (local_size_x = 16, local_size_y = 16) in;
const ivec3 workGroups = ivec3(32, 32, 1);

layout (rgba32f) uniform image2D usStandardAtmosphere_image;

uniform sampler3D noise3D;

#include "/lib/universal/universal.glsl"
#include "/lib/rng/pcg.glsl"

#include "/lib/atmosphere/constants.glsl"
#include "/lib/atmosphere/usStandardAtmosphere1976.glsl"

void US_StandardAtmosphereLookupUVReverse(vec2 coord, out float R) {
	float uvR  = RemoveUvMargin(coord.y, 512);

	const float H = sqrt(lutRadiusSquared - (planetRadius * planetRadius));

	float rho = H * uvR;
	R = sqrt(rho * rho + (planetRadius * planetRadius));
}

const float OzoneAltitude[40] = {
    0.0,
    500.0,
    1000.0,
    2000.0,
    4000.0,
    6000.0,
    8000.0,
    10000.0,
    12000.0,
    14000.0,
    16000.0,
    18000.0,
    20000.0,
    22000.0,
    24000.0,
    26000.0,
    28000.0,
    30000.0,
    32000.0,
    34000.0,
    36000.0,
    38000.0,
    40000.0,
    42000.0,
    44000.0,
    46000.0,
    48000.0,
    50000.0,
    52000.0,
    54000.0,
    56000.0,
    58000.0,
    60000.0,
    62000.0,
    64000.0,
    66000.0,
    68000.0,
    70000.0,
    72000.0,
    74000.0
};

const float US_OzoneNumberDensity[40] = {
    7.6e17, //The US Standard Atmosphere model from 1976 doesn't have a sea level value for ozone number density, so I set it based on how graphs of the density look. The following two density values are for more stable interpolation
    7.4e17,
    7.2e17,
    6.8e17,
    5.8e17,
    5.7e17,
    6.5e17,
    1.13e18,
    2.02e18,
    2.35e18,
    2.95e18,
    4.04e18,
    4.77e18,
    4.86e18,
    4.54e18,
    4.03e18,
    3.24e18,
    2.52e18,
    2.03e18,
    1.58e18,
    1.22e18,
    8.73e17,
    6.07e17,
    3.98e17,
    2.74e17,
    1.69e17,
    1.03e17,
    6.64e16,
    3.84e16,
    2.55e16,
    1.61e16,
    1.12e16,
    7.33e15,
    4.81e15,
    3.17e15,
    1.72e15,
    7.5e14,
    5.4e14,
    2.2e14,
    1.7e14
};

int BinarySearch_Ozone(int lowIndex, int highIndex, float toFind) {
    while (lowIndex < highIndex) {
        int midIndex = (lowIndex + highIndex) >> 1;
        if (OzoneAltitude[midIndex] < toFind) {
            lowIndex = midIndex + 1;
        } else if (OzoneAltitude[midIndex] > toFind) {
            highIndex = midIndex;
        } else {
            return midIndex;
        }
    }
    return highIndex - 1;
}

float O3_NumberDensity(in float altitude) {
	float altitudeMin = altitude;
	float start = max(altitudeMin, OzoneAltitude[0]);

    if (altitude > OzoneAltitude[39]) {
        return 0.0;
    }

    //*
	int k = BinarySearch_Ozone(0, 39, start);

    float ozoneK   =             US_OzoneNumberDensity[k];
    float ozoneK1  =         US_OzoneNumberDensity[k + 1];
    float ozoneK2  =         US_OzoneNumberDensity[k + 2];
    float ozoneKn1 = US_OzoneNumberDensity[max(k - 1, 0)];
    float altitudeK   =             OzoneAltitude[k];
    float altitudeK1  =         OzoneAltitude[k + 1];
    float altitudeK2  =         OzoneAltitude[k + 2];
    float altitudeKn1 = OzoneAltitude[max(k - 1, 0)];
    float mk = CardinalSplineM(ozoneK1, ozoneKn1, altitudeK1, altitudeKn1, 0.0);
    float mk1 = CardinalSplineM(ozoneK2, ozoneK, altitudeK2, altitudeK, 0.0);
    return CubicHermiteSpline(ozoneK, ozoneK1, mk, mk1, altitude, altitudeK, altitudeK1);
    /*/
	int idx = BinarySearch_Ozone(0, 39, start);

    return Remap(
                              altitude, 
                    OzoneAltitude[idx], 
                OzoneAltitude[idx + 1], 
            US_OzoneNumberDensity[idx], 
        US_OzoneNumberDensity[idx + 1]
    );
    //*/
}

const float AerosolAltitude[51] = {
    0.0,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    11.0,
    12.0,
    13.0,
    14.0,
    15.0,
    16.0,
    17.0,
    18.0,
    19.0,
    20.0,
    21.0,
    22.0,
    23.0,
    24.0,
    25.0,
    26.0,
    27.0,
    28.0,
    29.0,
    30.0,
    31.0,
    32.0,
    33.0,
    34.0,
    35.0,
    36.0,
    37.0,
    38.0,
    39.0,
    40.0,
    41.0,
    42.0,
    43.0,
    44.0,
    45.0,
    46.0,
    47.0,
    48.0,
    49.0,
    50.0,
};

const float AerosolCoefficient_Altitude[51] = {
    1.58e-1,
    6.95e-2,
    3.00e-2,
    1.26e-2,
    6.66e-3,
    5.02e-3,
    3.54e-3,
    3.29e-3,
    3.39e-3,
    3.25e-3,
    3.17e-3,
    2.97e-3,
    3.12e-3,
    2.88e-3,
    2.82e-3,
    2.65e-3,
    2.52e-3,
    2.49e-3,
    2.41e-3,
    2.03e-3,
    1.49e-3,
    1.08e-3,
    8.13e-4,
    6.22e-4,
    4.93e-4,
    4.15e-4,
    3.62e-4,
    2.77e-4,
    2.12e-4,
    1.63e-4,
    1.25e-4,
    9.55e-5,
    7.31e-5,
    5.60e-5,
    4.29e-5,
    3.29e-5,
    2.52e-5,
    1.93e-5,
    1.48e-5,
    1.13e-5,
    8.66e-6,
    6.64e-6,
    5.08e-6,
    3.89e-6,
    2.98e-6,
    2.28e-6,
    1.75e-6,
    1.34e-6,
    1.03e-6,
    7.86e-7,
    6.02e-7
};

int BinarySearch_Aerosol(int lowIndex, int highIndex, float toFind) {
    while (lowIndex < highIndex) {
        int midIndex = (lowIndex + highIndex) >> 1;
        if (AerosolAltitude[midIndex] < toFind) {
            lowIndex = midIndex + 1;
        } else if (AerosolAltitude[midIndex] > toFind) {
            highIndex = midIndex;
        } else {
            return midIndex;
        }
    }
    return highIndex - 1;
}

float AerosolDensity(in float altitude) {
	float altitudeMin = altitude;
	float start = max(altitudeMin, AerosolAltitude[0]);

    if (altitude > AerosolAltitude[50]) {
        return exp(-altitude / 1.2);
    }

    //*
	int k = BinarySearch_Aerosol(0, 50, start);

    float aerosolK   =             AerosolCoefficient_Altitude[k];
    float aerosolK1  =         AerosolCoefficient_Altitude[k + 1];
    float aerosolK2  =         AerosolCoefficient_Altitude[k + 2];
    float aerosolKn1 = AerosolCoefficient_Altitude[max(k - 1, 0)];
    float altitudeK   =             AerosolAltitude[k];
    float altitudeK1  =         AerosolAltitude[k + 1];
    float altitudeK2  =         AerosolAltitude[k + 2];
    float altitudeKn1 = AerosolAltitude[max(k - 1, 0)];
    float mk = CardinalSplineM(aerosolK1, aerosolKn1, altitudeK1, altitudeKn1, 0.0);
    float mk1 = CardinalSplineM(aerosolK2, aerosolK, altitudeK2, altitudeK, 0.0);
    return CubicHermiteSpline(aerosolK, aerosolK1, mk, mk1, altitude, altitudeK, altitudeK1) * 6.32911392405;
    /*/
	int idx = BinarySearch_Aerosol(0, 50, start);

    return Remap(
                                    altitude, 
                        AerosolAltitude[idx], 
                    AerosolAltitude[idx + 1], 
            AerosolCoefficient_Altitude[idx], 
        AerosolCoefficient_Altitude[idx + 1]
    ) * 6.32911392405; // Normalization constant calculated from the Preetham aerosol coefficient function with a turbidity of 1.44
    //*/
}

void main() {
    ivec2 fragPos = ivec2(gl_GlobalInvocationID.xy);

    float R;
    US_StandardAtmosphereLookupUVReverse(vec2(fragPos) / 512.0, R);

    usStandardAtmosphere1976LU us1976LU = usStandardAtmosphere1976Lookup(R - planetRadius);
    usStandardAtmosphere1976LU seaLevel = usStandardAtmosphere1976Lookup(0.0);

    if((R - planetRadius) > 85e3) {
        us1976LU.layerLocalDensity = 1.225 * exp(-(R - planetRadius) / 10e3);
    }

    imageStore(usStandardAtmosphere_image, fragPos, vec4(us1976LU.layerLocalDensity, AerosolDensity((R - planetRadius) / 1000.0), max(us1976LU.layerLocalTemp, 0.0), clamp(O3_NumberDensity(R - planetRadius), 0.0, 4.86e18)));
}