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

int BinarySearch(int lowIndex, int highIndex, float toFind) {
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
	float altitudeMin = altitude - 0.5, altitudeMax = altitude + 0.5;
	float start = max(altitudeMin, OzoneAltitude[0].x);
	float end = min(altitudeMax, OzoneAltitude[39].x);

	int idx = int(max(distance(0.0, float(BinarySearch(0, 39, start))), 1.0) - 1.0);

	if(end <= start) {
		return 0.0;
	}

	float numberDensity = 0.0;
	for(int entry = idx; entry < 40 && end >= OzoneAltitude[entry].x; ++entry) {
		float a = OzoneAltitude[entry],
		      b = OzoneAltitude[entry + 1],
			  ca = max(a, start),
			  cb = min(b, end),
			  fa = US_OzoneNumberDensity[entry],
			  fb = US_OzoneNumberDensity[entry + 1],
			  invAB = 1.0 / (b - a);

		if(ca >= cb) {
			continue;
		}

		float interpA = mix(fa, fb, (ca - a) * invAB);
		float interpB = mix(fa, fb, (cb - a) * invAB);

		numberDensity += 0.5 * (interpA + interpB) * (cb - ca);
	}
	numberDensity = numberDensity / (altitudeMax - altitudeMin);

	return numberDensity;
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

    imageStore(usStandardAtmosphere_image, fragPos, vec4(us1976LU.layerLocalDensity, max(us1976LU.layerLocalPressure, 0.0), max(us1976LU.layerLocalTemp, 0.0), clamp(O3_NumberDensity(R - planetRadius), 0.0, 4.86e18)));
}