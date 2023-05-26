#version 450

/*
    Copyright (C) 2023  Jessie

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

#define SHUTTER_SPEED 60.0 //[500.0 250.0 125.0 60.0 30.0 15.0 8.0]
#define ISO 100.0 //[50.0 100.0 200.0 400.0 600.0 800.0 1200.0 1600.0 2400.0 3200.0]
#define F_STOPS 16.0 //[1.4 2.0 2.8 4.0 5.6 8.0 11.0 16.0 22.0 32.0]

#define CAMERA_RESPONSE 4 //[0 1 2 3 4]

//#define DISPLAY_PHASE_FUNCTIONS

layout(location = 0) out vec3 finalColor;

/*
    const int colortex0Format = RGBA32F;
*/

uniform sampler2D colortex0;

uniform sampler3D noise3D;

uniform sampler1D phaseTextureRayleigh;

uniform sampler2D phaseTextureAerosol;
uniform sampler2D phaseTextureRainbow;
uniform sampler2D phaseTextureCloud;

uniform sampler2D cameraResponseLUT;

uniform sampler2D usStandardAtmosphere;

layout(location = 0) in vec2 textureCoordinate;

#include "/lib/universal/universal.glsl"
#include "/lib/rng/bayer.glsl"

#include "/lib/tonemap/sampledResponse.glsl"

const float K = 40.0;
const float calibration = 2.0 * K / 100.0;
float ComputeEV100() {
    return log2(square(F_STOPS) / (1.0 / SHUTTER_SPEED) * 100.0 / ISO);
}
float ComputeISO(float aperture, float shutterSpeed, float ev) {
    return (square(aperture) * 100.0) / (shutterSpeed * pow(2.0, ev));
}
float EV100ToExposure(in float EV100) {
    float maxLuminance = calibration * pow(2.0, EV100);

    return 1.0 / maxLuminance;
}

void main() {
    float exposure = EV100ToExposure(ComputeEV100()) * rcp(pi);
    vec3 color = texture(colortex0, textureCoordinate).rgb * exposure;
    #ifdef DISPLAY_PHASE_FUNCTIONS
        vec3 phaseXYZ = vec3(0.0);
        vec3 D65XYZ = vec3(0.0);
        for(int i = 0; i < 441; ++i) {
            float wavelength = (float(i + 1) / 441.0) * 441.0 + 390.0;

            float D65Spectrum = Plancks(6504.0, wavelength);

            vec2 lookupCoord = vec2(textureCoordinate.x, (wavelength - 390.0) / 441.0);
            if(textureCoordinate.y < 0.08 && textureCoordinate.y > 0.06) {
                phaseXYZ += SpectrumToXYZExact_CIE2012((texture(phaseTextureAerosol, lookupCoord).r * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
            } else if(textureCoordinate.y < 0.06 && textureCoordinate.y > 0.04) {
                phaseXYZ += SpectrumToXYZExact_CIE2012((texture(phaseTextureCloud, lookupCoord).r * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
            } else if(textureCoordinate.y < 0.04 && textureCoordinate.y > 0.02) {
                phaseXYZ += SpectrumToXYZExact_CIE2012((texture(phaseTextureRainbow, lookupCoord).r * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
            } else if(textureCoordinate.y < 0.02) {
                phaseXYZ += SpectrumToXYZExact_CIE2012((texture(phaseTextureRayleigh, textureCoordinate.x).r * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
            }

            D65XYZ += SpectrumToXYZExact_CIE2012(D65Spectrum / (1.0 / 441.0), wavelength) / 441.0;
        }
        vec3 phaseRGB = (phaseXYZ * xyzToRGBMatrix_D65) / (D65XYZ * xyzToRGBMatrix_D65);
        if(textureCoordinate.y < 0.08) {
            color = phaseRGB * pi;
        }
    #endif
    finalColor = LinearToSrgb(CameraTonemap(color, ISO));

    //finalColor = texture(usStandardAtmosphere, textureCoordinate).rrr / 1.0;

    finalColor += Bayer64(gl_FragCoord.xy) / 64.0;
}