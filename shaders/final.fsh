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

#define CIE_VERSION 0 //[0 1]

#define MAXIMUM_SAMPLE_COUNT -1 //[-1 100 200 300 400 500 1000 2000 3000 4000 5000 7000 8000 9000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000 65000 70000 75000 80000 85000 90000 95000 100000] Set to -1 for unlimited sample count

#define SHUTTER_SPEED 60.0 //[500.0 250.0 125.0 60.0 30.0 15.0 8.0]
#define ISO 100.0 //[50.0 100.0 200.0 400.0 600.0 800.0 1200.0 1600.0 2400.0 3200.0]
#define F_STOPS 16.0 //[1.4 2.0 2.8 4.0 5.6 8.0 11.0 16.0 22.0 32.0]

//#define USE_SPECTRAL_CAMERA_TONEMAP
#define CAMERA_RESPONSE 4 //[0 1 2 3 4]

//#define DISPLAY_PHASE_FUNCTIONS
//#define DISPLAY_GROUND_ALBEDO

#define GROUND_ALBEDO_R 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define GROUND_ALBEDO_G 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]
#define GROUND_ALBEDO_B 79 //[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255]

layout(location = 0) out vec3 finalColor;

/*
    const int colortex0Format = RGBA32F;
*/

uniform sampler2D colortex0;
uniform sampler2D grainTexture;

uniform sampler3D noise3D;

uniform sampler2D PhaseFunction;

uniform sampler2D phaseTextureRayleigh;

uniform sampler2D phaseTextureAerosol;
uniform sampler2D phaseTextureRainbow;
uniform sampler2D phaseTextureCloud;

uniform sampler2D cameraResponseLUT;

uniform sampler2D usStandardAtmosphere;

uniform sampler3D SunCoords;

uniform float frameTime;
uniform float viewWidth, viewHeight;
uniform int hideGUI;

const bool colortex0Clear  = false;

layout(location = 0) in vec2 textureCoordinate;

#include "/lib/universal/universal.glsl"
#include "/lib/rng/bayer.glsl"

#ifndef USE_SPECTRAL_CAMERA_TONEMAP
    #include "/lib/tonemap/sampledResponseRGB.glsl"
#else
    #include "/lib/tonemap/sampledResponseSpectral.glsl"
#endif

const float K = 16.0;
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

int GetMonth(int day) {
    day = day % 365;
    int month = 0;
    if (day >= 1 && day <= 31) {
        month = 1;
    }
    else if (day > 31 && day <= 59) {
        month = 2;
    }
    else if (day > 59 && day <= 90) {
        month = 3;
    }
    else if (day > 90 && day <= 120) {
        month = 4;
    }
    else if (day > 120 && day <= 151) {
        month = 5;
    }
    else if (day > 151 && day <= 181) {
        month = 6;
    }
    else if (day > 181 && day <= 212) {
        month = 7;
    }
    else if (day > 212 && day <= 243) {
        month = 8;
    }
    else if (day > 243 && day <= 273) {
        month = 9;
    }
    else if (day > 273 && day <= 304) {
        month = 10;
    }
    else if (day > 304 && day <= 334) {
        month = 11;
    }
    else if (day > 334 && day <= 365) {
        month = 12;
    }

    return month;
}

void PrintDebugText(inout vec3 color) {
    int size = viewWidth < 1920 / 2 ? 2 : 3;
    beginText(ivec2(gl_FragCoord.xy / size), ivec2(0, viewHeight / size - (5 + size)));
    text.bgCol = vec4(0.0, 0.0, 0.0, 0.8);
    printString((_S,_a,_m,_p,_l,_e,_s,_colon,_space));
    printInt(int(texture(colortex0, vec2(0.5)).a));

    printLine();

    printString((_U,_n,_l,_i,_m,_i,_t,_e,_d,_space,_S,_a,_m,_p,_l,_e,_s,_colon,_space));
    printBool(MAXIMUM_SAMPLE_COUNT == -1);

    printLine();

    float FPS = float(int(1.0 / frameTime));
    float seconds = texture(colortex0, vec2(0.5)).a / FPS;
    float minutes = seconds / 60.0;
    float hours = minutes / 60.0;

    printString((_T,_i,_m,_e,_space,_E,_l,_a,_s,_p,_e,_d,_colon,_space));
    printInt(int(hours));
    printString((_colon));
    printInt(int(mod(minutes, 60.0)));
    printString((_colon));
    printFloat(mod(seconds, 60.0));

    printLine();

    printString((_C,_I,_E,_space,_V,_e,_r,_s,_i,_o,_n,_colon,_space));
    #if CIE_VERSION == 1
        text.fgCol = vec4(1.0, 1.0, 1.0, 1.0);
        printString((_1,_9,_3,_1));
    #elif CIE_VERSION == 0
        text.fgCol = vec4(1.0, 1.0, 1.0, 1.0);
        printString((_2,_0,_1,_2));
    #endif

    printLine();

    endText(color);
}

void main() {
    float exposure = EV100ToExposure(ComputeEV100()) * rcp(pi);
    vec3 color = texture(colortex0, textureCoordinate).rgb * exposure;
    #ifdef DISPLAY_PHASE_FUNCTIONS
        if(textureCoordinate.y < 0.08) {
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
                    phaseXYZ += SpectrumToXYZExact_CIE2012((texture(phaseTextureRayleigh, lookupCoord).r * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
                }

                D65XYZ += SpectrumToXYZExact_CIE2012(D65Spectrum / (1.0 / 441.0), wavelength) / 441.0;
            }
            vec3 phaseRGB = (phaseXYZ * xyzToRGBMatrix_D65) / (D65XYZ * xyzToRGBMatrix_D65);
            color = phaseRGB * pi;
        }
    #endif
    #ifdef DISPLAY_GROUND_ALBEDO
        if(textureCoordinate.y > 0.9 && textureCoordinate.x > 0.9 && hideGUI == 0) {
            //color = SrgbToLinear(vec3(GROUND_ALBEDO_R, GROUND_ALBEDO_G, GROUND_ALBEDO_B) / 255.0);
            vec3 albedoXYZ = vec3(0.0);
            vec3 D65XYZ = vec3(0.0);
            for(int i = 0; i < 441; ++i) {
                float wavelength = (float(i + 1) / 441.0) * 441.0 + 390.0;

                vec3 groundAlbedoRGB = SrgbToLinear(vec3(GROUND_ALBEDO_R,GROUND_ALBEDO_G,GROUND_ALBEDO_B) / 255.0);
                float groundAlbedo;
                RGBToSpectrum(groundAlbedo, wavelength, groundAlbedoRGB.r, groundAlbedoRGB.g, groundAlbedoRGB.b, 0);
                groundAlbedo = saturate(groundAlbedo);
                if(isnan(groundAlbedo)) {
                    groundAlbedo = 1.0;
                }

                float D65Spectrum = Plancks(6504.0, wavelength);
                if(isnan(D65Spectrum)) {
                    D65Spectrum = 1.0;
                }

                albedoXYZ += SpectrumToXYZExact_CIE2012((groundAlbedo * D65Spectrum) / (1.0 / 441.0), wavelength) / 441.0;
                D65XYZ += SpectrumToXYZExact_CIE2012(D65Spectrum / (1.0 / 441.0), wavelength) / 441.0;
            }

            color = (albedoXYZ * xyzToRGBMatrix_D65) / (D65XYZ * xyzToRGBMatrix_D65);
        }
    #endif
    //finalColor = LinearToSrgb(CameraTonemap(color, ISO));
    finalColor = 1.0 - exp(-color);
    finalColor = LinearToSrgb(color);

    if(hideGUI == 0) {
        PrintDebugText(finalColor);
    }

    //float rayleigh = RayleighCrossSection[int(textureCoordinate.x * 441.0)] / 2.6325e-30;

    //finalColor = textureCoordinate.x > texture(usStandardAtmosphere, textureCoordinate).w / 4.86e18 ? vec3(0.0) : vec3(1.0);
    //finalColor  = vec3(0.0);
    //finalColor += textureCoordinate.x > texture(usStandardAtmosphere, textureCoordinate).x / 1.225 ? vec3(0.0) : vec3(1.0, 0.5, 0.25);
    //finalColor += textureCoordinate.x > texture(usStandardAtmosphere, textureCoordinate).y / 1.225 ? vec3(0.0) : vec3(0.25, 0.5, 1.0);

    finalColor += Bayer64(gl_FragCoord.xy) / 64.0;
}