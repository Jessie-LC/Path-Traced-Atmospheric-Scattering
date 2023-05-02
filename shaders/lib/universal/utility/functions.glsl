float square(in float x) {
    return _square(x);
}
int square(in int x) {
    return _square(x);
}
vec2 square(in vec2 x) {
    return _square(x);
}
vec3 square(in vec3 x) {
    return _square(x);
}
vec4 square(in vec4 x) {
    return _square(x);
}

float cube(in float x) {
    return _cube(x);
}
int cube(in int x) {
    return _cube(x);
}
vec2 cube(in vec2 x) {
    return _cube(x);
}
vec3 cube(in vec3 x) {
    return _cube(x);
}
vec4 cube(in vec4 x) {
    return _cube(x);
}

float pow5(in float x) {
    return _pow5(x);
}
int pow5(in int x) {
    return _pow5(x);
}
vec2 pow5(in vec2 x) {
    return _pow5(x);
}
vec3 pow5(in vec3 x) {
    return _pow5(x);
}
vec4 pow5(in vec4 x) {
    return _pow5(x);
}

float saturate(in float x) {
    return _saturate(x);
}
int saturate(in int x) {
    return _saturateInt(x);
}
vec2 saturate(in vec2 x) {
    return _saturate(x);
}
vec3 saturate(in vec3 x) {
    return _saturate(x);
}
vec4 saturate(in vec4 x) {
    return _saturate(x);
}

float minof(vec2 x) { 
    return min(x.x, x.y); 
}
float minof(vec3 x) { 
    return min(min(x.x, x.y), x.z); 
}
float minof(vec4 x) { 
    x.xy = min(x.xy, x.zw); return min(x.x, x.y); 
}

float maxof(vec2 x) { 
    return max(x.x, x.y); 
}
float maxof(vec3 x) { 
    return max(max(x.x, x.y), x.z); 
}
float maxof(vec4 x) { 
    x.xy = max(x.xy, x.zw); return max(x.x, x.y); 
}

float rcp(in int x) {
    float y = float(x);
    return _rcp(y);
}
float rcp(in float x) {
    return _rcp(x);
}
vec2 rcp(in vec2 x) {
    return _rcp(x);
}
vec3 rcp(in vec3 x) {
    return _rcp(x);
}
vec4 rcp(in vec4 x) {
    return _rcp(x);
}

float log10(in float x) {
    return _log10(x, 10.0);
}

int log10(in int x) {
    return int(_log10(x, 10.0));
}

vec2 log10(in vec2 x) {
    return _log10(x, 10.0);
}

vec3 log10(in vec3 x) {
    return _log10(x, 10.0);
}

vec4 log10(in vec4 x) {
    return _log10(x, 10.0);
}

float linearstep(in float x, float low, float high) {
    float data = x;
    float mapped = (data-low)/(high-low);

    return saturate(mapped);
}

vec2 linearstep(in vec2 x, float low, float high) {
    vec2 data = x;
    vec2 mapped = (data-low)/(high-low);

    return saturate(mapped);
}

vec3 linearstep(in vec3 x, float low, float high) {
    vec3 data = x;
    vec3 mapped = (data-low)/(high-low);

    return saturate(mapped);
}

vec4 linearstep(in vec4 x, float low, float high) {
    vec4 data = x;
    vec4 mapped = (data-low)/(high-low);

    return saturate(mapped);
}

vec2 sincos(float x) { return vec2(sin(x), cos(x)); }

vec2 GetEquirectangularCoord(in vec3 rayDirection) {
    float lon = atan(rayDirection.z, rayDirection.x);
    if(rayDirection.z < 0) {
        lon = 2 * pi - atan(-rayDirection.z, rayDirection.x);
    }

    float lat = acos(rayDirection.y);

    const vec2 rads = vec2(1.0 / (pi * 2.0), 1.0 / pi);
    vec2 sphereCoords = vec2(lon, lat) * rads;

    return sphereCoords;
}

vec2 RSI(vec3 pos, vec3 dir, float radius) {
    float radiusSquared = square(radius);
    float posDotDir = dot(pos, dir);
    float endDist = posDotDir*posDotDir + radiusSquared - dot(pos, pos);

    if(endDist < 0.0) return vec2(-1.0);

    endDist = sqrt(endDist);
    vec2 ret = -posDotDir + vec2(-endDist, endDist);

    return ret;
}

bool IntersectSphere(in vec3 rayPosition, in vec3 rayDirection, in float radius, out float dist) {
    vec2 sphereDists = RSI(rayPosition, rayDirection, radius);
    
    if(sphereDists.x > 0.001) {
        dist = sphereDists.x;
        return true;
    } else {
        dist = sphereDists.y;
        if(dist < 0.001) {
            return false;
        } else {
            return true;
        }
    }
}

int LinePlaneIntersection(
    vec3 Lo, vec3 Lv,
    vec3 Po, vec3 Pn,
    out float t
) {
    float NdotV = dot(Pn, Lv);
    float NdotO = dot(Pn, Lo - Po);
    if (NdotV == 0.0) {
        if (NdotO == 0.0) {
            return 2;
        } else {
            return 0;
        }
    }

    t = NdotO / -NdotV;

    return 1;
}

int RayPlaneIntersection(
    vec3 Ro, vec3 Rv,
    vec3 Po, vec3 Pn,
    out float t
) {
    int tmp = LinePlaneIntersection(Ro, Rv, Po, Pn, t);
    if (tmp == 1 && t < 0.0) {
        return 0;
    }
    
    if (dot(Rv, Pn) > 0.0) { return 0; }

    return tmp;
}

int LineSphereIntersect(
	vec3 Lo, vec3 Lv,
	vec3 So, float Sr,
	out float t0, out float t1
) {
	float Rsquared = dot(Lo - So, Lo - So);
	float RMu      = dot(Lv, Lo - So) / length(Lv);

	float discriminant = RMu * RMu - Rsquared + Sr * Sr;

	if (discriminant < 0.0) {
		return 0;
	}

	if (discriminant == 0.0) {
		t0 = -RMu;
		return 1;
	}

	t0 = (-RMu - sqrt(discriminant)) / length(Lv);
	t1 = (-RMu + sqrt(discriminant)) / length(Lv);

	return 2;
}

vec3 GenerateUnitVector(vec2 hash) {
    hash.x *= tau; hash.y = hash.y * 2.0 - 1.0;
    return vec3(vec2(sin(hash.x), cos(hash.x)) * sqrt(saturate(1.0 - hash.y * hash.y)), hash.y);
}

vec3 GenerateCosineVector(vec3 vector, vec2 xy) {
    vec3 dir = GenerateUnitVector(xy);

    if(vector + dir == vec3(0.0)) {
        return vector;
    } else {
        return fNormalize(vector + dir);
    }
}

vec3 GenerateConeVector(vec3 vector, vec2 xy, float angle) {
    xy.x *= radians(360.0);
    float cosAngle = cos(angle);
    xy.y = xy.y * (1.0 - cosAngle) + cosAngle;
    vec3 sphereCap = vec3(vec2(cos(xy.x), sin(xy.x)) * sqrt(1.0 - xy.y * xy.y), xy.y);
    return Rotate(sphereCap, vec3(0, 0, 1), vector);
}

float Plancks(float t, float lambda) {
    const float h = 6.62607015e-16;
    const float c = 2.9e17;
    const float k = 1.38e-5;

    float p1 = 2.0 * h * pow(c, 2.0) * pow(lambda, -5.0);
    float p2 = exp((h * c) / (lambda * k * t)) - 1.0;

    return p1 / p2;
}

uniform sampler2D CIELUT;
uniform sampler2D colorToSpectrumLUT;

vec3 SpectrumToXYZExact_CIE2012(in float spectrum, in float w) {
    float n = (w - 390.0);
    int i = int(n);

    int n0 = int(n);
    int n1 = min(n0 + 1, 440);

    vec3 cie0 = texelFetch(CIELUT, ivec2(n0, 1), 0).xyz;
    vec3 cie1 = texelFetch(CIELUT, ivec2(n1, 1), 0).xyz;

    vec3 xyz = mix(cie0, cie1, fract(n));

    xyz = vec3(spectrum * 683.368) * xyz;

    return xyz * float(441.0) / 113.042;
}

float InterpolateSpectrum(float wavelength, int lutRow) {
	if (wavelength < 380.000000 || wavelength > 720.000000) {
		return 0.0;
	}

	float fractionalIndex = float(textureSize(colorToSpectrumLUT, 0).x - 1) * (wavelength - 380.0) / (720.0 - 380.0);
	int indexBelow = int(floor(fractionalIndex));
	int indexAbove = int(ceil (fractionalIndex));
	return mix(
		texelFetch(colorToSpectrumLUT, ivec2(indexBelow, lutRow), 0).x,
		texelFetch(colorToSpectrumLUT, ivec2(indexAbove, lutRow), 0).x,
		fract(fractionalIndex)
	);
}

//This code is based directly on code from Mitsuba
/*
    Copyright (c) 2017 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
void RGBToSpectrum(out float result, in float wl, in float r, in float g, in float b, in int intent) {
    result = 0.0;

    if(intent == 0) {
        float rgbToSpecWhite = InterpolateSpectrum(wl, 0);
        float rgbToSpecCyan = InterpolateSpectrum(wl, 1);
        float rgbToSpecMagenta = InterpolateSpectrum(wl, 2);
        float rgbToSpecYellow = InterpolateSpectrum(wl, 3);
        float rgbToSpecRed = InterpolateSpectrum(wl, 4);
        float rgbToSpecGreen = InterpolateSpectrum(wl, 5);
        float rgbToSpecBlue = InterpolateSpectrum(wl, 6);
        if (r <= g && r <= b) {
            result += r * rgbToSpecWhite;

            if (g <= b) {
                result += (g - r) * rgbToSpecCyan;
                result += (b - g) * rgbToSpecBlue;
            } else {
                result += (b - r) * rgbToSpecCyan;
                result += (g - b) * rgbToSpecGreen;
            }
        } else if (g <= r && g <= b) {
            result += g * rgbToSpecWhite;

            if (r <= b) {
                result += (r - g) * rgbToSpecMagenta;
                result += (b - r) * rgbToSpecBlue;
            } else {
                result += (b - g) * rgbToSpecMagenta;
                result += (r - b) * rgbToSpecRed;
            }
        } else {
            result += b * rgbToSpecWhite;

            if (r <= g) {
                result += (r - b) * rgbToSpecYellow;
                result += (g - r) * rgbToSpecGreen;
            } else {
                result += (g - b) * rgbToSpecYellow;
                result += (r - g) * rgbToSpecRed;
            }
        }

        result *= .94;
    } else {
        float rgbToIllumWhite = InterpolateSpectrum(wl, 7);
        float rgbToIllumCyan = InterpolateSpectrum(wl, 8);
        float rgbToIllumMagenta = InterpolateSpectrum(wl, 9);
        float rgbToIllumYellow = InterpolateSpectrum(wl, 10);
        float rgbToIllumRed = InterpolateSpectrum(wl, 11);
        float rgbToIllumGreen = InterpolateSpectrum(wl, 12);
        float rgbToIllumBlue = InterpolateSpectrum(wl, 13);
        if (r <= g && r <= b) {
            result += r * rgbToIllumWhite;

            if (g <= b) {
                result += (g - r) * rgbToIllumCyan;
                result += (b - g) * rgbToIllumBlue;
            } else {
                result += (b - r) * rgbToIllumCyan;
                result += (g - b) * rgbToIllumGreen;
            }
        } else if (g <= r && g <= b) {
            result += g * rgbToIllumWhite;

            if (r <= b) {
                result += (r - g) * rgbToIllumMagenta;
                result += (b - r) * rgbToIllumBlue;
            } else {
                result += (b - g) * rgbToIllumMagenta;
                result += (r - b) * rgbToIllumRed;
            }
        } else {
            result += b * rgbToIllumWhite;

            if (r <= g) {
                result += (r - b) * rgbToIllumYellow;
                result += (g - r) * rgbToIllumGreen;
            } else {
                result += (g - b) * rgbToIllumYellow;
                result += (r - g) * rgbToIllumRed;
            }
        }

        result *= .86445;
    }
    
    if(isnan(result)) {
        result = 1.0;
    }
}
//End of Mistuba code

float Air(in float wavelength) {
    return 1.0+8.06051E-5+2.480990E-2/(132.274-pow(wavelength,-2.0))+1.74557E-4/(39.32957-pow(wavelength,-2.0));
}
float Water(in float wavelength) {
    const float n1 = 1.779e-4;
    const float n2 = -1.05e-6;
    const float n3 = 1.6e-8;
    const float n4 = -2.02e-6;
    const float n5 = 15.868;
    const float n6 = 0.01155;
    const float n7 = -0.00423;
    const float n8 = -4382.0;
    const float n9 = 1.1455e6;
    const float T = 20.0;
    const float S = 0.0;
    return 1.31405 + (n1 + n2 * T + n3 * pow(T, 2.0)) * S + n4 * pow(T, 2.0) + ((n5 + n6 * S + n7 * T) / wavelength) + (n8 / pow(wavelength, 2.0)) + (n9 / pow(wavelength, 3.0));
}
float Glass(in float wl) {
    return sqrt(1.0+1.03961212/(1.0-0.00600069867/pow(wl,2.0))+0.231792344/(1.0-0.0200179144/pow(wl,2.0))+1.01046945/(1.0-103.560653/pow(wl,2.0)));
}