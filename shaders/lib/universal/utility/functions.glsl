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

float CubicHermite00(float x) { return (1.0 + (2.0 * x)) * square(1.0 - x); }
float CubicHermite01(float x) { return square(x) * (3.0 - (2.0 * x)); }
float CubicHermite10(float x) { return x * square(1.0 - x); }
float CubicHermite11(float x) { return square(x) * (x - 1.0); }
float CubicHermiteSpline(in float pk, in float pk1, in float mk, in float mk1, in float x, in float xk, in float xk1) {
    float t = (x - xk) / (xk1 - xk);
    float scale = xk1 - xk;
    return CubicHermite00(t) * pk + CubicHermite10(t) * scale * mk + CubicHermite01(t) * pk1 + CubicHermite11(t) * scale * mk1;
}

float CardinalSplineM(in float pk1, in float pkn1, in float xk1, in float xkn1) {
    float c = 0.5;
    return (1.0 - c) * ((pk1 - pkn1) / (xk1 - xkn1));
}

float AddUvMargin(float uv, int   resolution) { return uv * (1.0 - 1.0 / resolution) + (0.5 / resolution); }
vec2  AddUvMargin(vec2 uv,  ivec2 resolution) { return uv * (1.0 - 1.0 / resolution) + (0.5 / resolution); }
vec3  AddUvMargin(vec3 uv,  ivec3 resolution) { return uv * (1.0 - 1.0 / resolution) + (0.5 / resolution); }
vec4  AddUvMargin(vec4 uv,  ivec4 resolution) { return uv * (1.0 - 1.0 / resolution) + (0.5 / resolution); }
float RemoveUvMargin(float uv, int   resolution) { return (uv - 0.5 / resolution) / (1.0 - 1.0 / resolution); }
vec2  RemoveUvMargin(vec2 uv,  ivec2 resolution) { return (uv - 0.5 / resolution) / (1.0 - 1.0 / resolution); }
vec3  RemoveUvMargin(vec3 uv,  ivec3 resolution) { return (uv - 0.5 / resolution) / (1.0 - 1.0 / resolution); }
vec4  RemoveUvMargin(vec4 uv,  ivec4 resolution) { return (uv - 0.5 / resolution) / (1.0 - 1.0 / resolution); }

// The following code is licensed under the MIT license: https://gist.github.com/TheRealMJP/bc503b0b87b643d3505d41eab8b332ae
vec4 textureCatmullRom(in sampler2D tex, vec2 uv) {
    vec2 texSize = textureSize(tex, 0);
    // We're going to sample a a 4x4 grid of texels surrounding the target UV coordinate. We'll do this by rounding
    // down the sample location to get the exact center of our "starting" texel. The starting texel will be at
    // location [1, 1] in the grid, where [0, 0] is the top left corner.
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5) + 0.5;

    // Compute the fractional offset from our starting texel to our original sample location, which we'll
    // feed into the Catmull-Rom spline function to get our filter weights.
    vec2 f = samplePos - texPos1;

    // Compute the Catmull-Rom weights using the fractional offset that we calculated earlier.
    // These equations are pre-expanded based on our knowledge of where the texels will be located,
    // which lets us avoid having to evaluate a piece-wise function.
    vec2 w0 = f * ( -0.5 + f * (1.0 - 0.5*f));
    vec2 w1 = 1.0 + f * f * (-2.5 + 1.5*f);
    vec2 w2 = f * ( 0.5 + f * (2.0 - 1.5*f) );
    vec2 w3 = f * f * (-0.5 + 0.5 * f);
    
    // Work out weighting factors and sampling offsets that will let us use bilinear filtering to
    // simultaneously evaluate the middle 2 samples from the 4x4 grid.
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);

    // Compute the final UV coordinates we'll use for sampling the texture
    vec2 texPos0 = texPos1 - vec2(1.0);
    vec2 texPos3 = texPos1 + vec2(2.0);
    vec2 texPos12 = texPos1 + offset12;

    texPos0 /= texSize;
    texPos3 /= texSize;
    texPos12 /= texSize;

    vec4 result = vec4(0.0);
    result += texture(tex, vec2(texPos0.x,  texPos0.y)) * w0.x * w0.y;
    result += texture(tex, vec2(texPos12.x, texPos0.y)) * w12.x * w0.y;
    result += texture(tex, vec2(texPos3.x,  texPos0.y)) * w3.x * w0.y;

    result += texture(tex, vec2(texPos0.x,  texPos12.y)) * w0.x * w12.y;
    result += texture(tex, vec2(texPos12.x, texPos12.y)) * w12.x * w12.y;
    result += texture(tex, vec2(texPos3.x,  texPos12.y)) * w3.x * w12.y;

    result += texture(tex, vec2(texPos0.x,  texPos3.y)) * w0.x * w3.y;
    result += texture(tex, vec2(texPos12.x, texPos3.y)) * w12.x * w3.y;
    result += texture(tex, vec2(texPos3.x,  texPos3.y)) * w3.x * w3.y;

    return result;
}
//End of MIT code

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

float RPI(in vec3 ro, in vec3 rd, in float y) {
    float t = -(ro.y - y) / rd.y;
    if (sign(t) == -1.0) t = -1.0;
    return t;
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
    // Returns radiance in units Watts per steradian per square meter surface
    const float h = 6.62607015e-16;
    const float c = 2.99792458e17;
    const float k = 1.380649e-5;

    float p1 = 2.0 * h * pow(c, 2.0) * pow(lambda, -5.0);
    float p2 = exp((h * c) / (lambda * k * t)) - 1.0;

    return p1 / p2;
}

float PlancksLawRadiance(float temperature) {
    // Returns radiance in units Watts per steradian per square meter surface
	const float h = 6.62607015e-16;
	const float c = 2.99792458e17;
	const float k_B = 1.380649e-5;
	const float constp = 2.0 * pow(3.1415926535 * k_B, 4.0) / (15.0 * c * c * h * h * h);
	float temperatureSquared = temperature * temperature;
	return constp * temperatureSquared * temperatureSquared;
}

float NormalDistribution(in float x, in float mean, in float standardDeviation) {
    return (1.0 / (standardDeviation * sqrt(2.0 * pi))) * exp(-0.5 * square((x - mean) / standardDeviation));
}

vec2 SampleCircle(in vec2 rng, in float radius) {
    float r = radius * sqrt(rng.x);
    float t = 2.0 * pi * rng.y;
    
    return r * vec2(cos(t), sin(t));
}

uniform sampler2D CIELUT;
uniform sampler2D CIELUT_1931;
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

    return xyz;
}
vec3 SpectrumToXYZExact_CIE1931(in float spectrum, in float w) {
    float n = (w - 360.0);
    int i = int(n);

    int n0 = int(n);
    int n1 = min(n0 + 1, 470);

    vec3 cie0 = texelFetch(CIELUT_1931, ivec2(n0, 1), 0).xyz;
    vec3 cie1 = texelFetch(CIELUT_1931, ivec2(n1, 1), 0).xyz;

    vec3 xyz = mix(cie0, cie1, fract(n));

    xyz = vec3(spectrum * 683.0) * xyz;

    return xyz;
}
vec3 SpectrumToXYZExact_CIE1964(in float spectrum, in float w) {
    float n = (w - 360.0);
    int i = int(n);

    int n0 = int(n);
    int n1 = min(n0 + 1, 470);

    vec3 cie0 = CIE_1964[n0];
    vec3 cie1 = CIE_1964[n1];

    vec3 xyz = mix(cie0, cie1, fract(n));

    xyz = vec3(spectrum * 683.0) * xyz;

    return xyz;
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

//From: https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py
float Z(in float T, in float p, in float xw) {
    const float a0 = 1.58123e-6;
    const float a1 = -2.9331e-8;
    const float a2 = 1.1043e-10;
    const float b0 = 5.707e-6;
    const float b1 = -2.051e-8;
    const float c0 = 1.9898e-4;
    const float c1 = -2.376e-6;
    const float d  = 1.83e-11;
    const float e  = -0.765e-8;
    return 1.0 - (p / T) * (a0 + a1 * T + a2 * pow(T, 2.0) + (b0 + b1 * T) * xw + (c0 + c1 * T) * pow(xw, 2.0)) + pow(p / T, 2.0) * (d + e * pow(xw, 2.0));
}
//From: https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py
float AirN(in float lambda, in float T, in float p, in float h) {
    const float xc = 280.0;
    const float R = 8.314510;
    const float k0 = 238.0185;
    const float k1 = 5792105.0;
    const float k2 = 57.362;
    const float k3 = 167917.0;
    const float w0 = 295.235;
    const float w1 = 2.6422;
    const float w2 = -0.032380;
    const float w3 = 0.004028;
    const float A = 1.2378847e-5;
    const float B = -1.9121316e-2;
    const float C = 33.93711047;
    const float D = 6.3431645e3;
    const float alpha = 1.00062;
    const float beta = 3.14e-8;
    const float gamma = 5.6e-7;

    float sigma = 1.0 / lambda;
    float sigma2 = pow(sigma, 2.0);

    float svp = 0.0;
    if(T - 273.15 >= 0.0) {
        svp = exp(A * pow(T, 2.0) + B * T + C + D / T);
    } else {
        svp = pow(10.0, -2663.5 / T + 12.537);
    }
    
    float f = alpha + beta * p + gamma * pow(T - 273.15, 2.0);
    
    float xw = f * h * svp / p;

    float nas = 1.0 + (k1 / (k0 - sigma2) + k3 / (k2 - sigma2)) * 1e-8;
    float naxs = 1.0 + (nas - 1.0) * (1.0 + 0.534e-6 * (xc - 450.0));
    float nws = 1.0 + 1.022 * (w0 + w1 * sigma2 + w2 * pow(sigma, 4.0) + w3 * pow(sigma, 6.0)) * 1e-8;

    float Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400.0));
    float Mw = 0.018015;

    float Za = Z(288.15, 101325.0, 0.0);
    float Zw = Z(293.15, 1333.0, 1.0);

    float paxs = 101325.0 * Ma / (Za * R * 288.15);
    float pws = 1333.0 * Mw / (Zw * R * 293.15);

    float pa = p * Ma / (Z(T, p, xw) * R * T) * (1.0 - xw);
    float pw = p * Mw / (Z(T, p, xw) * R * T) * xw;

    return 1.0 + (pa / paxs) * (naxs - 1.0) + (pw / pws) * (nws - 1.0);
}