#if !defined LIB_ATMOSPHERE_PHASE
#define LIB_ATMOSPHERE_PHASE
    float BinarySearchPhase(sampler1D cdfTexture, float lowBound, float highBound, float toFind) {
        float midIndex = (lowBound + highBound) / 2.0;

        int lookupRes = textureSize(cdfTexture, 0);
        for (int x = 0; x < 12; ++x) {
            float coordinate = midIndex / pi;
            float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / float(lookupRes);

            float value = saturate(texture(cdfTexture, coordTweaked).r); //The CDF value should never go above 1, but clamp it anyway just in case

            if (value < toFind) {
                lowBound = midIndex;
            }
            else if (value > toFind) {
                highBound = midIndex;
            }
            else {
                return midIndex;
            }

            midIndex = (lowBound + highBound) / 2.0;
        }

        return midIndex;
    }
    float BinarySearchPhase(sampler2D cdfTexture, float wavelength, float lowBound, float highBound, float toFind) {
        float midIndex = (lowBound + highBound) / 2.0;

        int lookupRes = textureSize(cdfTexture, 0).x;
        for (int x = 0; x < 12; ++x) {
            float coordinate = midIndex / pi;
            float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / float(lookupRes);

            float value = saturate(texture(cdfTexture, vec2(coordTweaked, (wavelength - 390.0) / 441.0)).r); //The CDF value should never go above 1, but clamp it anyway just in case

            if (value < toFind) {
                lowBound = midIndex;
            }
            else if (value > toFind) {
                highBound = midIndex;
            }
            else {
                return midIndex;
            }

            midIndex = (lowBound + highBound) / 2.0;
        }

        return midIndex;
    }

    vec3 SampleCloudPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureCloud, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleAerosolPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureAerosol, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleLowAltitudeAerosolPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureAerosolLowAltitude, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleRainbowPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRainbow, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleRayleighPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRayleigh, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }

    float CloudPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureCloud, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float AerosolPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureAerosol, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float LowAltitudeAerosolPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureAerosolLowAltitude, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RainbowPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureRainbow, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RayleighPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureRayleigh, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }

    float HenyeyGreensteinPhase(float cosTheta, float g) {
        const float norm = 0.25/pi;

        float gg = g * g;
        return norm * ((1.0-gg) / pow(1.0+gg-2.0*g*cosTheta, 3.0/2.0));
    }
	vec3 SampleHenyeyGreensteinPhase(float P, float g) {
		float s = 2.0 * P - 1.0;
		float u = (0.5 / g) * (1.0 + g * g - pow((1.0 - g * g) / (1.0 + g * s), 2.0));

		return GenerateUnitVector(vec2(RandNextF(), u * 0.5 + 0.5));
	}

    float KleinNishinaPhase(float cosTheta, float g) {
        float e = 1.0;
        for (int i = 0; i < 8; ++i) {
            float gFromE = 1.0 / e - 2.0 / log(2.0 * e + 1.0) + 1.0;
            float deriv = 4.0 / ((2.0 * e + 1.0) * square(log(2.0 * e + 1.0))) - 1.0 / square(e);
            if (abs(deriv) < 0.00000001) break;
            e = e - (gFromE - g) / deriv;
        }

        return e / (2.0 * pi * (e * (1.0 - cosTheta) + 1.0) * log(2.0 * e + 1.0));
    }
    vec3 SampleKleinNishinaPhase(float g) {
        float e = 1.0;
        for (int i = 0; i < 8; ++i) {
            float gFromE = 1.0 / e - 2.0 / log(2.0 * e + 1.0) + 1.0;
            float deriv = 4.0 / ((2.0 * e + 1.0) * square(log(2.0 * e + 1.0))) - 1.0 / square(e);
            if (abs(deriv) < 0.00000001) break;
            e = e - (gFromE - g) / deriv;
        }

        float cosTheta = (-pow(2.0 * e + 1.0, 1.0 - RandNextF()) + e + 1.0) / e;
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
#endif