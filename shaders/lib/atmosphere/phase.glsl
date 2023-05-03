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
    vec3 SampleRainbowPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRainbow, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleRayleighPhase() {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRayleigh, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }

    float CloudPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureCloud, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float AerosolPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureAerosol, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RainbowPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureRainbow, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RayleighPhase(in float cosTheta) {
        return texture(phaseTextureRayleigh, acos(cosTheta) / pi).r;
    }
#endif