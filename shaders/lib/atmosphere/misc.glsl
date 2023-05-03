#if !defined LIB_ATMOSPHERE_MISC
#define LIB_ATMOSPHERE_MISC
    float RayleighDensity(in float h) {
        return exp(-h * inverseScaleHeights.x + scaledPlanetRadius.x);
    }

    float AerosolDensity(in float h) {
        return exp(-h * inverseScaleHeights.y + scaledPlanetRadius.y);
    }

    float OzoneDensity(in float altitudeKm) {
        float o1 = 25.0 *     exp(( 0.0 - altitudeKm) /   8.0) * 0.5;
        float o2 = 30.0 * pow(exp((18.0 - altitudeKm) /  80.0), altitudeKm - 18.0);
        float o3 = 75.0 * pow(exp((25.3 - altitudeKm) /  35.0), altitudeKm - 25.3);
        float o4 = 50.0 * pow(exp((30.0 - altitudeKm) / 150.0), altitudeKm - 30.0);
        return (o1 + o2 + o3 + o4) / 134.628;
    }

    vec3 CalculateAtmosphereDensity(in float h) {
        vec2 rm;
        rm.x = RayleighDensity(h);
        rm.y = AerosolDensity(h);

        float altitudeKm = (h - planetRadius) / kilometer;

        float ozone = OzoneDensity(altitudeKm);

        return vec3(rm, ozone);
    }

    float AverageBetaM() {
        const float junge = 4.0;

        float mie = 0.0;
        for(int i = 0; i < 441; ++i) {
            float wavelength = (float(i + 1) / 441.0) * 441.0 + 390.0;
            float c = (0.6544 * TURBIDITY - 0.6510) * 4e-18;
            float K = (0.773335 - 0.00386891 * wavelength) / (1.0 - 0.00546759 * wavelength);
            mie += (0.434 * c * pi * pow(tau / (wavelength * 1e-9), junge - 2.0) * K) / 441.0;
        }

        return mie;
    }

    #ifdef NISHITA_MIE
        float BetaM(in float wavelength) {
            float B = TURBIDITY < 1.1 ? 0.0009 : 0.166 * TURBIDITY - 0.164;
            return pow(wavelength, -1.0) * B;
        }
    #else
        float BetaM(in float wavelength) {
            const float junge = 4.0;

            float c = (0.6544 * TURBIDITY - 0.6510) * 4e-18;
            float K = (0.773335 - 0.00386891 * wavelength) / (1.0 - 0.00546759 * wavelength);
            return 0.434 * c * pi * pow(tau / (wavelength * 1e-9), junge - 2.0) * K;
        }
    #endif

    float BetaR(in float wavelength) {
        float nanometers = wavelength * 1e-9;
        float micrometers = wavelength * 1e-6;

        float F_N2 = 1.034 + 3.17e-4 * (1.0 / pow(wavelength, 2.0));
        float F_O2 = 1.096 + 1.385e-3 * (1.0 / pow(wavelength, 2.0)) + 1.448e-4 * (1.0 / pow(wavelength, 4.0));
        float CCO2 = 0.045;
        float kingFactor = (78.084 * F_N2 + 20.946 * F_O2 + 0.934 + CCO2 * 1.15) / (78.084 + 20.946 + 0.934 + CCO2);
        float n = square(Air(wavelength * 1e-3)) - 1.0;

        return ((8.0 * pow(pi, 3.0) * pow(n, 2.0)) / (3.0 * airNumberDensity * pow(nanometers, 4.0))) * kingFactor;
    }

    float BetaO(float wavelength) {
        if(wavelength < 390.0 || wavelength > 830.0) return 0.0;

        return 0.0001 * ozoneNumberDensity * ozoneCrossSection[int(wavelength - 390.0)];
    }

    vec4 GetNoise(sampler3D noiseSampler, vec3 position) {
        return texture(noiseSampler, fract(position));
    }

    float remap(float value, float oldLow, float oldHight, float newLow, float newHight) {
        return newLow + (value - oldLow) * (newHight - newLow) / (oldHight - oldLow);
    }

    float HeightAlter(float percentHeight, float localizedCoverage) {
        float returnValue =  saturate(remap(percentHeight, 0.0, 0.07, 0.0, 1.0));
        float stopHeight = saturate(localizedCoverage + 0.6); 
        returnValue *= saturate(remap(percentHeight, stopHeight * 0.2, stopHeight, 1.0, 0.0));

        returnValue = pow(returnValue, saturate(remap(percentHeight, 0.65, 0.95, 1.0, (1 - cloudAnvilAmount * globalCoverage))));

        return returnValue;
    }

    float DensityAlter(float percentHeight, float localizedCoverage) {
        float returnValue = percentHeight;
        returnValue *= saturate(remap(percentHeight, 0.0, 0.2, 0.0, 1.0));
        returnValue *= localizedCoverage * 2.0;

        return returnValue;
    }

    float CalculateCloudShape(in vec3 position) {
        if(length(position) - planetRadius > cloudsMaxAltitude || length(position) - planetRadius < cloudsAltitude) return 0.0;
        vec3 cloudPosition = position / (cloudsThickness*32.0);

        float coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.1).r;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.2).r * 0.8 + coverageNoise;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.4).r * 0.6 + coverageNoise;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.8).r * 0.4 + coverageNoise;

        float localizedCoverage = coverageNoise * 3.0 - 2.0;
        float coverage = globalCoverage * localizedCoverage;

        float height = (length(position) - planetRadius) - cloudsAltitude;
        float normalizeHeight = height * (1.0 / cloudsThickness);
        float heightAlteration = HeightAlter(normalizeHeight, coverage);
        float densityAlteration = DensityAlter(normalizeHeight, coverage);
            densityAlteration = sqrt(densityAlteration * 0.25);

        float shapeNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 1.0).r;
              shapeNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 2.0).r * 0.8 + shapeNoise;
              shapeNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 3.0).r * 0.6 + shapeNoise;

        float mainNoise  = shapeNoise * 0.5;
              mainNoise -= coverage;
              mainNoise += heightAlteration * 1.1;

        float detailNoise  = GetNoise(noise3D, cloudPosition * pi).g * 0.8;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi).g * 0.7;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi * pi).g * 0.6;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi * pi * pi).g * 0.5;

        float cloudNoise = mainNoise - 0.5 * detailNoise;

        float density = saturate(1.0 * cloudNoise);

        density *= clamp(normalizeHeight * 2.0, 0.0, 1.0);

        if(density < 0.01) {
            return 0.0;
        }

        return 1.0 - pow(1.0 - density, 8.0);
    }
#endif