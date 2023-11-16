#if !defined LIB_ATMOSPHERE_MISC
#define LIB_ATMOSPHERE_MISC
    float S(in float x) {
        return 10.0 * pow(x, 3.0) - 15.0 * pow(x, 4.0) + 6.0 * pow(x, 5.0);
    }
    float E(in float h) {
        return exp(-h / 40000.0) * S(min(max((h - (2000.0 + 500.0 / 2.0)) / ((2000.0 - 500.0 / 2.0) - (2000.0 + 500.0 / 2.0)), 0.0), 1.0));
    }
    float AerosolDensity(in float h) {
        return ((h < 11000.0 ? exp(-h / 80000.0) : exp(-11000.0 / 80000.0) * exp(-(h - 11000.0) / 2000.0)) + E(h) * ((0.166 * TURBIDITY - 0.164) * 1000.0)) / 2.0;
    }

    float RayleighDensityExp(in float h) {
        return exp(-h * inverseScaleHeights.x + scaledPlanetRadius.x);
    }
    float AerosolDensityExp(in float h) {
        return exp(-h * inverseScaleHeights.y + scaledPlanetRadius.y);
    }

    float OzoneDensity(in float altitudeKm) {
        float o1 = 25.0 *     exp(( 0.0 - altitudeKm) /   8.0) * 0.5;
        float o2 = 30.0 * pow(exp((18.0 - altitudeKm) /  80.0), altitudeKm - 18.0);
        float o3 = 75.0 * pow(exp((25.3 - altitudeKm) /  35.0), altitudeKm - 25.3);
        float o4 = 50.0 * pow(exp((30.0 - altitudeKm) / 150.0), altitudeKm - 30.0);
        return (o1 + o2 + o3 + o4) / 134.628;
    }

    vec2 US_StandardAtmosphereLookupUV(float R) {
        const float H = sqrt(atmosphereRadiusSquared - atmosphereLowerLimitSquared);

        float rho = sqrt(max(R * R - atmosphereLowerLimitSquared, 0.0));

        float uvR = AddUvMargin(rho / H, 512);

        return vec2(0.0, uvR);
    }
    vec3 GetUS_StandardAtmosphereLUT(in float coreDistance) {
        return texture(usStandardAtmosphere, US_StandardAtmosphereLookupUV(coreDistance)).rgb;
    }

    vec3 CalculateAtmosphereDensity(in float coreDistance) {
        float altitudeKm = (coreDistance - planetRadius) / kilometer;

        vec2 rm;
        #ifndef EXPONENTIAL_DENSITY
            rm.x = GetUS_StandardAtmosphereLUT(coreDistance).r;
            rm.y = AerosolDensity(coreDistance - planetRadius);
        #else
            rm.x = RayleighDensityExp(coreDistance);
            rm.y = AerosolDensityExp(coreDistance);
        #endif

        float ozone = OzoneDensity(altitudeKm);

        return max(vec3(rm, ozone), 0.0);
    }

    #ifndef EXPONENTIAL_DENSITY
        float BetaM(in float wavelength) {
            float B = 0.0009;
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

    float PreethamBetaO_Fit(in float wavelength) {
        wavelength = wavelength - 390.0;
        float p1 = NormalDistribution(wavelength, 202.0, 15.0) * 14.4;
        float p2 = NormalDistribution(wavelength, 170.0, 10.0) * 6.5;
        float p3 = NormalDistribution(wavelength, 50.0, 20.0) * 3.0;
        float p4 = NormalDistribution(wavelength, 100.0, 25.0) * 7.0;
        float p5 = NormalDistribution(wavelength, 140.0, 30.0) * 20.0;
        float p6 = NormalDistribution(wavelength, 150.0, 10.0) * 3.0;
        float p7 = NormalDistribution(wavelength, 290.0, 30.0) * 12.0;
        float p8 = NormalDistribution(wavelength, 330.0, 80.0) * 10.0;
        float p9 = NormalDistribution(wavelength, 240.0, 20.0) * 13.0;
        float p10 = NormalDistribution(wavelength, 220.0, 10.0) * 2.0;
        float p11 = NormalDistribution(wavelength, 186.0, 8.0) * 1.3;
        return 0.0001 * ozoneNumberDensity * ((p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11) / 1e20);
    }

    int BinarySearch(int lowIndex, int highIndex, float toFind) {
        while (lowIndex < highIndex) {
            int midIndex = (lowIndex + highIndex) >> 1;
            if (preethamWavelengths[midIndex] < toFind) {
                lowIndex = midIndex + 1;
            } else if (preethamWavelengths[midIndex] > toFind) {
                highIndex = midIndex;
            } else {
                return midIndex;
            }
        }
        return highIndex - 1;
    }

    float PreethamBetaO_Interpolated(in float wavelength) {
        if(wavelength < 455.0 || wavelength > 765.0) return 0.0;
        wavelength = clamp(wavelength, 450.0, 770.0);
        int k = BinarySearch(0, 32, wavelength);
        float pk = preethamOzone[k];
        float pk1 = preethamOzone[k + 1];
        float pk2 = preethamOzone[k + 2];
        float pkn1 = preethamOzone[k - 1];
        float xk = preethamWavelengths[k];
        float xk1 = preethamWavelengths[k + 1];
        float xk2 = preethamWavelengths[k + 2];
        float xkn1 = preethamWavelengths[k - 1];
        float mk = CardinalSplineM(pk1, pkn1, xk1, xkn1);
        float mk1 = CardinalSplineM(pk2, pk, xk2, xk);
        return 0.0001 * ozoneNumberDensity * (CubicHermiteSpline(pk, pk1, mk, mk1, wavelength, xk, xk1) / 2.5e21);
    }

    #ifdef USER_DEFINED_COEFFICIENTS
        float BetaM_Arbitrary(in float wavelength) {
            vec3 color = vec3(AEROSOL_COLOR_R, AEROSOL_COLOR_G, AEROSOL_COLOR_B) / 255.0;

            float spectrum;
            RGBToSpectrum(spectrum, wavelength, color.r, color.g, color.b, 0);
            return spectrum * 3e-6;
        }
        float BetaR_Arbitrary(in float wavelength) {
            vec3 color = vec3(RAYLEIGH_COLOR_R, RAYLEIGH_COLOR_G, RAYLEIGH_COLOR_B) / 255.0;

            float spectrum;
            RGBToSpectrum(spectrum, wavelength, color.r, color.g, color.b, 0);
            return spectrum * 3e-5;
        }
        float BetaO_Arbitrary(in float wavelength) {
            vec3 color = vec3(OZONE_COLOR_R, OZONE_COLOR_G, OZONE_COLOR_B) / 255.0;

            float spectrum;
            RGBToSpectrum(spectrum, wavelength, color.r, color.g, color.b, 0);
            return 5e-21 * spectrum * ozoneNumberDensity * 0.0001;
        }
    #endif

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

        position += vec3(1500.0, 0.0, 0.0);

        vec3 cloudPosition = position / (cloudsThickness*64.0);

        float coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.1).r;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.3).r * 0.8 + coverageNoise;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.4).r * 0.6 + coverageNoise;
              coverageNoise = GetNoise(noise3D, vec3(cloudPosition.x, 0.0, cloudPosition.z) * 0.5).r * 0.4 + coverageNoise;

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
              mainNoise += heightAlteration * 1.3;

        float detailNoise  = GetNoise(noise3D, cloudPosition * pi).g * 0.8;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi).g * 0.7;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi * pi).g * 0.6;
              detailNoise += GetNoise(noise3D, cloudPosition * pi * pi * pi * pi).g * 0.5;

        float cloudNoise = mainNoise - 0.4 * detailNoise;

        float density = saturate(3.0 * cloudNoise);

        density *= clamp(normalizeHeight * 2.0, 0.0, 5.0);

        return 1.0 - pow(1.0 - density, 2.0);
    }
#endif