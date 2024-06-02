#if !defined LIB_ATMOSPHERE_RAYMARCHED
#define LIB_ATMOSPHERE_RAYMARCHED
    #define SCATTERING_STEPS 128 //[8 16 32 64 128 256 512]
    #define TRANSMITTANCE_STEPS 64 //[8 16 32 64 128 256 512]

    float RaymarchAtmosphereTransmittance(in vec3 rayVector, in vec3 position, in vec4 baseAttenuationCoefficients) {
        float rayLength = dot(position, rayVector);
              rayLength = sqrt(rayLength * rayLength + square(atmosphereRadius) - dot(position, position)) - rayLength;
        float stepSize  = rayLength / float(TRANSMITTANCE_STEPS);
        vec3  increment = rayVector * stepSize;
        position += increment * RandNextF();

        vec3 thickness = vec3(0.0);
        for(int i = 0; i < TRANSMITTANCE_STEPS; ++i, position += increment) {
            vec3 density = CalculateAtmosphereDensity(length(position));
            if(density.x > 1e35) break;
            if(density.y > 1e35) break;
            if(density.z > 1e35) break;
            thickness += density;
        }

        float opticalDepth = (baseAttenuationCoefficients.x * stepSize * thickness.x) + (baseAttenuationCoefficients.y * stepSize * thickness.y) + (baseAttenuationCoefficients.z * stepSize * thickness.z);

        float transmittance = exp(-opticalDepth);
        if(isnan(transmittance)) {
            transmittance = 0.0;
        }
        if(isinf(transmittance)) {
            transmittance = 1.0;
        }

        return transmittance;
    }

    float RaymarchAtmosphereScattering(in vec3 viewPosition, in vec3 viewVector, in vec3 lightVector, in vec4 baseAttenuationCoefficients, in float irradiance, in float wavelength) {
        vec2 upperSphere = RSI(viewPosition, viewVector, atmosphereRadius);
        vec2 lowerSphere = RSI(viewPosition, viewVector, atmosphereLowerLimit);

        bool lowerSphereIntersected = lowerSphere.y >= 0.0;
        bool upperSphereIntersected = upperSphere.y >= 0.0;

        vec2 sd = vec2((lowerSphereIntersected && lowerSphere.x < 0.0) ? lowerSphere.y : max(upperSphere.x, 0.0), (lowerSphereIntersected && lowerSphere.x > 0.0) ? lowerSphere.x : upperSphere.y);

        float stepSize = length(sd.y - sd.x) / float(SCATTERING_STEPS);
        vec3 increment = viewVector * stepSize;
        vec3 position = viewVector * sd.x + (increment * RandNextF() + viewPosition);

        vec3 scatteringCoefficients = vec3(baseAttenuationCoefficients.x, baseAttenuationCoefficients.y * aerosolScatteringAlbedo, 0.0);

        vec3 sunDirection = GenerateConeVector(lightVector, RandNext2F(), sunAngularRadius);

        #if PHASE_FUNCTION_RAYLEIGH == 0
            float phaseR = RayleighPhase(dot(viewVector, sunDirection), wavelength);
        #elif PHASE_FUNCTION_RAYLEIGH == 1
            float phaseR = 0.25 / pi;
        #endif
        #if PHASE_FUNCTION_AEROSOL == 0
            float phaseM = AerosolPhase(dot(viewVector, sunDirection), wavelength);
        #elif PHASE_FUNCTION_AEROSOL == 1
            float phaseM = HenyeyGreensteinPhase(dot(viewVector, sunDirection), aerosol_g);
        #elif PHASE_FUNCTION_AEROSOL == 2
            float phaseM = KleinNishinaPhase(dot(viewVector, sunDirection), aerosol_g);
        #elif PHASE_FUNCTION_AEROSOL == 3
            float phaseM = ApproximateMiePhase(dot(viewVector, sunDirection), meanAerosolParticleDiameter);
        #endif

        float scattering = 0.0;
        float transmittance = 1.0;
        for(int i = 0; i < SCATTERING_STEPS; ++i, position += increment) {
            vec3 density = CalculateAtmosphereDensity(length(position));
            if(density.x > 1e35) break;
            if(density.y > 1e35) break;
            if(density.z > 1e35) break;
            vec3 airMass = stepSize * density;
            if(any(isnan(airMass))) {
                airMass = vec3(0.0);
            }
            float stepOpticalDepth = (baseAttenuationCoefficients.x * airMass.x) + (baseAttenuationCoefficients.y * airMass.y) + (baseAttenuationCoefficients.z * airMass.z);

            float stepTransmittance = saturate(exp(-stepOpticalDepth));
            float visibleScattering = transmittance * saturate((stepTransmittance - 1.0) / -stepOpticalDepth);

            float shadow = float(RSI(position, lightVector, planetRadius).x < 0.0);
            scattering += (scatteringCoefficients.x * airMass.x * phaseR + scatteringCoefficients.y * airMass.y * phaseM) * visibleScattering * RaymarchAtmosphereTransmittance(lightVector, position, baseAttenuationCoefficients) * shadow;
            transmittance *= stepTransmittance;

            if(RSI(position, viewVector, atmosphereLowerLimit).x <= 0.0) {
                //Break out of the loop early if this ray intersects the lower sphere
                break;
            }
        }

        if(isnan(scattering)) {
            return irradiance;
        }

        return scattering * irradiance;
    }
#endif