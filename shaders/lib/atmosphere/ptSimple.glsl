#if !defined LIB_ATMOSPHERE_PATHTRACED
#define LIB_ATMOSPHERE_PATHTRACED
    bool IntersectionDeltaTracking(in vec3 position, in vec3 direction, in vec3 baseAttenuationCoefficients, out float volumeDistance, inout int component) {
        vec2 aid = RSI(position, direction, atmosphereRadius);
        if (aid.y < 0.0) { return false; }
        vec2 pid = RSI(position, direction, atmosphereLowerLimit);
        bool planetIntersected = pid.y >= 0.0;
        bool atmosphereIntersected = aid.y >= 0.0;

        float tMax = (planetIntersected && pid.x > 0.0) ? pid.x : aid.y;

        vec3 maxAttenuationCoefficients = baseAttenuationCoefficients * CalculateAtmosphereDensity(atmosphereLowerLimit);
        float maxAttenuationCoefficient = saturate(maxAttenuationCoefficients.x + maxAttenuationCoefficients.y + maxAttenuationCoefficients.z);

        vec3 previousPosition = position;

        float t = 0.0;
        while(t < tMax) {
            float stepSize = -log(RandNextF()) / maxAttenuationCoefficient;

            position = previousPosition + direction * stepSize;

            t += stepSize;

            if (length(position) < atmosphereLowerLimit) {
                return false;
                break;
            }
            if(length(position) > atmosphereRadius) break;

            vec3 density = CalculateAtmosphereDensity(length(position));
            vec3 stepAttenuationCoeffcients = baseAttenuationCoefficients * density;

            float rand = RandNextF();
            if (rand < ((stepAttenuationCoeffcients.x + stepAttenuationCoeffcients.y + stepAttenuationCoeffcients.z) / maxAttenuationCoefficient)) {
                float fraction = 0.0;
                component = 0;
                while (component < 2) {
                    fraction += stepAttenuationCoeffcients[component];
                    if (rand * maxAttenuationCoefficient < fraction) { break; }
                    ++component;
                }
                volumeDistance = t;
                return true;
            }

            previousPosition = position;
        }

        return false;
    }
    float TransmittanceRatioTracking(in vec3 position, in vec3 direction, in vec3 baseAttenuationCoefficients) {
        vec2 aid = RSI(position, direction, atmosphereRadius);
        if (aid.y < 0.0) { return 0.0; }
        vec2 pid = RSI(position, direction, atmosphereLowerLimit);
        bool planetIntersected = pid.y >= 0.0;
        if(planetIntersected) { return 0.0; } //This is to fix an issue with the transmittance below the horizon.
        bool atmosphereIntersected = aid.y >= 0.0;

        float tMax = (planetIntersected && pid.x > 0.0) ? pid.x : aid.y;

        vec3 maxAttenuationCoefficients = baseAttenuationCoefficients * CalculateAtmosphereDensity(atmosphereLowerLimit);
        float maxAttenuationCoefficient = saturate(maxAttenuationCoefficients.x + maxAttenuationCoefficients.y + maxAttenuationCoefficients.z);

        vec3 previousPosition = position;

        float transmittance = 1.0;
        float t = 0.0;
        while(t < tMax) {
            float stepSize = -log(RandNextF()) / maxAttenuationCoefficient;

            position = previousPosition + direction * stepSize;

            t += stepSize;

            if (length(position) < atmosphereLowerLimit) {
                return 0.0;
                break;
            }
            if(length(position) > atmosphereRadius) break;

            vec3 density = CalculateAtmosphereDensity(length(position));
            vec3 stepAttenuationCoeffcients = baseAttenuationCoefficients * density;

            transmittance *= 1.0 - ((stepAttenuationCoeffcients.x + stepAttenuationCoeffcients.y + stepAttenuationCoeffcients.z) / maxAttenuationCoefficient);

            previousPosition = position;
        }

        return transmittance;
    }

    float PathtraceAtmosphereScattering(in vec3 viewPosition, in vec3 viewVector, in vec3 lightVector, in vec4 baseAttenuationCoefficients, in float irradiance, in float wavelength) {
        int component;

        float scattering = 0.0;

        vec3 position = viewPosition;
        vec3 rayDirection = viewVector;
        vec3 sunDirection = GenerateConeVector(lightVector, RandNext2F(), sunAngularRadius);

        float throughput = 1.0;
        int bounces = 0;
        while(bounces < SCATTERING_EVENTS) {
            if(RandNextF() > throughput) {
                throughput = 0.0;
                break;
            }
            throughput /= saturate(throughput);

            vec3 oldRayDirection = normalize(rayDirection);

            float hitDistance;
            bool hitAtmosphere = IntersectionDeltaTracking(position, rayDirection, baseAttenuationCoefficients.xyz, hitDistance, component);

            if(hitAtmosphere) {
                position = position + rayDirection * hitDistance;

                bool doScatter = false;
                switch (component) {
                    case 0: {
                        doScatter = true;
                        break; 
                    }
                    case 1: {
                        doScatter = RandNextF() < aerosolScatteringAlbedo;
                        break; 
                    }
                    case 2: { 
                        doScatter = false; 
                        break; 
                    }
                    case 3: {
                        doScatter = RandNextF() < cloudScatteringAlbedo;
                        break;
                    }
                }

                if(doScatter) {
                    float phase = 1.0;
                    switch (component) {
                        case 0: {
                            phase = RayleighPhase(dot(rayDirection, sunDirection));
                            break;
                        }
                        case 1: {
                            phase = AerosolPhase(dot(rayDirection, sunDirection), wavelength);
                            break;
                        }
                        case 2: {
                            phase = 0.0;
                            break;
                        }
                        case 3: {
                            phase = CloudPhase(dot(rayDirection, sunDirection), wavelength);
                            break;
                        }
                    }

                    float transmittance = TransmittanceRatioTracking(position, sunDirection, baseAttenuationCoefficients.xyz);

                    scattering += throughput * phase * irradiance * transmittance;

                    switch (component) {
                        case 0: {
                            rayDirection = Rotate(SampleRayleighPhase(), vec3(0.0, 0.0, 1.0), rayDirection);
                            throughput *= RayleighPhase(dot(oldRayDirection, rayDirection)) /
                                          RayleighPhase(dot(oldRayDirection, rayDirection));
                            break; 
                        }
                        case 1: {
                            rayDirection = Rotate(SampleAerosolPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                            throughput *= AerosolPhase(dot(oldRayDirection, rayDirection), wavelength) / 
                                          AerosolPhase(dot(oldRayDirection, rayDirection), wavelength);
                            break; 
                        }
                        case 2: { 
                            rayDirection = viewVector; 
                            break; 
                        }
                        case 3: { 
                            rayDirection = Rotate(SampleCloudPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                            throughput *= CloudPhase(dot(oldRayDirection, rayDirection), wavelength) / 
                                          CloudPhase(dot(oldRayDirection, rayDirection), wavelength);
                            break; 
                        }
                    }
                } else {
                    break;
                }
            } else {
                break;
            }

            ++bounces;
        }

        if(isinf(scattering)) scattering = 1.0;
        if(isnan(scattering)) scattering = 0.0;

        return scattering;
    }
#endif