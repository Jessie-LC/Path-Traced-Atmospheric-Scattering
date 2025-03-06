#if !defined LIB_ATMOSPHERE_PATHTRACED
#define LIB_ATMOSPHERE_PATHTRACED
    float SampleNullDistance(
        vec3  rayPosition,
        vec3  rayVector,
        AttenuationCoefficients coefficients,
        out float localMajorant
    ) {
        const vec3 planetPosition = vec3(0.0);

        float tScale = length(rayVector);

        float altitude = distance(rayPosition, planetPosition) - planetRadius;

        if (dot(rayPosition - planetPosition, rayVector) < 0.0) {
            // Ray is (locally) moving downwards
            // Find minimum altitude along the ray

            float tMinAltitude = dot(planetPosition - rayPosition, rayVector) / dot(rayVector, rayVector);
            vec3  pMinAltitude = rayPosition + rayVector * tMinAltitude;
            float minAltitude = distance(planetPosition, pMinAltitude) - planetRadius;

            minAltitude = max(minAltitude, 0.0);
            pMinAltitude = planetPosition + vec3(0.0, planetRadius + minAltitude, 0.0);

            localMajorant  = coefficients.rayleigh * CalculateAtmosphereDensity(planetRadius + minAltitude).x;
            localMajorant += coefficients.aerosol * CalculateAtmosphereDensity(planetRadius + minAltitude).y;
            localMajorant += coefficients.mist * CalculateMistDensity(planetRadius + minAltitude);

            if (altitude <= cloudsAltitude) {
                // Below cloud layer
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(planetRadius + minAltitude).z;
            } else if (altitude < cloudsMaxAltitude) {
                // Inside cloud layer
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(planetRadius + minAltitude).z;
                localMajorant += coefficients.cloud;
            } else if (altitude < 25e3) {
                // Below ozone peak
                localMajorant += max(
                    coefficients.ozone * CalculateAtmosphereDensity(planetRadius + minAltitude).z,
                    coefficients.ozone * CalculateAtmosphereDensity(altitude + planetRadius).z
                );
            } else {
                // Above ozone peak
                localMajorant += max(
                    coefficients.ozone * CalculateAtmosphereDensity(planetRadius + minAltitude).z,
                    coefficients.ozone * CalculateAtmosphereDensity(atmosphereLowerLimit + 25e3).z
                );
            }

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            if (altitude >= cloudsMaxAltitude) {
                float t0, t1;
                if (LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsMaxAltitude, t0, t1) != 0 && t0 < tSampled) {
                    localMajorant += coefficients.cloud;
                    tSampled = t0 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            if (altitude <= cloudsAltitude) {
                float t0, t1;
                LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsAltitude, t0, t1);
                if (t1 < tSampled) {
                    localMajorant += coefficients.cloud;
                    tSampled = t1 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            return tSampled;
        } else {
            // Ray is (locally) moving upwards (or horizontally)
            localMajorant  = coefficients.rayleigh * CalculateAtmosphereDensity(altitude + planetRadius).x;
            localMajorant += coefficients.aerosol * CalculateAtmosphereDensity(altitude + planetRadius).y;
            localMajorant += coefficients.mist * CalculateMistDensity(altitude + planetRadius);

            if (altitude <= cloudsAltitude) {
                // Below cloud layer
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(altitude + planetRadius).z;
            } else if (altitude < cloudsMaxAltitude) {
                // Inside cloud layer
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(altitude + planetRadius).z;
                localMajorant += coefficients.cloud;
            } else if (altitude < 25e3) {
                // Below ozone peak
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(atmosphereLowerLimit+25e3).z;
            } else {
                // Above ozone peak
                localMajorant += coefficients.ozone * CalculateAtmosphereDensity(altitude + planetRadius).z;
            }

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            if (altitude <= cloudsAltitude) {
                float t0, t1;
                LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsAltitude, t0, t1);
                if (t1 < tSampled) {
                    localMajorant += coefficients.cloud;
                    tSampled = t1 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            return tSampled;
        }
    }

    int FindNextInteraction(
        in vec3 position,
        in vec3 rayDirection, 
        in AttenuationCoefficients coefficients,
        out int component, 
        out float volumeDistance, 
        in float tMax
    ) {
        // Find the interaction point with the atmosphere and output the distance to said point alongside the component of the atmosphere that was interacted with
        rayDirection = normalize(rayDirection);
        vec2 aid = RSI(position, rayDirection, atmosphereRadius);

        float t = max(aid.x, 0.0);
        float t1 = min(aid.y, tMax);

        float maxAttenuationCoefficient;

        vec3 startPosition = position;

        t += SampleNullDistance(
            startPosition,
            rayDirection,
            coefficients,
            maxAttenuationCoefficient
        );

        while (t < t1) {
            position = startPosition + rayDirection * t;

            if (length(position) < atmosphereLowerLimit) {
                return -1;
                break;
            }
            if (length(position) > atmosphereRadius) {
                return -1;
                break;
            }

            float densityRayleigh = CalculateAtmosphereDensity(length(position)).x;
            float densityMie      = CalculateAtmosphereDensity(length(position)).y;
            float densityOzo      = CalculateAtmosphereDensity(length(position)).z;
            float densityMist     = CalculateMistDensity(length(position));
            if (densityOzo > 1e35) { break; }
            if (densityMie > 1e35) { break; }
            if (densityRayleigh > 1e35) { break; }

            float rayleighAttenuation = coefficients.rayleigh * densityRayleigh;
            float mieAttenuation      = coefficients.aerosol * densityMie;
            float ozoneAttenuation    = coefficients.ozone * densityOzo;
            float cloudAttenuation    = 0.0;
            float mistAttenuation     = coefficients.mist * densityMist;
            if(length(position) > cloudsAltitude || length(position) < cloudsMaxAltitude) {
                cloudAttenuation = coefficients.cloud * CalculateCloudShape(position);
            }
            float[5] stepAttenuationCoeffcients = float[](
                rayleighAttenuation,
                mieAttenuation, 
                ozoneAttenuation, 
                cloudAttenuation,
                mistAttenuation
            );

            float rand = RandNextF();
            if (rand < ((rayleighAttenuation + mieAttenuation + ozoneAttenuation + cloudAttenuation + mistAttenuation) / maxAttenuationCoefficient)) {
                float fraction = 0.0;
                component = 0;
                while (component < 5) {
                    fraction += stepAttenuationCoeffcients[component];
                    if (rand * maxAttenuationCoefficient < fraction) { break; }
                    ++component;
                }
                volumeDistance = t;
                return 1;
            }

            t += SampleNullDistance(
                position,
                rayDirection,
                coefficients,
                maxAttenuationCoefficient
            );
        }

        return 0;
    }

    float EstimateTransmittance(in vec3 position, in vec3 rayDirection, in AttenuationCoefficients coefficients) {
        // Estimate the atmosphere transmittance
        vec2 aid = RSI(position, rayDirection, atmosphereRadius);
        if (aid.y < 0.0) { return 0.0; }
        vec2 pid = RSI(position, rayDirection, atmosphereLowerLimit);
        bool planetIntersected = pid.x >= 0.0;
        if(planetIntersected) { return 0.0; } //This is to fix an issue with the transmittance below the horizon.
        bool atmosphereIntersected = aid.y >= 0.0;

        float tMax = (planetIntersected && pid.x > 0.0) ? pid.x : aid.y;

        float maxAttenuationCoefficient;

        vec3 previousPosition = position;

        float transmittance = 1.0;
        float t = 0.0;
        int i = 0;
        while (t < tMax) {
            ++i;
            if (i > 1000) {
                break;
            }

            if (RandNextF() > transmittance) {
                transmittance = 0.0;
                break;
            }
            transmittance /= transmittance;

            float stepSize = SampleNullDistance(
                position,
                rayDirection,
                coefficients,
                maxAttenuationCoefficient
            );

            t += stepSize;

            position = previousPosition + rayDirection * stepSize;

            if (length(position) < atmosphereLowerLimit) {
                transmittance *= 0.0;
                return transmittance;
                break; 
            }
            if (length(position) > atmosphereRadius) { break; }

            float densityRayleigh = CalculateAtmosphereDensity(length(position)).x;
            float densityMie      = CalculateAtmosphereDensity(length(position)).y;
            float densityOzo      = CalculateAtmosphereDensity(length(position)).z;
            float densityMist     = CalculateMistDensity(length(position));
            if (densityOzo > 1e35) { break; }
            if (densityMie > 1e35) { break; }
            if (densityRayleigh > 1e35) { break; }

            float rayleighAttenuation = coefficients.rayleigh * densityRayleigh;
            float mieAttenuation      = coefficients.aerosol * densityMie;
            float ozoneAttenuation    = coefficients.ozone * densityOzo;
            float cloudAttenuation    = 0.0;
            float mistAttenuation     = coefficients.mist * densityMist;
            if(length(position) > cloudsAltitude || length(position) < cloudsMaxAltitude) {
                cloudAttenuation = coefficients.cloud * CalculateCloudShape(position);
            }

            transmittance *= 1.0 - ((rayleighAttenuation + mieAttenuation + ozoneAttenuation + cloudAttenuation + mistAttenuation) / maxAttenuationCoefficient);

            previousPosition = position;
        }

        return saturate(transmittance);
    }

    float PathtraceAtmosphereScattering(
        in vec3 viewPosition, 
        in vec3 viewVector, 
        in vec3 lightVector, 
        in AttenuationCoefficients coefficients, 
        in float irradiance, 
        in float wavelength
    ) {
        const float ballRadius = 1e6;
        const vec3 ballPosition = vec3(1e3, planetRadius + ballRadius, 1e3);
        int component;

        float estimate = 0.0;

        vec3 groundAlbedoRGB = SrgbToLinear(vec3(GROUND_ALBEDO_R,GROUND_ALBEDO_G,GROUND_ALBEDO_B) / 255.0);

        float groundAlbedo;
        RGBToSpectrum(groundAlbedo, wavelength, groundAlbedoRGB.r, groundAlbedoRGB.g, groundAlbedoRGB.b, 0);
        groundAlbedo = saturate(groundAlbedo);

        vec3 position = viewPosition;
        vec3 rayDirection = viewVector;
        vec3 sunDirection = GenerateConeVector(lightVector, RandNext2F(), sunAngularRadius);

        float throughput = 1.0;
        int bounces = 0;
        while(bounces < SCATTERING_EVENTS) {
            vec2 aid = RSI(position, rayDirection, atmosphereRadius);
            vec2 pid = RSI(position, rayDirection, atmosphereLowerLimit);
            bool planetIntersected = pid.x > 0.0;

            if(RandNextF() > throughput) {
                throughput = 0.0;
                break;
            }
            throughput /= saturate(throughput);

            vec3 normal = normalize(position + rayDirection * pid.x);

            vec3 oldRayDirection = normalize(rayDirection);

            float hitDistance;
            bool hitAtmosphere;
            if(planetIntersected) {
                int interaction = FindNextInteraction(
                    position, 
                    rayDirection, 
                    coefficients, 
                    component, 
                    hitDistance, 
                    pid.x
                );

                hitAtmosphere = interaction == 1;
            } else {
                int interaction = FindNextInteraction(
                    position, 
                    rayDirection, 
                    coefficients, 
                    component, 
                    hitDistance, 
                    aid.y
                );

                hitAtmosphere = interaction == 1;
            }

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
                    case 4: {
                        doScatter = RandNextF() < cleanAerosolAlbedo;
                        break;
                    }
                }

                if(doScatter) {
                    float phase = 1.0;
                    switch (component) {
                        case 0: {
                            #if PHASE_FUNCTION_RAYLEIGH == 0
                                phase = RayleighPhase(dot(rayDirection, sunDirection), wavelength);
                            #elif PHASE_FUNCTION_RAYLEIGH == 1
                                phase = 0.25 / pi;
                            #endif
                            break;
                        }
                        case 1: {
                            #if PHASE_FUNCTION_AEROSOL == 0
                                phase = AerosolPhase(dot(rayDirection, sunDirection), wavelength);
                            #elif PHASE_FUNCTION_AEROSOL == 1
                                phase = HenyeyGreensteinPhase(dot(rayDirection, sunDirection), aerosol_g);
                            #elif PHASE_FUNCTION_AEROSOL == 2
                                phase = KleinNishinaPhase(dot(rayDirection, sunDirection), aerosol_g);
                            #elif PHASE_FUNCTION_AEROSOL == 3
                                phase = ApproximateMiePhase(dot(rayDirection, sunDirection), meanAerosolParticleDiameter);
                            #elif PHASE_FUNCTION_AEROSOL == 4
                                phase = RainbowPhase(dot(rayDirection, sunDirection), wavelength);
                            #endif
                            break;
                        }
                        case 2: {
                            phase = 0.0;
                            break;
                        }
                        case 3: {
                            #if PHASE_FUNCTION_CLOUD == 0
                                phase = CloudPhase(dot(rayDirection, sunDirection), wavelength);
                            #elif PHASE_FUNCTION_CLOUD == 1
                                phase = ApproximateMiePhase(dot(rayDirection, sunDirection), meanCloudParticleDiameter);
                            #endif
                            break;
                        }
                        case 4: {
                            phase = RainbowPhase(dot(rayDirection, sunDirection), wavelength);
                            break;
                        }
                    }

                    float transmittance = EstimateTransmittance(position, sunDirection, coefficients);

                    vec2 bid = RSI(position - ballPosition, sunDirection, ballRadius);
                    bool hitBall = bid.x >= 0.0;

                    estimate += throughput * phase * irradiance * transmittance;

                    switch (component) {
                        case 0: {
                            #if PHASE_FUNCTION_RAYLEIGH == 0
                                rayDirection = Rotate(SampleRayleighPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= RayleighPhase(dot(oldRayDirection, rayDirection), wavelength) /
                                              RayleighPhase(dot(oldRayDirection, rayDirection), wavelength);
                            #elif PHASE_FUNCTION_RAYLEIGH == 1
                                rayDirection = Rotate(GenerateUnitVector(RandNext2F()), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= (0.25 / pi) /
                                              (0.25 / pi);
                            #endif
                            break; 
                        }
                        case 1: {
                            #if PHASE_FUNCTION_AEROSOL == 0
                                rayDirection = Rotate(SampleAerosolPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= AerosolPhase(dot(oldRayDirection, rayDirection), wavelength) /
                                              AerosolPhase(dot(oldRayDirection, rayDirection), wavelength);
                            #elif PHASE_FUNCTION_AEROSOL == 1
                                rayDirection = Rotate(SampleHenyeyGreensteinPhase(RandNextF(), aerosol_g), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= HenyeyGreensteinPhase(dot(oldRayDirection, rayDirection), aerosol_g) /
                                              HenyeyGreensteinPhase(dot(oldRayDirection, rayDirection), aerosol_g);
                            #elif PHASE_FUNCTION_AEROSOL == 2
                                rayDirection = Rotate(SampleKleinNishinaPhase(aerosol_g), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= KleinNishinaPhase(dot(oldRayDirection, rayDirection), aerosol_g) /
                                              KleinNishinaPhase(dot(oldRayDirection, rayDirection), aerosol_g);
                            #elif PHASE_FUNCTION_AEROSOL == 3
                                rayDirection = Rotate(SampleApproximateMiePhase(meanAerosolParticleDiameter), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= ApproximateMiePhase(dot(oldRayDirection, rayDirection), meanAerosolParticleDiameter) /
                                              ApproximateMiePhase(dot(oldRayDirection, rayDirection), meanAerosolParticleDiameter);
                            #elif PHASE_FUNCTION_AEROSOL == 4
                                rayDirection = Rotate(SampleRainbowPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= RainbowPhase(dot(oldRayDirection, rayDirection), wavelength) /
                                              RainbowPhase(dot(oldRayDirection, rayDirection), wavelength);
                            #endif
                            break; 
                        }
                        case 2: { 
                            rayDirection = viewVector; 
                            break; 
                        }
                        case 3: { 
                            #if PHASE_FUNCTION_CLOUD == 0
                                rayDirection = Rotate(SampleCloudPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= CloudPhase(dot(oldRayDirection, rayDirection), wavelength) / 
                                              CloudPhase(dot(oldRayDirection, rayDirection), wavelength);
                            #elif PHASE_FUNCTION_CLOUD == 1
                                rayDirection = Rotate(SampleApproximateMiePhase(meanCloudParticleDiameter), vec3(0.0, 0.0, 1.0), rayDirection);
                                throughput *= ApproximateMiePhase(dot(oldRayDirection, rayDirection), meanCloudParticleDiameter) / 
                                              ApproximateMiePhase(dot(oldRayDirection, rayDirection), meanCloudParticleDiameter);
                            #endif
                            break; 
                        }
                        case 4: { 
                            rayDirection = Rotate(SampleRainbowPhase(wavelength), vec3(0.0, 0.0, 1.0), rayDirection);
                            throughput *= RainbowPhase(dot(oldRayDirection, rayDirection), wavelength) /
                                          RainbowPhase(dot(oldRayDirection, rayDirection), wavelength);
                            break; 
                        }
                    }
                } else {
                    break;
                }
            } else if(planetIntersected) {
                position = position + rayDirection * pid.x;

                {
                    float transmittance = EstimateTransmittance(position, sunDirection, coefficients);

                    float bsdf = groundAlbedo * saturate(dot(sunDirection, normal)) / pi;

                    estimate += throughput * bsdf * irradiance * transmittance;
                }

                {
                    rayDirection = GenerateCosineVector(normal, RandNext2F());
                    throughput *= groundAlbedo;

                    if(dot(rayDirection, normal) < 0.0) {
                        break;
                    }
                }
            } else {
                vec2 bid = RSI(position - ballPosition, rayDirection, ballRadius);
                bool hitBall = bid.x >= 0.0;
                if(bounces < 1) {
                    estimate += float(!planetIntersected) * PhysicalSun(rayDirection, lightVector, wavelength, irradiance / ConeAngleToSolidAngle(sunAngularRadius)) * throughput;
                }
                break;
            }

            ++bounces;
        }

        if(isinf(estimate)) estimate = 0.0;
        if(isnan(estimate)) estimate = 0.0;

        return estimate;
    }
#endif