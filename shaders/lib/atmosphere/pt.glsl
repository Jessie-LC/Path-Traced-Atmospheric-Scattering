#if !defined LIB_ATMOSPHERE_PATHTRACED
#define LIB_ATMOSPHERE_PATHTRACED
    struct BaseAttenuationCoefficients {
        float rayleigh;
        float aerosol;
        float ozone;
        float cloud;
        float mist;
    };

    float SampleNullDistance(
        vec3  rayPosition,
        vec3  rayVector,
        vec4 baseAttenuationCoefficients,
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

            localMajorant  = baseAttenuationCoefficients.x * CalculateAtmosphereDensity(planetRadius + minAltitude).x;
            localMajorant += baseAttenuationCoefficients.y * CalculateAtmosphereDensity(planetRadius + minAltitude).y;

            if (altitude <= cloudsAltitude) {
                // Below cloud layer
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(planetRadius + minAltitude).z;
            } else if (altitude < cloudsMaxAltitude) {
                // Inside cloud layer
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(planetRadius + minAltitude).z;
                localMajorant += baseAttenuationCoefficients.w;
            } else if (altitude < 25e3) {
                // Below ozone peak
                localMajorant += max(
                    baseAttenuationCoefficients.z * CalculateAtmosphereDensity(planetRadius + minAltitude).z,
                    baseAttenuationCoefficients.z * CalculateAtmosphereDensity(altitude + planetRadius).z
                );
            } else {
                // Above ozone peak
                localMajorant += max(
                    baseAttenuationCoefficients.z * CalculateAtmosphereDensity(planetRadius + minAltitude).z,
                    baseAttenuationCoefficients.z * CalculateAtmosphereDensity(atmosphereLowerLimit + 25e3).z
                );
            }

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            if (altitude >= cloudsMaxAltitude) {
                float t0, t1;
                if (LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsMaxAltitude, t0, t1) != 0 && t0 < tSampled) {
                    localMajorant += baseAttenuationCoefficients.w;
                    tSampled = t0 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            if (altitude <= cloudsAltitude) {
                float t0, t1;
                LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsAltitude, t0, t1);
                if (t1 < tSampled) {
                    localMajorant += baseAttenuationCoefficients.w;
                    tSampled = t1 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            return tSampled;
        } else {
            // Ray is (locally) moving upwards (or horizontally)
            localMajorant  = baseAttenuationCoefficients.x * CalculateAtmosphereDensity(altitude + planetRadius).x;
            localMajorant += baseAttenuationCoefficients.y * CalculateAtmosphereDensity(altitude + planetRadius).y;

            if (altitude <= cloudsAltitude) {
                // Below cloud layer
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(altitude + planetRadius).z;
            } else if (altitude < cloudsMaxAltitude) {
                // Inside cloud layer
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(altitude + planetRadius).z;
                localMajorant += baseAttenuationCoefficients.w;
            } else if (altitude < 25e3) {
                // Below ozone peak
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(atmosphereLowerLimit+25e3).z;
            } else {
                // Above ozone peak
                localMajorant += baseAttenuationCoefficients.z * CalculateAtmosphereDensity(altitude + planetRadius).z;
            }

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            if (altitude <= cloudsAltitude) {
                float t0, t1;
                LineSphereIntersect(rayPosition, rayVector, planetPosition, planetRadius + cloudsAltitude, t0, t1);
                if (t1 < tSampled) {
                    localMajorant += baseAttenuationCoefficients.w;
                    tSampled = t1 + -log(RandNextF()) / (tScale * localMajorant);
                }
            }

            return tSampled;
        }
    }

    bool FindNextInteraction(in vec3 position,in  vec3 rayDirection, in vec4 baseAttenuationCoefficients, out int component, out float volumeDistance, in float tMax) {
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
            baseAttenuationCoefficients,
            maxAttenuationCoefficient
        );

        while (t < t1) {
            position = startPosition + rayDirection * t;

            if (length(position) < atmosphereLowerLimit) {
                return false;
                break;
            }

            float baseDensity = CalculateAtmosphereDensity(length(position)).x;
            vec3 volumeGradient = normalize(vec3(
                CalculateAtmosphereDensity(length(position + vec3(0.01, 0.0, 0.0))).x - baseDensity,
                CalculateAtmosphereDensity(length(position + vec3(0.0, 0.01, 0.0))).x - baseDensity,
                CalculateAtmosphereDensity(length(position + vec3(0.0, 0.0, 0.01))).x - baseDensity
            ));

            float densityRayleigh = CalculateAtmosphereDensity(length(position)).x;
            float densityMie      = CalculateAtmosphereDensity(length(position)).y;
            float densityOzo      = CalculateAtmosphereDensity(length(position)).z;
            if (densityOzo > 1e35) { break; }
            if (densityMie > 1e35) { break; }
            if (densityRayleigh > 1e35) { break; }

            float rayleighAttenuation = baseAttenuationCoefficients.x * densityRayleigh;
            float mieAttenuation      = baseAttenuationCoefficients.y * densityMie;
            float ozoneAttenuation    = baseAttenuationCoefficients.z * densityOzo;
            float cloudAttenuation    = 0.0;
            if(length(position) > cloudsAltitude || length(position) < cloudsMaxAltitude) {
                cloudAttenuation = CalculateCloudShape(position) * baseAttenuationCoefficients.w;
            }
            vec4 stepAttenuationCoeffcients = vec4(rayleighAttenuation, mieAttenuation, ozoneAttenuation, cloudAttenuation);

            float rand = RandNextF();
            if (rand < ((rayleighAttenuation + mieAttenuation + ozoneAttenuation + cloudAttenuation) / maxAttenuationCoefficient)) {
                float fraction = 0.0;
                component = 0;
                while (component < 3) {
                    fraction += stepAttenuationCoeffcients[component];
                    if (rand * maxAttenuationCoefficient < fraction) { break; }
                    ++component;
                }
                volumeDistance = t;
                return true;
            }

            t += SampleNullDistance(
                position,
                rayDirection,
                baseAttenuationCoefficients,
                maxAttenuationCoefficient
            );
        }

        return false;
    }

    float EstimateTransmittance(in vec3 position, in vec3 rayDirection, in vec4 baseAttenuationCoefficients) {
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
        while (t < tMax) {
            if (RandNextF() > transmittance) {
                transmittance = 0.0;
                break;
            }
            transmittance /= transmittance;

            float stepSize = SampleNullDistance(
                position,
                rayDirection,
                baseAttenuationCoefficients,
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
            if (densityOzo > 1e35) { break; }
            if (densityMie > 1e35) { break; }
            if (densityRayleigh > 1e35) { break; }

            float rayleighAttenuation = baseAttenuationCoefficients.x * densityRayleigh;
            float mieAttenuation      = baseAttenuationCoefficients.y * densityMie;
            float ozoneAttenuation    = baseAttenuationCoefficients.z * densityOzo;
            float cloudAttenuation    = 0.0;
            if(length(position) > cloudsAltitude || length(position) < cloudsMaxAltitude) {
                cloudAttenuation = CalculateCloudShape(position) * baseAttenuationCoefficients.w;
            }

            transmittance *= 1.0 - ((rayleighAttenuation + mieAttenuation + ozoneAttenuation + cloudAttenuation) / maxAttenuationCoefficient);

            previousPosition = position;
        }

        return saturate(transmittance);
    }

    float PathtraceAtmosphereScattering(in vec3 viewPosition, in vec3 viewVector, in vec3 lightVector, in vec4 baseAttenuationCoefficients, in float irradiance, in float wavelength) {
        int component;

        float scattering = 0.0;

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
                hitAtmosphere = FindNextInteraction(position, rayDirection, baseAttenuationCoefficients, component, hitDistance, pid.x);
            } else {
                hitAtmosphere = FindNextInteraction(position, rayDirection, baseAttenuationCoefficients, component, hitDistance, aid.y);
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

                    float transmittance = EstimateTransmittance(position, sunDirection, baseAttenuationCoefficients);

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
            } else if(planetIntersected) {
                position = position + rayDirection * pid.x;

                {
                    float transmittance = EstimateTransmittance(position, sunDirection, baseAttenuationCoefficients);

                    float bsdf = groundAlbedo * saturate(dot(sunDirection, normal)) / pi;

                    scattering += throughput * bsdf * irradiance * transmittance;
                }

                {
                    rayDirection = GenerateCosineVector(normal, RandNext2F());
                    throughput *= groundAlbedo;

                    if(dot(rayDirection, normal) < 0.0) {
                        break;
                    }
                }
            } else {
                if(bounces < 1) scattering += PhysicalSun(rayDirection, lightVector, wavelength, Plancks(5778.0, wavelength)) * throughput;
                break;
            }

            ++bounces;
        }

        if(isinf(scattering)) scattering = 1.0;
        if(isnan(scattering)) scattering = 0.0;

        return scattering;
    }
#endif