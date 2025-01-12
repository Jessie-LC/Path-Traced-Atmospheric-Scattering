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

        float altitude = abs(rayPosition.y - planetPosition.y) - planetRadius;

        if (dot(vec3(0.0, rayPosition.y - planetPosition.y, 0.0), rayVector) < 0.0) {
            // Ray is (locally) moving downwards
            // Find minimum altitude along the ray
            float minAltitude = 0.0;

            localMajorant  = coefficients.rayleigh * RayleighDensityExp(planetRadius + minAltitude);
            localMajorant += coefficients.aerosol * AerosolDensityExp(planetRadius + minAltitude);

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            return tSampled;
        } else {
            // Ray is (locally) moving upwards (or horizontally)
            localMajorant  = coefficients.rayleigh * RayleighDensityExp(altitude + planetRadius);
            localMajorant += coefficients.aerosol * AerosolDensityExp(altitude + planetRadius);

            float tSampled = -log(RandNextF()) / (tScale * localMajorant);

            return tSampled;
        }
    }

    bool IntersectionDeltaTracking(in vec3 position, in vec3 direction, in AttenuationCoefficients coefficients, out float volumeDistance, inout int component, in float tMax) {
        direction = normalize(direction);
        float aid;
        int atmosplane = RayPlaneIntersection(
            position, direction,
            vec3(0.0, atmosphereRadius, 0.0), vec3(0.0, -1.0, 0.0),
            aid
        );

        float t1 = min(aid, tMax);

        vec3 startPosition = position;

        float maxAttenuationCoefficient;
        float t = SampleNullDistance(
            startPosition,
            direction,
            coefficients,
            maxAttenuationCoefficient
        );

        int steps = 0;
        while(t < t1) {
            position = startPosition + direction * t;

            if (position.y < atmosphereLowerLimit) {
                return false;
            }
            if(position.y > atmosphereRadius) {
                return false;
            }

            float stepAttenuationRayleigh = coefficients.rayleigh * RayleighDensityExp(position.y);
            float stepAttenuationAerosol = coefficients.aerosol * AerosolDensityExp(position.y);

            float[2] stepAttenuationCoeffcients = float[](
                stepAttenuationRayleigh,
                stepAttenuationAerosol
            );

            float rand = RandNextF();
            if (rand < ((stepAttenuationRayleigh + stepAttenuationAerosol) / maxAttenuationCoefficient)) {
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

            steps += 1;

            t += SampleNullDistance(
                position,
                direction,
                coefficients,
                maxAttenuationCoefficient
            );
        }

        return false;
    }

    float TransmittanceAnalytical(in vec3 position, in vec3 direction, in AttenuationCoefficients coefficients) {
        float cosViewZenith = dot(direction, vec3(0.0, 1.0, 0.0));
        if (cosViewZenith <= 0.0) {
            return 0.0;
        }

        position.y -= planetRadius;

        float scaleSlopeViewRayleigh = scaleHeights.x / cosViewZenith;
        float scaleSlopeViewAerosol  = scaleHeights.y / cosViewZenith;
        float opticalDepthRayleigh   = coefficients.rayleigh * scaleSlopeViewRayleigh * exp(-position.y / scaleHeights.x);
        float opticalDepthAerosol    = coefficients.aerosol  *  scaleSlopeViewAerosol * exp(-position.y / scaleHeights.y);
        return exp(-(opticalDepthRayleigh + opticalDepthAerosol));
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
            float aid;
            int atmosplane = RayPlaneIntersection(
                position, rayDirection,
                vec3(0.0, atmosphereRadius, 0.0), vec3(0.0, -1.0, 0.0),
                aid
            );
            float pid;
            int groundPlane = RayPlaneIntersection(
                position, rayDirection,
                vec3(0.0, atmosphereLowerLimit, 0.0), vec3(0.0, 1.0, 0.0),
                pid
            );

            if(RandNextF() > throughput) {
                throughput = 0.0;
                break;
            }
            throughput /= saturate(throughput);

            vec3 normal = normalize(vec3(0.0, 1.0, 0.0));

            vec3 oldRayDirection = normalize(rayDirection);

            float hitDistance;
            bool hitGround = groundPlane == 1;
            bool hitAtmosphere = false;
            if (hitGround) {
                hitAtmosphere = IntersectionDeltaTracking(position, rayDirection, coefficients, hitDistance, component, pid);
            } else {
                hitAtmosphere = IntersectionDeltaTracking(position, rayDirection, coefficients, hitDistance, component, aid);
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

                    float transmittance = TransmittanceAnalytical(position, sunDirection, coefficients);

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
            } else if(hitGround) {
                position = position + rayDirection * pid;

                {
                    float transmittance = TransmittanceAnalytical(position, sunDirection, coefficients);

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
                    estimate += float(!hitGround) * PhysicalSun(rayDirection, lightVector, wavelength, irradiance / ConeAngleToSolidAngle(sunAngularRadius)) * throughput;
                }
                break;
            }

            ++bounces;
        }

        if(isinf(estimate)) estimate = 1.0;
        if(isnan(estimate)) estimate = 0.0;

        return estimate;
    }
#endif