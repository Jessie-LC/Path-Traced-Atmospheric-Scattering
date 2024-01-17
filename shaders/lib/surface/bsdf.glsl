#if !defined PLANET_BSDF
#define PLANET_BSDF
    #define MAX_BSDF_SCATTERING_EVENTS 20

    const float gamma = 3.0;

    float StudentT_D(in vec3 omegaH, in float gamma, in vec2 alpha) {
        //https://www.shadertoy.com/view/sdjXWw
        float cosTheta = omegaH.z;
        if(cosTheta <= 0.0) {
            return 0.0;
        }
        
        gamma = clamp(gamma, 1.501, 3.4e38); //This is to prevent the physically impossible case that occurs when gamma equals 1.5
        float cosThetaSquare = cosTheta * cosTheta;
        float exponent = ((square(omegaH.x) / square(alpha.x)) + (square(omegaH.y) / square(alpha.y))) / cosThetaSquare;
        float root = 1.0 + exponent / (gamma - 1.0);
        float powRootGamma = pow(root, gamma);
        
        return 1.0 / (pi * alpha.x * alpha.y * square(cosThetaSquare) * powRootGamma);
    }

    vec3 StudentT_SampleD(in vec2 U, in float gamma, in vec2 alpha) {
        /*
            From: https://mribar03.bitbucket.io/projects/eg_2017/
        */
        //The closer to infinity gamma goes, the closer to Beckmann this becomes.
        gamma = clamp(gamma, 1.501, 3.4e38); //This is to prevent the physically impossible case that occurs when gamma equals 1.5
        float phi = atan((alpha.x / alpha.y) * tan(2.0 * pi * U.x)) + pi * round(2.0 * U.x);
        float A = ((square(cos(phi)) / square(alpha.x)) + (square(sin(phi)) / square(alpha.y))) / (gamma - 1.0);
        float theta = atan(sqrt(
            (pow(1.0 - U.y, 1.0 / (1.0 - gamma)) - 1.0)
            / A
        ));

        vec3 wm = normalize(vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)));

        return wm;
    }

    float P22(in vec2 slope, in vec2 alpha) {
        if (any(equal(alpha, vec2(0.0)))) { return 0.0; }

        slope /= alpha;

        float tmp1 = slope.x*slope.x;
        float tmp2 = slope.y*slope.y;

        float value = (1.0 / (pi * alpha.x * alpha.y)) * exp(-tmp1 - tmp2);
        return value;
    }
    float D(in vec3 wm, in vec2 alpha) {
        if(wm.z <= 0.0) {
            return 0.0;
        }
        vec2 slope = -vec2(wm.x, wm.y) / wm.z;

        float value = P22(slope, alpha);
        return value / pow(normalize(wm).z, 4.0);
    }
    vec3 SampleD(in vec2 U, in vec2 alpha) {
        float phi = 2.0 * pi * U.y;
        float theta = atan(sqrt(-square(alpha.x) * log(1.0 - U.x)));

        vec3 wm = vec3(-cos(phi) * sin(theta), -sin(phi) * sin(theta), cos(theta));

        return wm;
    }

    float LambertSphere(in float cosTheta) {
        return (2.0 * (sqrt(1.0 - pow(cosTheta, 2.0)) - cosTheta * acos(cosTheta))) / (3.0 * pow(pi, 2.0));
    }
    vec3 SampleLambertSphere() {
        float p = RandNextF();
        float u = 1.0 - 2.0 * pow(1.0 - pow(p, 0.0401885 * p + 1.01938), 0.397225);

        return GenerateUnitVector(vec2(RandNextF(), u * 0.5 + 0.5));
    }

    vec3 SampleLambertSphereBSDF(
        in vec3 direction,
        in float albedo,
        in float roughness,
        inout float throughput
    ) {
        vec3 rayPosition = vec3(0.0);

        int scatteringOrder = 0;
        while(scatteringOrder <= MAX_BSDF_SCATTERING_EVENTS) {
            float tDist = -log(RandNextF());
            rayPosition = rayPosition + direction * tDist;
            if(rayPosition.z >= 0.0) {
                break;
            }

            scatteringOrder++;

            throughput *= albedo;

            direction = Rotate(SampleLambertSphere(), vec3(0.0, 0.0, 1.0), direction);

            if(direction.z != direction.z) {
                return vec3(0.0);
            }
        }

        return GenerateCosineVector(vec3(0.0, 0.0, 1.0), RandNext2F());
    }

    float EstimateLambertSphereBSDF(
        in vec3 incomingDirection, 
        in vec3 outgoingDirection, 
        in float albedo,
        in float roughness
    ) {
        if(outgoingDirection.z < 0.0) return 0.0;

        float estimate = 0.0;
        float throughput = 1.0;

        vec3 rayPosition = vec3(0.0);

        int scatteringOrder = 0;
        while(scatteringOrder <= MAX_BSDF_SCATTERING_EVENTS) {
            float tDist = -log(RandNextF());
            rayPosition = rayPosition + incomingDirection * tDist;
            if(rayPosition.z >= 0.0) {
                break;
            }

            scatteringOrder++;

            float outgoingDistance = -rayPosition.z / outgoingDirection.z;

            float phase = albedo * LambertSphere(dot(incomingDirection, outgoingDirection));
            estimate += phase * exp(-outgoingDistance) * throughput;        
            throughput *= albedo;

            incomingDirection = Rotate(SampleLambertSphere(), vec3(0.0, 0.0, 1.0), incomingDirection);

            if(incomingDirection.z != incomingDirection.z) {
                return 0.0;
            }
        }

        return albedo / pi;
    }

    const complexFloat n1 = complexFloat(1.00027, 0.0);
    const complexFloat n2 = complexFloat(1.33333, 0.0);

    float FacetProjectedWidth(in vec3 normal, in vec3 direction) {
        return (1.0 / normal.z) * dot(normal, direction);
    }
    float FacetProjectedVisibleWidthPrimary(in vec3 primaryNormal, in vec3 alternateNormal, in vec3 direction) {
        float   projectedWidthPrimary = FacetProjectedWidth(  primaryNormal, direction);
        float projectedWidthAlternate = FacetProjectedWidth(alternateNormal, direction);
        return max(projectedWidthPrimary + min(projectedWidthAlternate, 0.0), 0.0);
    }
    float FacetProjectedVisibleWidthAlternate(in vec3 primaryNormal, in vec3 alternateNormal, in vec3 direction) {
        float   projectedWidthPrimary = FacetProjectedWidth(  primaryNormal, direction);
        float projectedWidthAlternate = FacetProjectedWidth(alternateNormal, direction);
        return max(projectedWidthAlternate + min(projectedWidthPrimary, 0.0), 0.0);
    }
    float LambdaPrimary(in vec3 primaryNormal, in vec3 alternateNormal, in vec3 direction) {
        float   projectedVisibleWidthPrimary =   FacetProjectedVisibleWidthPrimary(primaryNormal, alternateNormal, direction);
        float projectedVisibleWidthAlternate = FacetProjectedVisibleWidthAlternate(primaryNormal, alternateNormal, direction);
        return projectedVisibleWidthPrimary / (projectedVisibleWidthPrimary + projectedVisibleWidthAlternate);
    }
    float LambdaAlternate(in vec3 primaryNormal, in vec3 alternateNormal, in vec3 direction) {
        float   projectedVisibleWidthPrimary =   FacetProjectedVisibleWidthPrimary(primaryNormal, alternateNormal, direction);
        float projectedVisibleWidthAlternate = FacetProjectedVisibleWidthAlternate(primaryNormal, alternateNormal, direction);
        return projectedVisibleWidthAlternate / (projectedVisibleWidthPrimary + projectedVisibleWidthAlternate);
    }

    float Height(in vec3 primaryNormal, in vec3 alternateNormal, in vec3 primaryOrigin, in vec3 alternateOrigin, in vec3 direction, in float position, in bool facet) {
        float HP = direction.z * (dot(  primaryNormal,   primaryOrigin - mix(primaryOrigin, alternateOrigin, position)) / dot(  primaryNormal, direction));
        float HA = direction.z * (dot(alternateNormal, alternateOrigin - mix(primaryOrigin, alternateOrigin, position)) / dot(alternateNormal, direction));
        return facet ? HP : HA;
    }

    vec3 SampleOpaqueDielectricFacetBSDF(in vec3 direction, in vec3 normal, in float albedo, inout float throughput) {
        complexFloat n1 = n1;
        complexFloat n2 = n2;

        float fresnel_R = FresnelNonPolarized_R(dot(normal, direction), n1, n2);
        float fresnel_T = FresnelNonPolarized_T(dot(normal, direction), n1, n2);

        float totalLum = albedo * fresnel_T + fresnel_R;
        float specBounceProbability = fresnel_R / totalLum;
        bool specularBounce = specBounceProbability > RandNextF();

        if(specularBounce) {
            throughput /= specBounceProbability;

            throughput *= fresnel_R;

            return reflect(-direction, normal);
        } else {
            throughput /= 1.0 - specBounceProbability;

            float energyConservationFactor = 1.0;
            uint N = 128u;
            for(uint i = 0u; i < N; ++i) {
                //Calculate the energy conservation factor via brute force integration.
                float nDotV = 1.0 - RandNext2F().x; //This is equivalent to nDotV here.
                float fresnel = FresnelNonPolarized_R(
                                    nDotV, 
                                    n1, 
                                    n2
                                );
                energyConservationFactor -= fresnel * (2.0 * nDotV) / float(N);
            }

            throughput *= fresnel_T;

            throughput /= energyConservationFactor;

            vec3 newDirection = GenerateCosineVector(normal, RandNext2F());

            throughput *= albedo;

            float fresnelL_T = FresnelNonPolarized_T(dot(normal, newDirection), n1, n2);

            throughput *= fresnelL_T;

            if(isnan(throughput)) throughput = 1.0;

            return newDirection;
        }
    }

    vec3 SampleOpaqueDielectric(in vec3 direction, in float albedo, in float roughness, inout float throughput) {
        vec3 geometricNormal = vec3(0.0, 0.0, 1.0);
        vec3 priFacetNormal = StudentT_SampleD(RandNext2F(), gamma, vec2(roughness));
        vec3 altFacetNormal = -reflect(priFacetNormal, geometricNormal);
        vec3 cavityWidthVector = normalize(geometricNormal * dot(geometricNormal, priFacetNormal) / dot(geometricNormal, geometricNormal) - priFacetNormal);
        vec3 priFacetOrigin =  cavityWidthVector;
        vec3 altFacetOrigin = -cavityWidthVector;

        float cavityPosition = RandNextF();
        bool isPriFacet = cavityPosition < LambdaPrimary(priFacetNormal, altFacetNormal, -direction);

        float height = Height(priFacetNormal, altFacetNormal, priFacetOrigin, altFacetOrigin, direction, cavityPosition, isPriFacet);

        vec3 rayEntry = mix(priFacetOrigin, altFacetOrigin, cavityPosition);
        vec3 rayPositionPriFacet = rayEntry + direction * (dot(priFacetNormal, priFacetOrigin - rayEntry) / dot(priFacetNormal, direction));
        vec3 rayPositionAltFacet = rayEntry + direction * (dot(altFacetNormal, altFacetOrigin - rayEntry) / dot(altFacetNormal, direction));
        vec3 rayPosition = isPriFacet ? rayPositionPriFacet : rayPositionAltFacet;

        int scatteringOrder = 0;
        vec3 currentFacetNormal = isPriFacet ? priFacetNormal : altFacetNormal;
        while(++scatteringOrder <= MAX_BSDF_SCATTERING_EVENTS) {
            if(RandNextF() > throughput) {
                throughput = 0.0;
                break;
            }
            throughput /= saturate(throughput);

            direction = SampleOpaqueDielectricFacetBSDF(-direction, currentFacetNormal, albedo, throughput);
            if(dot(currentFacetNormal, direction) < 0.0) break;

            bool escaped;
            if(isPriFacet) {
                bool awayFromAltFacet = dot(altFacetNormal, direction) > 0.0;
                if(awayFromAltFacet) {
                    escaped = true;
                } else {
                    rayPosition += (dot(altFacetNormal, altFacetOrigin - rayPosition) / dot(altFacetNormal, direction)) * direction;
                    height = dot(geometricNormal, rayPosition);
                    escaped = height > 0.0;
                }
            } else {
                bool awayFromPriFacet = dot(priFacetNormal, direction) > 0.0;
                if(awayFromPriFacet) {
                    escaped = true;
                } else {
                    rayPosition += (dot(priFacetNormal, priFacetOrigin - rayPosition) / dot(priFacetNormal, direction)) * direction;
                    height = dot(geometricNormal, rayPosition);
                    escaped = height > 0.0;
                }
            }

            if(escaped) {
                break;
            } else {
                isPriFacet = !isPriFacet;

                currentFacetNormal = isPriFacet ? priFacetNormal : altFacetNormal;
            }
        }

        return direction;
    }

    float VCavity_MaskingShadowng(in vec3 incomingDirection, in vec3 outgoingDirection) {
        vec2 U = RandNext2F();

        vec3 microsurfacePlaneNormal[2];
        vec3 sampledNormal = normalize(incomingDirection + outgoingDirection);
        for (int i = 0; i < 2; i++) {
            //Calculate the microfacet normal for each plane
            float angle = float(i) * (tau / float(2));

            microsurfacePlaneNormal[i] = GetRotationMatrix(vec3(0.0, 0.0, 1.0), angle) * sampledNormal;
        }

        float height = -tan(acos(sampledNormal.z));

        vec3 rayPosition = vec3(sampledNormal.xy == vec2(0.0) ? vec2(0.0) : -normalize(sampledNormal.xy) * RandNextF(), 0.0);

        vec3 oldIncomingDirection = incomingDirection;

        vec3 point = rayPosition + incomingDirection * 3.4e38;

        bool hit = false;
        //Update the hit point on the cavity, and whether we hit anything or not
        float t = 3.4e38;
        int p = RayPlaneIntersection(
            rayPosition, vec3(0, 0, -1),
            vec3(0.0, 0.0, height), sampledNormal,
            t
        );

        if(p != 0) {
            vec3 currentPoint = rayPosition + vec3(0, 0, -1) * t;
            if(currentPoint.z <= 0.0 && currentPoint.z >= height) {
                point = currentPoint;
                hit = true;
            }
        }

        bool inShadow = false;
        bool inView = true;

        rayPosition = point;

        for (int i = 0; i < 2; ++i) {
            //Calculate the visibility of the facet
            float t;
            int p = RayPlaneIntersection(
                rayPosition, incomingDirection,
                vec3(0.0, 0.0, height), microsurfacePlaneNormal[i],
                t
            );

            if(p != 0) {
                vec3 currentPoint = rayPosition + incomingDirection * t;
                if(currentPoint.z <= 0.0 && currentPoint.z >= height) {
                    inView = false;
                }
            }
        }

        for (int i = 0; i < 2; ++i) {
            //Calculate the shadowing of the facet
            float t;
            int p = RayPlaneIntersection(
                rayPosition, outgoingDirection,
                vec3(0.0, 0.0, height), microsurfacePlaneNormal[i],
                t
            );

            if(p != 0) {
                vec3 currentPoint = rayPosition + outgoingDirection * t;
                if(currentPoint.z <= 0.0 && currentPoint.z >= height) {
                    inShadow = true;
                }
            }
        }

        float maskingShadowing = 0.0;
        if(!inShadow) {
            maskingShadowing = 1.0;
        }
        if(!inView) {
            maskingShadowing = 0.0;
        }

        return maskingShadowing;
    }

    float EstimateFacetBSDF(in vec3 incomingDirection, in vec3 outgoingDirection, in vec3 normal, in float albedo) {
        float energyConservationFactor = 1.0;
        uint N = 128u;
        for(uint i = 0u; i < N; ++i) {
            //Calculate the energy conservation factor via brute force integration.
            float nDotV = 1.0 - RandNext2F().x; //This is equivalent to nDotV here.
            float fresnel = FresnelNonPolarized_R(
                                nDotV, 
                                n1, 
                                n2
                            );
            energyConservationFactor -= fresnel * (2.0 * nDotV) / float(N);
        }


        float diffuse  = albedo / pi;
              diffuse *= FresnelNonPolarized_T(dot(incomingDirection, normal), n1, n2);
              diffuse *= FresnelNonPolarized_T(dot(outgoingDirection, normal), n1, n2);
              diffuse /= energyConservationFactor;

        return diffuse;
    }

    float EstimateOpaqueDielectric(in vec3 incomingDirection, in vec3 outgoingDirection, in float albedo, in float roughness) {
        if(outgoingDirection.z < 0.0) {
            return 0.0;
        }
        vec3 oldIncomingDirection = -incomingDirection;

        vec3 geometricNormal = vec3(0.0, 0.0, 1.0);
        vec3 priFacetNormal = StudentT_SampleD(RandNext2F(), gamma, vec2(roughness));
        vec3 altFacetNormal = -reflect(priFacetNormal, geometricNormal);
        vec3 cavityWidthVector = normalize(geometricNormal * dot(geometricNormal, priFacetNormal) / dot(geometricNormal, geometricNormal) - priFacetNormal);
        vec3 priFacetOrigin =  cavityWidthVector;
        vec3 altFacetOrigin = -cavityWidthVector;

        float cavityPosition = RandNextF();
        bool isPriFacet = cavityPosition < LambdaPrimary(priFacetNormal, altFacetNormal, -incomingDirection);

        float height = Height(priFacetNormal, altFacetNormal, priFacetOrigin, altFacetOrigin, incomingDirection, cavityPosition, isPriFacet);

        vec3 rayEntry = mix(priFacetOrigin, altFacetOrigin, cavityPosition);
        vec3 rayPositionPriFacet = rayEntry + incomingDirection * (dot(priFacetNormal, priFacetOrigin - rayEntry) / dot(priFacetNormal, incomingDirection));
        vec3 rayPositionAltFacet = rayEntry + incomingDirection * (dot(altFacetNormal, altFacetOrigin - rayEntry) / dot(altFacetNormal, incomingDirection));
        vec3 rayPosition = isPriFacet ? rayPositionPriFacet : rayPositionAltFacet;

        float shadowHeightPriFacet = dot(priFacetNormal, outgoingDirection) < 0.0 ? 0.0 : dot(geometricNormal, outgoingDirection) * (dot(priFacetNormal, priFacetOrigin - altFacetOrigin) / dot(priFacetNormal, outgoingDirection));
        float shadowHeightAltFacet = dot(altFacetNormal, outgoingDirection) < 0.0 ? 0.0 : dot(geometricNormal, outgoingDirection) * (dot(altFacetNormal, altFacetOrigin - priFacetOrigin) / dot(altFacetNormal, outgoingDirection));

        int scatteringOrder = 0;
        float estimate = 0.0;
        float throughput = 1.0;

        vec3 currentFacetNormal = isPriFacet ? priFacetNormal : altFacetNormal;

        bool inShadow = roughness == 0.0 ? false : height < (isPriFacet ? shadowHeightPriFacet : shadowHeightAltFacet);
        if(!inShadow) {
            float NoL = saturate(dot(outgoingDirection, currentFacetNormal));
            estimate += NoL * EstimateFacetBSDF(-incomingDirection, outgoingDirection, currentFacetNormal, albedo) * throughput;
        }

        while(++scatteringOrder < MAX_BSDF_SCATTERING_EVENTS) {
            if(RandNextF() > throughput) {
                throughput = 0.0;
                break;
            }
            throughput /= saturate(throughput);

            incomingDirection = SampleOpaqueDielectricFacetBSDF(-incomingDirection, currentFacetNormal, albedo, throughput);
            if(dot(currentFacetNormal, incomingDirection) < 0.0) break;

            bool escaped;
            if(isPriFacet) {
                bool awayFromAltFacet = dot(altFacetNormal, incomingDirection) > 0.0;
                if(awayFromAltFacet) {
                    escaped = true;
                } else {
                    rayPosition += (dot(altFacetNormal, altFacetOrigin - rayPosition) / dot(altFacetNormal, incomingDirection)) * incomingDirection;
                    height = dot(geometricNormal, rayPosition);
                    escaped = height > 0.0;
                }
            } else {
                bool awayFromPriFacet = dot(priFacetNormal, incomingDirection) > 0.0;
                if(awayFromPriFacet) {
                    escaped = true;
                } else {
                    rayPosition += (dot(priFacetNormal, priFacetOrigin - rayPosition) / dot(priFacetNormal, incomingDirection)) * incomingDirection;
                    height = dot(geometricNormal, rayPosition);
                    escaped = height > 0.0;
                }
            }

            if(escaped) {
                break;
            } else {
                isPriFacet = !isPriFacet;

                currentFacetNormal = isPriFacet ? priFacetNormal : altFacetNormal;

                bool inShadow = height < (isPriFacet ? shadowHeightPriFacet : shadowHeightAltFacet);
                if(!inShadow) {
                    float NoL = saturate(dot(outgoingDirection, currentFacetNormal));
                    estimate += NoL * EstimateFacetBSDF(-incomingDirection, outgoingDirection, currentFacetNormal, albedo) * throughput;
                }
            }
        }

        vec3 halfway = normalize(oldIncomingDirection + outgoingDirection);
        float F = FresnelNonPolarized_R(dot(oldIncomingDirection, halfway), n1, n2);
        float G = VCavity_MaskingShadowng(oldIncomingDirection, outgoingDirection);
        float numerator = StudentT_D(halfway, gamma, vec2(roughness)) * F * G;
        float denominator = 4.0 * outgoingDirection.z * oldIncomingDirection.z;
        float specular = numerator / denominator;

        return (estimate / outgoingDirection.z) + max(specular, 0.0);
    }
#endif