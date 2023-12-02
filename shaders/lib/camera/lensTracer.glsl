struct LensInterface {
    float axisPosition; // Position along the lens axis
    float axisRadius; // Radius around the lens axis

    int shapeType; // 0 = Ring, 1 = Sphere, 2 = Disc
    float curvatureRadius; // Defines the curvature radius, lower is higher curvature

	bool hasThinCoating; // Enable thin film coating
	float coatThickness; // Thickness of the thin film
};

//#define IDEAL_FILM_IOR
#define USE_POLYGONAL_APERTURE

#define WIDE_FOV_LENS_APERTURE_SIZE 5.0 //[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0]

//#define ENABLE_ANTI_REFLECTIVE_COATING

#ifdef ENABLE_ANTI_REFLECTIVE_COATING
bool arCoating = true;
#else
bool arCoating = false;
#endif

const int num_interfaces = 5;

const vec2 sensorMaxRadii = vec2(127.0, 101.0) / 2.0;

void WideFoVLensInterfaces(in float wavelength, out LensInterface lens[5], out float N[6]) {
    lens = LensInterface[](
        LensInterface(  0.0, 20.0, 1,  21.0, arCoating, 500.0),
        LensInterface(  4.0, 20.0, 1,  22.5, arCoating, 600.0),

        LensInterface( 20.0, WIDE_FOV_LENS_APERTURE_SIZE, 2,1.0/0.0, arCoating, 650.0),

        LensInterface( 36.0, 20.0, 1, -22.5, arCoating, 400.0),
        LensInterface( 40.0, 20.0, 1, -21.0, arCoating, 300.0)
    );

    N = float[](
        Air(wavelength * 1e-3),

        F2(wavelength * 1e-3),

        Air(wavelength * 1e-3),
        Air(wavelength * 1e-3),

        F2(wavelength * 1e-3),

        Air(wavelength * 1e-3)
    );
}

// Gets the 2x2 ray transfer matrix of a single lens interface.
mat2 GetInterfaceTransferMatrix(LensInterface lensInterface, float n1, float n2) {
	mat2 transfer_matrix;
	if (lensInterface.shapeType == 1) {
		transfer_matrix = mat2(
			1.0, (n1 - n2) / (lensInterface.curvatureRadius * n2),
			0.0, n1 / n2
		);
	} else {
		transfer_matrix = mat2(
			1.0, 0.0,
			0.0, n1 / n2
		);
	}
	return transfer_matrix;
}
// Gets the 2x2 ray transfer matrix of the entire lens system.
// This can be used to solve for, for example, the focus plane.
mat2 GetLensTransferMatrix(in LensInterface interfaces[num_interfaces], in float refractiveIndices[num_interfaces + 1]) {
	// first interface
	float n1 = refractiveIndices[0];
	float n2 = refractiveIndices[1];
	mat2 lens_transfer_matrix = GetInterfaceTransferMatrix(interfaces[0], n1, n2);

	for (int i = 1; i < num_interfaces; ++i) {
		// Account for separation
		float separation = interfaces[i].axisPosition - interfaces[i - 1].axisPosition;
		mat2 separation_transfer_matrix = mat2(
			1.0, 0.0,
			separation, 1.0
		);
		lens_transfer_matrix = separation_transfer_matrix * lens_transfer_matrix;

		// Account for interface
		float n1 = refractiveIndices[i];
		float n2 = refractiveIndices[i + 1];
		mat2 interface_transfer_matrix = GetInterfaceTransferMatrix(interfaces[i], n1, n2);
		lens_transfer_matrix = interface_transfer_matrix * lens_transfer_matrix;
	}

	return lens_transfer_matrix;
}
mat2 GetLensTransferMatrix(in float wavelength) {
    LensInterface[num_interfaces] interfaces;
    float N[num_interfaces + 1];

    WideFoVLensInterfaces(wavelength, interfaces, N);

    return GetLensTransferMatrix(interfaces, N);
}
mat2 GetLensTransferMatrix() {
    LensInterface[num_interfaces] interfaces;
    float N[num_interfaces + 1];

    WideFoVLensInterfaces(587.6, interfaces, N);

    return GetLensTransferMatrix(interfaces, N);
}

// Using the ray transfer matrix, finds the plane that objects at a specified distance are focused onto.
float FindFocusPlanePosition(in mat2 transferMatrix, in float focusDepth, in float axisPosition) {
	// Simpler equivalent, for reference:
	// vec2 outRay = transferMatrix * vec2(1.0, 1.0 / focusDepth);
	// float t = outRay.x / outRay.y;

	float t;
	// Use one of 2 versions to avoid issues.
	// Each becomes an inf/inf depending on the value of focusDepth.
	// One at focusDepth = 0, the other at focusDepth = infinity.
	// (0/0 is also possible depending on details of the matrix)
	if (focusDepth < 1.0) {
		t = (focusDepth * transferMatrix[0].x + transferMatrix[1].x)
		  / (focusDepth * transferMatrix[0].y + transferMatrix[1].y);
	} else {
		float invFocusDepth = 1.0 / focusDepth;
		t = (transferMatrix[0].x + transferMatrix[1].x * invFocusDepth)
		  / (transferMatrix[0].y + transferMatrix[1].y * invFocusDepth);
	}
	return axisPosition - t;
}

vec3 GeneratePointOnFilm(in float sensorPosition) {
    vec2 pixel = vec2(gl_FragCoord.xy) + RandNext2F();
    vec2 sensorPoint = (vec2(viewWidth, viewHeight) - 2.0 * pixel)
                     * max(sensorMaxRadii.x / viewWidth, sensorMaxRadii.y / viewHeight);
    return vec3(2.0 * sensorPoint, sensorPosition);
}

vec3 SampleBackLensInterface(LensInterface lensInterface, vec2 xy) {
    if(lensInterface.shapeType == 1) {
        float sinAngle = lensInterface.axisRadius / -lensInterface.curvatureRadius;
        float cosAngle = cos(asin(sinAngle));

        xy.x *= radians(360.0);
        xy.y = xy.y * (1.0 - cosAngle) + cosAngle;

        vec3 shape = vec3(vec2(cos(xy.x), sin(xy.x)) * sqrt(1.0 - xy.y * xy.y), xy.y - 1.0);
        return -lensInterface.curvatureRadius * shape + vec3(0.0, 0.0, lensInterface.axisPosition);
    } else {
        return vec3(lensInterface.axisRadius * (vec2(cos(tau * xy.x), sin(tau * xy.x)) * sqrt(xy.y)), lensInterface.axisPosition);
    }
}

void CalculateLensSystem(out vec3 rayDirection, out vec3 rayPosition, out float throughput, out bool invalid, in float focusDistance, in float wavelength) {
    LensInterface[num_interfaces] interfaces;
    float N_Lens[num_interfaces + 1];

    WideFoVLensInterfaces(wavelength, interfaces, N_Lens);

    mat2 T = GetLensTransferMatrix();

    float sensorPosition = FindFocusPlanePosition(T, focusDistance * 1000.0, interfaces[num_interfaces - 1].axisPosition);

    rayPosition = GeneratePointOnFilm(sensorPosition);

    LensInterface backInterface = interfaces[num_interfaces - 1];

    vec3 interfacePoint = SampleBackLensInterface(backInterface, RandNext2F());
    vec3 interfaceNormal = normalize(interfacePoint - vec3(0.0, 0.0, backInterface.axisPosition - -backInterface.curvatureRadius));
    vec3 interfaceFromSensor = interfacePoint - rayPosition;
    rayDirection = normalize(interfaceFromSensor);

    float interfaceArea = 2.0 * pi * (1.0 - cos(asin(backInterface.axisRadius / -backInterface.curvatureRadius))) * -backInterface.curvatureRadius * -backInterface.curvatureRadius;
    float interfaceWeight = interfaceArea * dot(-rayDirection, interfaceNormal) / dot(interfaceFromSensor, interfaceFromSensor);
    float sensorWeight = dot(rayDirection, vec3(0.0, 0.0, -1.0));
    throughput = sensorWeight * interfaceWeight;

    invalid = false;
    int bounces = 0;
    int interfaceIndex = num_interfaces - 1;
    int direction = -1;
    do {
        if(bounces > num_interfaces * 20) break;
        LensInterface interfaceData = interfaces[interfaceIndex];

        vec3 interfaceNormal;
        if (interfaceData.shapeType == 1) {
            vec3 sphereOrigin = vec3(0.0, 0.0, interfaceData.axisPosition - -interfaceData.curvatureRadius);

            float t0, t1;
            int intersectionType = RaySphereIntersection(
                rayPosition, rayDirection,
                sphereOrigin, -interfaceData.curvatureRadius,
                t0, t1
            );

            if (intersectionType == 0) {
                invalid = true;
                break;
            } else if (intersectionType == 1) {
                if (-direction * -interfaceData.curvatureRadius > 0.0) {
                    invalid = true;
                    break;
                }

                rayPosition += rayDirection * t0;
            } else if (intersectionType == 2) {
                if (-direction * -interfaceData.curvatureRadius < 0.0) {
                    rayPosition += rayDirection * t1;
                } else {
                    rayPosition += rayDirection * t0;
                }
            } else {
                // This should be unreachable.
                invalid = true;
                break;
            }

            if(length(rayPosition.xy) > interfaceData.axisRadius) invalid = true;

            interfaceNormal = normalize(rayPosition - sphereOrigin);
        } else if (interfaceData.shapeType == 0) {
            interfaceNormal = vec3(0.0, 0.0, 1.0);

            float t;
            int intersectionType = RayPlaneIntersection(
                rayPosition, rayDirection,
                vec3(0.0, 0.0, interfaceData.axisPosition), interfaceNormal,
                t
            );

            if (intersectionType == 0) {
                invalid = true;
                break;
            } else if (intersectionType == 1) {
                rayPosition += rayDirection * t;
            } else if (intersectionType == 2) {
                invalid = true;
                break;
            } else {
                // This should be unreachable.
                invalid = true;
                break;
            }
        } else if (interfaceData.shapeType == 2) {
            interfaceNormal = vec3(0.0, 0.0, 1.0);

            #ifndef USE_POLYGONAL_APERTURE
                float t = RDI(rayPosition - vec3(0.0, 0.0, interfaceData.axisPosition), rayDirection, interfaceData.axisRadius);
                bool intersection = t > 0.0001;
            #else
                float t;
                bool intersection = IntersectAperture(rayPosition - vec3(0.0, 0.0, interfaceData.axisPosition), rayDirection, interfaceData.axisRadius, t);
            #endif
            int intersectionType = intersection ? 1 : 0;

            if (intersectionType == 0) {
                invalid = true;
                break;
            } else if (intersectionType == 1) {
                rayPosition += rayDirection * t;
            } else if (intersectionType == 2) {
                invalid = true;
                break;
            } else {
                // This should be unreachable.
                invalid = true;
                break;
            }
        }

        float n1;
        float n2;
        if(direction == -1) {
            n1 = N_Lens[interfaceIndex + 1];
            n2 = N_Lens[interfaceIndex];
        } else if(direction == 1) {
            n1 = N_Lens[interfaceIndex];
            n2 = N_Lens[interfaceIndex + 1];
        }

        if (dot(rayDirection, interfaceNormal) > 0.0) {
            interfaceNormal = -interfaceNormal;
        }

        float fresnel_R;
        float fresnel_T;
        if(interfaceData.hasThinCoating) {
            #ifdef IDEAL_FILM_IOR
                vec2 filmIOR = vec2(sqrt(n1 * n2), 0.0);
            #else
                vec2 filmIOR = Silica(wavelength);
            #endif
            //The divide by 2 should hopefully result in a pi/2 phase change
            fresnel_R = FresnelThinFilmInterferenceReflected(dot(interfaceNormal,-rayDirection), interfaceData.coatThickness / 2.0, wavelength, vec2(n1, 0.0), filmIOR, vec2(n2, 0.0));
            fresnel_T = FresnelThinFilmInterferenceTransmitted(dot(interfaceNormal,-rayDirection), interfaceData.coatThickness / 2.0, wavelength, vec2(n1, 0.0), filmIOR, vec2(n2, 0.0));
        } else {
            fresnel_R = FresnelNonPolarized_R(dot(interfaceNormal,-rayDirection), vec2(n1, 0.0), vec2(n2, 0.0));
            fresnel_T = FresnelNonPolarized_T(dot(interfaceNormal,-rayDirection), vec2(n1, 0.0), vec2(n2, 0.0));
        }

        float specBounceProbability = fresnel_R;

        bool specularBounce = specBounceProbability > RandNextF();
    
        if(specularBounce) {
            direction = -direction;

            throughput /= specBounceProbability;

            throughput *= fresnel_R;

            rayDirection = reflect(rayDirection, interfaceNormal);        
        } else {
            throughput /= 1.0 - specBounceProbability;

            throughput *= fresnel_T;

            rayDirection = refract(rayDirection, interfaceNormal, n1 / n2);
        }

        //Afaik this branch shouldn't do anything if fresnel is working properly for TIR, but for some reason this fixes an error 
        if (rayDirection == vec3(0.0)) {
            invalid = true;
            break;
        }

        if(isinf(throughput)) {
            invalid = true;
            break;
        }
        if(isnan(throughput)) {
            invalid = true;
            break;
        }
        if(any(isnan(rayDirection))) {
            invalid = true;
            break;
        }
        if(any(isinf(rayDirection))) {
            invalid = true;
            break;
        }

        interfaceIndex += direction;
        bounces += 1;
    } while(interfaceIndex >= 0 && interfaceIndex < num_interfaces);

    if (direction > 0 || interfaceIndex >= 0) {
        invalid = true;
    }

    if(isnan(throughput)) throughput = 0.0;
    if(isinf(throughput)) throughput = 1.0;
}