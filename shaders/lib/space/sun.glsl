#if !defined LIB_SPACE_SUN
#define LIB_SPACE_SUN
    float PhysicalSun(in vec3 sceneDirection, in vec3 lightDirection, in float wavelength, in float radiance) {
        /*
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_u.pro
            https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_v.pro
        */
        const float au = -8.9829751;
        const float bu = 0.0069093916;
        const float cu = -1.8144591e-6;
        const float du = 2.2540875e-10;
        const float eu = 2.2540875e-10;
        const float fu = 3.0453572e-19;

        const float av = 9.2891180;
        const float bv = -0.0062212632;
        const float cv = 1.5788029e-6;
        const float dv = -1.9359644e-10;
        const float ev = 1.1444469e-14;
        const float fv = -2.599494e-19;

        float cosTheta = dot(sceneDirection, lightDirection);
        float cosSunRadius = cos(sunAngularRadius);
        if(cosTheta < cosSunRadius) return 0.0;

        float lambda = wavelength * 10.0;
        float u = au + (bu * lambda) + (cu * pow(lambda, 2.0)) + (du * pow(lambda, 3.0)) + (eu * pow(lambda, 4.0)) + (fu * pow(lambda, 5.0));
        float v = av + (bv * lambda) + (cv * pow(lambda, 2.0)) + (dv * pow(lambda, 3.0)) + (ev * pow(lambda, 4.0)) + (fv * pow(lambda, 5.0));

        float darkening = 1.0 - u - v + u*cos(asin(max(0.0, 1.0 - cosTheta))) + v*pow(cos(asin(max(0.0, 1.0 - cosTheta))), 2.0);
        float finalLuminance = radiance / darkening;

        return finalLuminance;
    }
#endif