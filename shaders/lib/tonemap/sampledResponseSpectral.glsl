#if !defined TONEMAP_SAMPLED_RESPONSE
#define TONEMAP_SAMPLED_RESPONSE
    float GetOrOne(in float wavelength, in ivec2 index) {
        vec4 LUT = texelFetch(cameraResponseLUT, index, 0);
        float spectrum;
        RGBToSpectrum(spectrum, wavelength, LUT.r, LUT.g, LUT.b, 0);
        return index.x < 1024.0 || index.y < 1024.0 ? spectrum : 1.0;
    }

    int BinarySearch(in float wavelength, in float target) {
        int first = 0;
        int last = 1024 - 1;

        while(first <= last) {
            int guess = int(floor(float(last + first) / 2.0));

            vec4 LUT = texelFetch(cameraResponseLUT, ivec2(guess, 0), 0);
            float spectrum;
            RGBToSpectrum(spectrum, wavelength, LUT.r, LUT.g, LUT.b, 1);
            if(spectrum == target) {
                return guess;
            } else if(spectrum < target) {
                first = guess + 1;
            } else {
                last = guess - 1;
            }
        }

        return -1;
    }

    float CameraGetIntensity(in float f, in float iso, in float wavelength) {
        f = clamp(f, 0.0, iso);
        #if CAMERA_RESPONSE == 1
            f /= 4.0;
        #elif CAMERA_RESPONSE == 0
            f /= 20.0;
        #elif CAMERA_RESPONSE == 2
            f /= 2.5;
        #elif CAMERA_RESPONSE == 3
            f /= 6.0;
        #elif CAMERA_RESPONSE == 4
            f /= 12.0;
        #endif

        if(isnan(f)) {
            f = 0.0;
        }
        if(isinf(f)) {
            f = 3.4e38;
        }

        vec4 LUT = texelFetch(cameraResponseLUT, ivec2(0, 0), 0);
        float spectrum;
        RGBToSpectrum(spectrum, wavelength, LUT.r, LUT.g, LUT.b, 0);

        float upper = float(BinarySearch(wavelength, f));
        float idx = distance(spectrum, upper);

        LUT = texelFetch(cameraResponseLUT, ivec2(int(idx), 0), 0);
        RGBToSpectrum(spectrum, wavelength, LUT.r, LUT.g, LUT.b, 1);

        float lowIrradiance = spectrum;
        float highIrradiance = GetOrOne(wavelength, ivec2(int(idx + 1), 0));
        float mixParameter = (f - lowIrradiance) / (highIrradiance - lowIrradiance);

        LUT = texelFetch(cameraResponseLUT, ivec2(0, int(idx)), 0);
        RGBToSpectrum(spectrum, wavelength, LUT.r, LUT.g, LUT.b, 0);

        float lowValue = spectrum;
        float highValue = GetOrOne(wavelength, ivec2(0, int(idx + 1)));

        return saturate(mix(lowValue, highValue, mixParameter));
    }

    vec3 CameraTonemap(in vec3 v, in float iso) {
        const int numberOfWavelengths = 32;
        vec3 xyz = vec3(0.0);
        for(int i = 0; i < numberOfWavelengths; ++i) {
            float wavelength = (float(i + 1) / numberOfWavelengths) * 441.0 + 390.0;

            float spectrum;
            RGBToSpectrum(spectrum, wavelength, v.r, v.g, v.b, 1);
            spectrum = CameraGetIntensity(spectrum, iso, wavelength);
            xyz += SpectrumToXYZExact_CIE2012(spectrum / 683.368, wavelength) / numberOfWavelengths;
        }
        return xyz * xyzToRGBMatrix_D65;
    }
#endif