#if !defined TONEMAP_SAMPLED_RESPONSE
#define TONEMAP_SAMPLED_RESPONSE
    float GetOrOne(in int component, in ivec2 index) {
        return index.x < 1024.0 || index.y < 1024.0 ? texelFetch(cameraResponseLUT, index, 0)[component] : 1.0;
    }

    int BinarySearch(in int component, in float target) {
        int first = 0;
        int last = 1024 - 1;

        while(first <= last) {
            int guess = int(floor(float(last + first) / 2.0));
            if(texelFetch(cameraResponseLUT, ivec2(guess, 0), 0)[component] == target) {
                return guess;
            } else if(texelFetch(cameraResponseLUT, ivec2(guess, 0), 0)[component] < target) {
                first = guess + 1;
            } else {
                last = guess - 1;
            }
        }

        return -1;
    }

    float CameraGetIntensity(in float f, in float iso, in int component) {
        f = clamp(f, 0.0, iso);
        #if CAMERA_RESPONSE == 1
            f /= 4.0;
        #elif CAMERA_RESPONSE == 0
            f /= 20.0;
        #elif CAMERA_RESPONSE == 2
            f /= 2.5;
        #elif CAMERA_RESPONSE == 3
            f /= 6.0;
        #endif

        float upper = float(BinarySearch(component, f));
        float idx = distance(texelFetch(cameraResponseLUT, ivec2(0, 0), 0)[component], upper);

        float lowIrradiance = texelFetch(cameraResponseLUT, ivec2(int(idx), 0), 0)[component];
        float highIrradiance = GetOrOne(component, ivec2(int(idx + 1), 0));
        float mixParameter = (f - lowIrradiance) / (highIrradiance - lowIrradiance);

        float lowValue = texelFetch(cameraResponseLUT, ivec2(0, int(idx)), 0)[component];
        float highValue = GetOrOne(component, ivec2(0, int(idx + 1)));

        return saturate(mix(lowValue, highValue, mixParameter));
    }

    vec3 CameraTonemap(in vec3 v, in float iso) {
        float r = CameraGetIntensity(v.r, iso, 0);
        float g = CameraGetIntensity(v.g, iso, 1);
        float b = CameraGetIntensity(v.b, iso, 2);
        return vec3(r, g, b);
    }
#endif