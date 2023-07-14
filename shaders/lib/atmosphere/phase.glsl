#if !defined LIB_ATMOSPHERE_PHASE
#define LIB_ATMOSPHERE_PHASE
    float BinarySearchPhase(sampler1D cdfTexture, float lowBound, float highBound, float toFind) {
        float midIndex = (lowBound + highBound) / 2.0;

        int lookupRes = textureSize(cdfTexture, 0);
        for (int x = 0; x < 12; ++x) {
            float coordinate = midIndex / pi;
            float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / float(lookupRes);

            float value = saturate(texture(cdfTexture, coordTweaked).r); //The CDF value should never go above 1, but clamp it anyway just in case

            if (value < toFind) {
                lowBound = midIndex;
            }
            else if (value > toFind) {
                highBound = midIndex;
            }
            else {
                return midIndex;
            }

            midIndex = (lowBound + highBound) / 2.0;
        }

        return midIndex;
    }
    float BinarySearchPhase(sampler2D cdfTexture, float wavelength, float lowBound, float highBound, float toFind) {
        float midIndex = (lowBound + highBound) / 2.0;

        int lookupRes = textureSize(cdfTexture, 0).x;
        for (int x = 0; x < 12; ++x) {
            float coordinate = midIndex / pi;
            float coordTweaked = (coordinate * float(lookupRes - 1) + 0.5) / float(lookupRes);

            float value = saturate(texture(cdfTexture, vec2(coordTweaked, (wavelength - 390.0) / 441.0)).r); //The CDF value should never go above 1, but clamp it anyway just in case

            if (value < toFind) {
                lowBound = midIndex;
            }
            else if (value > toFind) {
                highBound = midIndex;
            }
            else {
                return midIndex;
            }

            midIndex = (lowBound + highBound) / 2.0;
        }

        return midIndex;
    }

    vec3 SampleCloudPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureCloud, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleAerosolPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureAerosol, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleLowAltitudeAerosolPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureAerosolLowAltitude, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleRainbowPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRainbow, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }
    vec3 SampleRayleighPhase(in float wavelength) {
        float cosTheta = cos(BinarySearchPhase(cdfTextureRayleigh, wavelength, 0.0, pi, RandNextF()));
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }

    float CloudPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureCloud, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float AerosolPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureAerosol, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float LowAltitudeAerosolPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureAerosolLowAltitude, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RainbowPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureRainbow, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }
    float RayleighPhase(in float cosTheta, in float wavelength) {
        return texture(phaseTextureRayleigh, vec2(acos(cosTheta) / pi, (wavelength - 390.0) / 441.0)).r;
    }

    float HenyeyGreensteinPhase(float cosTheta, float g) {
        const float norm = 0.25/pi;

        float gg = g * g;
        return norm * ((1.0-gg) / pow(1.0+gg-2.0*g*cosTheta, 3.0/2.0));
    }
	vec3 SampleHenyeyGreensteinPhase(float P, float g) {
		float s = 2.0 * P - 1.0;
		float u = (0.5 / g) * (1.0 + g * g - pow((1.0 - g * g) / (1.0 + g * s), 2.0));

		return GenerateUnitVector(vec2(RandNextF(), u * 0.5 + 0.5));
	}

    float KleinNishinaPhase(float cosTheta, float g) {
        float e = 1.0;
        for (int i = 0; i < 8; ++i) {
            float gFromE = 1.0 / e - 2.0 / log(2.0 * e + 1.0) + 1.0;
            float deriv = 4.0 / ((2.0 * e + 1.0) * square(log(2.0 * e + 1.0))) - 1.0 / square(e);
            if (abs(deriv) < 0.00000001) break;
            e = e - (gFromE - g) / deriv;
        }

        return e / (2.0 * pi * (e * (1.0 - cosTheta) + 1.0) * log(2.0 * e + 1.0));
    }
    vec3 SampleKleinNishinaPhase(float g) {
        float e = 1.0;
        for (int i = 0; i < 8; ++i) {
            float gFromE = 1.0 / e - 2.0 / log(2.0 * e + 1.0) + 1.0;
            float deriv = 4.0 / ((2.0 * e + 1.0) * square(log(2.0 * e + 1.0))) - 1.0 / square(e);
            if (abs(deriv) < 0.00000001) break;
            e = e - (gFromE - g) / deriv;
        }

        float cosTheta = (-pow(2.0 * e + 1.0, 1.0 - RandNextF()) + e + 1.0) / e;
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }

    float DrainePhase(in float u, in float g, in float a) {
        return ((1 - g*g)*(1 + a*u*u))/(4.*(1 + (a*(1 + 2*g*g))/3.) * pi * pow(1 + g*g - 2*g*u,1.5));
    }
    float SampleDrainePhase(in float xi, in float g, in float a)
    {
        const float g2 = g * g;
        const float g3 = g * g2;
        const float g4 = g2 * g2;
        const float g6 = g2 * g4;
        const float pgp1_2 = (1 + g2) * (1 + g2);
        const float T1 = (-1 + g2) * (4 * g2 + a * pgp1_2);
        const float T1a = -a + a * g4;
        const float T1a3 = T1a * T1a * T1a;
        const float T2 = -1296 * (-1 + g2) * (a - a * g2) * (T1a) * (4 * g2 + a * pgp1_2);
        const float T3 = 3 * g2 * (1 + g * (-1 + 2 * xi)) + a * (2 + g2 + g3 * (1 + 2 * g2) * (-1 + 2 * xi));
        const float T4a = 432 * T1a3 + T2 + 432 * (a - a * g2) * T3 * T3;
        const float T4b = -144 * a * g2 + 288 * a * g4 - 144 * a * g6;
        const float T4b3 = T4b * T4b * T4b;
        const float T4 = T4a + sqrt(-4 * T4b3 + T4a * T4a);
        const float T4p3 = pow(T4, 1.0 / 3.0);
        const float T6 = (2 * T1a + (48 * pow(2, 1.0 / 3.0) *
            (-(a * g2) + 2 * a * g4 - a * g6)) / T4p3 + T4p3 / (3. * pow(2, 1.0 / 3.0))) / (a - a * g2);
        const float T5 = 6 * (1 + g2) + T6;
        return (1 + g2 - pow(-0.5 * sqrt(T5) + sqrt(6 * (1 + g2) - (8 * T3) / (a * (-1 + g2) * sqrt(T5)) - T6) / 2., 2)) / (2. * g);
    }
    vec3 SampleDrainePhase(in float g, in float a) {
        float cosTheta = SampleDrainePhase(RandNextF(), g, a);
        return GenerateUnitVector(vec2(RandNextF(), cosTheta * 0.5 + 0.5));
    }

    float ApproximateMiePhase(in float cosTheta, in float d) {
        float g_hg;
        float g_d;
        float a;
        float wd;
        if(50.0 > d && d >= 5.0) {
            g_hg = exp(-0.990567 / (d - 1.67154));
            g_d  = exp(-2.20679 / (d + 3.91029)) - 0.428934;
            a    = exp(3.62489) - (8.29288 / (d + 5.52825));
            wd   = exp(-0.599085 / (d - 0.641583)) - 0.665888;
        } else if(d >= 1.5 && d < 5.0) {
            g_hg = 0.0604931 * log(log(d)) + 0.940256;
            g_d  = 0.500411 - (0.081287 / (-2.0 * log(d) + tan(log(d)) + 1.27551));
            a    = 7.30354 * log(d) + 6.31675;
            wd   = 0.026914 * (log(d) - cos(5.68947 * (log(log(d))- 0.0292149))) + 0.376475;
        } else if(d < 1.5 && d > 0.1) {
            float logD = log(d);
            float num = (logD - 0.238604) * (logD + 1.00667);
            float denom = 0.507522 - logD * 0.15677;
            float a = 1.19692 * cos(num / denom) + logD * 1.37932 + 0.0625835;

            g_hg = 0.862 - 0.143 * pow(log(d), 2.0);
            g_d  = 0.379685 * cos(a) + 0.344213;
            a    = 250.0;
            wd   = 0.146209 * cos(3.38707 * log(d) + 2.11193) + 0.316072 + 0.0778917 * log(d);
        } else if(d <= 0.1) {
            g_hg = 13.8 * pow(d, 2.0);
            g_d  = 1.1456 * d * sin(9.29044 * d);
            a    = 250.0;
            wd   = 0.252977 - 312.983 * pow(d, 4.3);
        }

        float henyeyGreenstein = (1.0 - wd) * HenyeyGreensteinPhase(cosTheta, g_hg);
        float draine = wd * DrainePhase(cosTheta, g_d, a);

        return henyeyGreenstein + draine;
    }
    vec3 SampleApproximateMiePhase(in float d) {
        float g_hg;
        float g_d;
        float a;
        float wd;
        if(50.0 > d && d >= 5.0) {
            g_hg = exp(-0.990567 / (d - 1.67154));
            g_d  = exp(-2.20679 / (d + 3.91029)) - 0.428934;
            a    = exp(3.62489) - (8.29288 / (d + 5.52825));
            wd   = exp(-0.599085 / (d - 0.641583)) - 0.665888;
        } else if(d >= 1.5 && d < 5.0) {
            g_hg = 0.0604931 * log(log(d)) + 0.940256;
            g_d  = 0.500411 - (0.081287 / (-2.0 * log(d) + tan(log(d)) + 1.27551));
            a    = 7.30354 * log(d) + 6.31675;
            wd   = 0.026914 * (log(d) - cos(5.68947 * (log(log(d))- 0.0292149))) + 0.376475;
        } else if(d < 1.5 && d > 0.1) {
            float logD = log(d);
            float num = (logD - 0.238604) * (logD + 1.00667);
            float denom = 0.507522 - logD * 0.15677;
            float a = 1.19692 * cos(num / denom) + logD * 1.37932 + 0.0625835;

            g_hg = 0.862 - 0.143 * pow(log(d), 2.0);
            g_d  = 0.379685 * cos(a) + 0.344213;
            a    = 250.0;
            wd   = 0.146209 * cos(3.38707 * log(d) + 2.11193) + 0.316072 + 0.0778917 * log(d);
        } else if(d <= 0.1) {
            g_hg = 13.8 * pow(d, 2.0);
            g_d  = 1.1456 * d * sin(9.29044 * d);
            a    = 250.0;
            wd   = 0.252977 - 312.983 * pow(d, 4.3);
        }

        if(RandNextF() < 1.0 - wd) {
            return SampleHenyeyGreensteinPhase(RandNextF(), g_hg);
        }
        if(RandNextF() < wd) {
            return SampleDrainePhase(g_d, a);
        }
    }
#endif