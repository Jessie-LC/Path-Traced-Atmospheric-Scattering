#if !defined LIB_ATMOSPHERE_CONSTANTS
#define LIB_ATMOSPHERE_CONSTANTS
    const float kilometer = 1000.0;
    const float airNumberDensity   = 2.5035422e25;
    const float ozonePeakDensity   = 3e-6;
    const float ozoneNumberDensity = airNumberDensity * exp(-35e3 / 8e3) * (134.628/48.0) * ozonePeakDensity;
    const float planetRadius       = 6371e3;
    const float atmosphereHeight   = 100e3;

    const vec2 scaleHeights = vec2(8.0, 1.2) * kilometer;

    const float aerosol_g = 0.76;
    const float meanAerosolParticleDiameter = 0.8;

    const vec2 inverseScaleHeights   = 1.0 / (scaleHeights);
    const vec2 scaledPlanetRadius    = planetRadius * inverseScaleHeights;
    const float atmosphereRadius     = planetRadius + atmosphereHeight;
    const float atmosphereRadiusSquared = atmosphereRadius*atmosphereRadius; 
    const float atmosphereLowerLimit = planetRadius + -1.0;
    const float atmosphereLowerLimitSquared = atmosphereLowerLimit*atmosphereLowerLimit;

    const float aerosolScatteringAlbedo = 0.98;
    const float mistScatteringAlbedo = 0.99;
    const float cloudScatteringAlbedo = 0.99;

    const float cloudsAltitude = 600.0;
    const float cloudsThickness = cloudsAltitude * 0.1;
    const float cloudsMaxAltitude = cloudsAltitude + cloudsThickness;
    const float globalCoverage = 0.8;
    const float cloudAnvilAmount = 0.0;

    const float rAir = 287.053;
    const float gammaAir = 1.4;
#endif