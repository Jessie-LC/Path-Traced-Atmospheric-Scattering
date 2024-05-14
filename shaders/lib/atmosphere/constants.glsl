#if !defined LIB_ATMOSPHERE_CONSTANTS
#define LIB_ATMOSPHERE_CONSTANTS
    const float kilometer = 1000.0;
    const float airNumberDensity   = 2.5035422e25;
    const float planetRadius       = 6371e3;
    const float atmosphereHeight   = 400e3;

    const vec2 scaleHeights = vec2(8.0, 1.4) * kilometer;

    const float aerosol_g = 0.76;
    const float meanAerosolParticleDiameter = 0.8;

    const vec2 inverseScaleHeights   = 1.0 / (scaleHeights);
    const vec2 scaledPlanetRadius    = planetRadius * inverseScaleHeights;
    const float atmosphereRadius     = planetRadius + atmosphereHeight;
    const float atmosphereRadiusSquared = atmosphereRadius*atmosphereRadius; 
    const float atmosphereLowerLimit = planetRadius;
    const float atmosphereLowerLimitSquared = atmosphereLowerLimit*atmosphereLowerLimit;

    const float aerosolScatteringAlbedo = 0.90; //Clean aerosols scatter around 99.9% of light, however dirty aerosols scatter around 90%
    const float mistScatteringAlbedo = 0.99;
    const float cloudScatteringAlbedo = 0.99;

    const float cloudsAltitude = 600.0;
    const float cloudsThickness = cloudsAltitude * 0.8;
    const float cloudsMaxAltitude = cloudsAltitude + cloudsThickness;
    const float globalCoverage = 0.7;
    const float cloudAnvilAmount = 0.0;

    const float rAir = 287.053;
    const float gammaAir = 1.4;
#endif