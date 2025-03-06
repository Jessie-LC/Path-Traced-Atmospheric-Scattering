#if !defined LIB_ATMOSPHERE_CONSTANTS
#define LIB_ATMOSPHERE_CONSTANTS
    const float kilometer = 1000.0;
    const float earthRadius        = 6371e3;
    const float airNumberDensity   = 2.5035422e25;
    const float planetRadius       = earthRadius;
    const float atmosphereHeight   = 120e3;
    const float lutHeight          = 120e3; // This is maximum altitude of the LUT

    const vec2 scaleHeights = vec2(8.0, 1.4) * kilometer;

    const float aerosol_g = 0.76;
    const float meanAerosolParticleDiameter = 0.8;
    const float meanCloudParticleDiameter = 5.0;

    const vec2 inverseScaleHeights          = 1.0 / (scaleHeights);
    const vec2 scaledPlanetRadius           = planetRadius * inverseScaleHeights;
    const float atmosphereRadius            = planetRadius + atmosphereHeight;
    const float atmosphereRadiusSquared     = atmosphereRadius*atmosphereRadius;
    const float lutRadius                   = planetRadius + lutHeight;
    const float lutRadiusSquared            = lutRadius*lutRadius;
    const float atmosphereLowerLimit        = planetRadius;
    const float atmosphereLowerLimitSquared = atmosphereLowerLimit*atmosphereLowerLimit;

    const float cleanAerosolAlbedo = 0.99;
    const float dirtyAerosolAlbedo = 0.90;
    const float aerosolScatteringAlbedo = dirtyAerosolAlbedo; // Clean aerosols scatter around 99.9% of light, however dirty aerosols scatter around 90%
    const float mistScatteringAlbedo = 0.99;
    const float cloudScatteringAlbedo = 0.99;

    const float cloudsAltitude = 600.0;
    const float cloudsThickness = cloudsAltitude * 0.3;
    const float cloudsMaxAltitude = cloudsAltitude + cloudsThickness;
    const float globalCoverage = 0.5;
    const float cloudAnvilAmount = 0.0;

    const float rAir = 287.053;
    const float gammaAir = 1.4;

    const float g = 9.80665;

    // Units are inverse kilometers
    const float AerosolCoefficient[] = {
        2.60e-1,
        2.50e-1,
        2.40e-1,
        2.40e-1,
        2.30e-1,
        2.00e-1,
        1.80e-1,
        1.67e-1,
        1.58e-1,
        1.50e-1,
        1.42e-1,
        1.35e-1,
        1.27e-1
    };

    // Units are nanometers
    const float AerosolWavelengths[] = {
        300.0,
        320.0,
        340.0,
        360.0,
        380.0,
        400.0,
        450.0,
        500.0,
        550.0,
        600.0,
        650.0,
        700.0,
        800.0
    };
#endif