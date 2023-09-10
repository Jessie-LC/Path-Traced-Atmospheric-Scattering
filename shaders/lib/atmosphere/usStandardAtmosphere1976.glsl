#if !defined LIB_ATMOSPHERE_USSTANDARD1976
#define LIB_ATMOSPHERE_USSTANDARD1976
    const float g = 9.80665;
    /*
        Read this paper to understand the following model
        http://jimhawley.ca/downloads/Ballistics/Formulae_and_code_US_Standard_Atmosphere_1976.pdf
    */
    struct usStandardAtmosphere1976Data {
        float altAtBotOfLayer[7];
        float altAtTopOfLayer[7];
        float tempAtBotOfLayer[7];
        float pressureAtBotOfLayer[7];

        float lapseRateInLayer[7];
        float lapseExponentInLayer[7];
    } us1976Data;

    struct usStandardAtmosphere1976LU {
        float layerGeopotentialAlt;
        float layerGeometricAlt;

        float layerLocalTemp;
        float layerLocalPressure;
        float layerLocalDensity;

        int layer;
    };

    void fill1976LayerData(in float altitude) {
        us1976Data.tempAtBotOfLayer[0] = 15 + 273.15;
        us1976Data.pressureAtBotOfLayer[0] = 101325.0;

        us1976Data.altAtBotOfLayer[0] = 0.0;
        us1976Data.altAtBotOfLayer[1] = 11e3;
        us1976Data.altAtBotOfLayer[2] = 20e3;
        us1976Data.altAtBotOfLayer[3] = 32e3;
        us1976Data.altAtBotOfLayer[4] = 47e3;
        us1976Data.altAtBotOfLayer[5] = 51e3;
        us1976Data.altAtBotOfLayer[6] = 71e3;

        us1976Data.altAtTopOfLayer[0] = 11e3;
        us1976Data.altAtTopOfLayer[1] = 20e3;
        us1976Data.altAtTopOfLayer[2] = 32e3;
        us1976Data.altAtTopOfLayer[3] = 47e3;
        us1976Data.altAtTopOfLayer[4] = 51e3;
        us1976Data.altAtTopOfLayer[5] = 71e3;
        us1976Data.altAtTopOfLayer[6] = 85e3;

        us1976Data.lapseRateInLayer[0] = -6.5 / 1000.0;
        us1976Data.lapseRateInLayer[1] =           0.0;
        us1976Data.lapseRateInLayer[2] =  1.0 / 1000.0;
        us1976Data.lapseRateInLayer[3] =  2.8 / 1000.0;
        us1976Data.lapseRateInLayer[4] =           0.0;
        us1976Data.lapseRateInLayer[5] = -2.8 / 1000.0;
        us1976Data.lapseRateInLayer[6] = -2.0 / 1000.0;

        float geopotentialHeight = planetRadius * (1.0 - (planetRadius / (planetRadius + altitude)));

        for(int i = 0; i <= 6; ++i) {
            if(us1976Data.lapseRateInLayer[i] == 0.0) {
                us1976Data.lapseExponentInLayer[i] = -g / (rAir * us1976Data.tempAtBotOfLayer[i]);
            } else {
                us1976Data.lapseExponentInLayer[i] = -g / (rAir * us1976Data.lapseRateInLayer[i]);
            }

            if(i <= 5) {
                float thicknessOfLayer = min(us1976Data.altAtTopOfLayer[i], geopotentialHeight) - us1976Data.altAtBotOfLayer[i];
                us1976Data.tempAtBotOfLayer[i + 1] = us1976Data.tempAtBotOfLayer[i] + (thicknessOfLayer * us1976Data.lapseRateInLayer[i]);
            }
        }
        for(int i = 0; i <= 5; ++i) {
            float thicknessOfLayer = min(us1976Data.altAtTopOfLayer[i], geopotentialHeight) - us1976Data.altAtBotOfLayer[i];
            float temp = 0.0;
            if(us1976Data.lapseRateInLayer[i] == 0.0) {
                temp = us1976Data.lapseExponentInLayer[i] * thicknessOfLayer;
                us1976Data.pressureAtBotOfLayer[i + 1] = us1976Data.pressureAtBotOfLayer[i] * exp(temp);
            } else {
                temp = us1976Data.tempAtBotOfLayer[i + 1] / us1976Data.tempAtBotOfLayer[i];
                us1976Data.pressureAtBotOfLayer[i + 1] = us1976Data.pressureAtBotOfLayer[i] * pow(temp, us1976Data.lapseExponentInLayer[i]);
            }
        }
    }

    usStandardAtmosphere1976LU usStandardAtmosphere1976Lookup(in float altitude) {
        usStandardAtmosphere1976LU us1976LU;

        fill1976LayerData(altitude);

        us1976LU.layerGeometricAlt = altitude;
        us1976LU.layerGeopotentialAlt = planetRadius * (1.0 - (planetRadius / (planetRadius + us1976LU.layerGeometricAlt)));

        for(int i = 0; i <= 6; ++i) {
            if(us1976LU.layerGeopotentialAlt <= us1976Data.altAtTopOfLayer[i]) {
                us1976LU.layer = i;
                break;
            }
        }

        float differenceInAlt = us1976LU.layerGeopotentialAlt - us1976Data.altAtBotOfLayer[us1976LU.layer];

        us1976LU.layerLocalTemp = us1976Data.tempAtBotOfLayer[us1976LU.layer] + (differenceInAlt * us1976Data.lapseRateInLayer[us1976LU.layer]);

        float temp = 0.0;
        if(us1976Data.lapseRateInLayer[us1976LU.layer] == 0.0) {
            temp = us1976Data.lapseExponentInLayer[us1976LU.layer] * differenceInAlt;
            us1976LU.layerLocalPressure = us1976Data.pressureAtBotOfLayer[us1976LU.layer] * exp(temp);
        } else {
            temp = us1976LU.layerLocalTemp / us1976Data.tempAtBotOfLayer[us1976LU.layer];
            us1976LU.layerLocalPressure = us1976Data.pressureAtBotOfLayer[us1976LU.layer] * pow(temp, us1976Data.lapseExponentInLayer[us1976LU.layer]);
        }

        us1976LU.layerLocalDensity = us1976LU.layerLocalPressure / (rAir * us1976LU.layerLocalTemp);

        return us1976LU;
    }
#endif