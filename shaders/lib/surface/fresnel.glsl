#if !defined SURFACE_FRESNEL
#define SURFACE_FRESNEL
    complexFloat SnellsLaw(complexFloat theta_i, complexFloat n_i, complexFloat n_t) {
        return complexArcsin(complexMul(complexDiv(n_i, n_t), complexSin(theta_i)));
    }

    // complex wave amplitudes for reflection & transmission, for s & p polarizations
    complexFloat ReflectionSenkrecht(complexFloat theta_i, complexFloat theta_t, complexFloat n_i, complexFloat n_t) {
        return complexDiv(
            complexSub(complexMul(n_i, complexCos(theta_i)), complexMul(n_t, complexCos(theta_t))),
            complexAdd(complexMul(n_i, complexCos(theta_i)), complexMul(n_t, complexCos(theta_t)))
        );
    }
    complexFloat TransmittanceSenkrecht(complexFloat theta_i, complexFloat theta_t, complexFloat n_i, complexFloat n_t) {
        return complexDiv(
            complexMul(2.0, complexMul(n_i, complexCos(theta_i))),
            complexAdd(complexMul(n_i, complexCos(theta_i)), complexMul(n_t, complexCos(theta_t)))
        );
    }
    complexFloat ReflectionParallel(complexFloat theta_i, complexFloat theta_t, complexFloat n_i, complexFloat n_t) {
        return complexDiv(
            complexSub(complexMul(n_t, complexCos(theta_i)), complexMul(n_i, complexCos(theta_t))),
            complexAdd(complexMul(n_t, complexCos(theta_i)), complexMul(n_i, complexCos(theta_t)))
        );
    }
    complexFloat TransmittanceParallel(complexFloat theta_i, complexFloat theta_t, complexFloat n_i, complexFloat n_t) {
        return complexDiv(
            complexMul(2.0, complexMul(n_i, complexCos(theta_i))),
            complexAdd(complexMul(n_t, complexCos(theta_i)), complexMul(n_i, complexCos(theta_t)))
        );
    }

    float FresnelNonPolarized_R(in float cosTheta, in complexFloat n1, in complexFloat n2) {
        complexFloat thetaI = complexFloat(acos(cosTheta), 0.0);
        complexFloat thetaT = SnellsLaw(thetaI, n1, n2);

        if(complexAbs(complexSin(thetaT)) > 1.0) {
            return 1.0;
        }

        float Rs = square(complexAbs(ReflectionSenkrecht(thetaI, thetaT, n1, n2)));
        float Rp = square(complexAbs(ReflectionParallel(thetaI, thetaT, n1, n2)));

        return saturate((Rs + Rp) * 0.5);
    }

    float FresnelNonPolarized_T(in float cosTheta, in complexFloat n1, in complexFloat n2) {
        complexFloat thetaI = complexFloat(acos(cosTheta), 0.0);
        complexFloat thetaT = SnellsLaw(thetaI, n1, n2);

        if(complexAbs(complexSin(thetaT)) > 1.0) {
            return 0.0;
        }

        float Ts = square(complexAbs(TransmittanceSenkrecht(thetaI, thetaT, n1, n2)));
        float Tp = square(complexAbs(TransmittanceParallel(thetaI, thetaT, n1, n2)));

        complexFloat cosThetaI = complexCos(thetaI);
        complexFloat cosThetaT = complexCos(thetaT);

        float beamRatio = complexAbs(complexDiv(complexMul(n2, cosThetaT), complexMul(n1, cosThetaI)));

        return saturate(beamRatio * (Ts + Tp) * 0.5);
    }

    complexFloat PhaseChange(in complexFloat wavelength, in float thickness, in complexFloat theta_2) {
        //Both of these are equivalent, but I have both of them so I can use one or the other if I ever feel like it.
        //*
        complexFloat pathLengthWavelengths = complexDiv(complexMul(complexFloat(thickness, 0.0), complexCos(theta_2)), wavelength);
        complexFloat pathLengthRadians = complexMul(tau, pathLengthWavelengths);
        complexFloat phaseChange = complexExp(complexMul(complexFloat(0.0, 1.0), pathLengthRadians));
        return phaseChange;
        //*/
        /*
        complexFloat v = complexDiv(complexFloat(1.0, 0.0), wavelength);
        complexFloat D = complexMul(thickness, complexCos(theta_2));
        return complexExp(complexMul(complexFloat(0.0, 1.0), complexMul(tau, complexMul(v, D))));
        //*/
    }

    void ThinFilmAmplitudes(
        float incident_angle,
        float incident_wavelength,
        float thickness,
        complexFloat n_1, // media on current side
        complexFloat n_2, // media of thin film
        complexFloat n_3, // media on other side
        out complexFloat apparentAmplitudeSReflected,
        out complexFloat apparentAmplitudePReflected,
        out complexFloat apparentAmplitudeSTransmitted,
        out complexFloat apparentAmplitudePTransmitted
    ) {
        // Wave propagation direction changes in the different media according to Snell's law
        complexFloat theta_1 = complexFloat(incident_angle, 0.0);
        complexFloat theta_2 = SnellsLaw(theta_1, n_1, n_2);
        complexFloat theta_3 = SnellsLaw(theta_2, n_2, n_3);

        // Wavelength also changes in the different media in a similar manner
        complexFloat wavelength_1 = complexFloat(incident_wavelength / Air(incident_wavelength * 1e-3), 0.0);
        complexFloat wavelength_2 = complexMul(complexDiv(n_1, n_2), wavelength_1);
        complexFloat wavelength_3 = complexMul(complexDiv(n_1, n_3), wavelength_1);

        // compute wave amplitudes
        complexFloat propagationPhaseChange = PhaseChange(wavelength_2, thickness, theta_2);

        complexFloat amplitudeSIncident_1_2    = complexFloat(1.0, 0.0);
        complexFloat amplitudeSReflected_1_2   = complexMul(amplitudeSIncident_1_2, ReflectionSenkrecht(theta_1, theta_2, n_1, n_2));
        complexFloat amplitudeSTransmitted_1_2 = complexMul(amplitudeSIncident_1_2, TransmittanceSenkrecht(theta_1, theta_2, n_1, n_2));
        complexFloat amplitudePIncident_1_2    = complexFloat(1.0, 0.0);
        complexFloat amplitudePReflected_1_2   = complexMul(amplitudePIncident_1_2, ReflectionParallel(theta_1, theta_2, n_1, n_2));
        complexFloat amplitudePTransmitted_1_2 = complexMul(amplitudePIncident_1_2, TransmittanceParallel(theta_1, theta_2, n_1, n_2));

        complexFloat r_s_2_3 = ReflectionSenkrecht(theta_2, theta_3, n_2, n_3);
        complexFloat t_s_2_3 = TransmittanceSenkrecht(theta_2, theta_3, n_2, n_3);
        complexFloat r_p_2_3 = ReflectionParallel(theta_2, theta_3, n_2, n_3);
        complexFloat t_p_2_3 = TransmittanceParallel(theta_2, theta_3, n_2, n_3);

        complexFloat amplitudeSIncident_2_3    = complexMul(amplitudeSTransmitted_1_2, propagationPhaseChange);
        complexFloat amplitudeSReflected_2_3   = complexMul(amplitudeSIncident_2_3, r_s_2_3);
        complexFloat amplitudeSTransmitted_2_3 = complexMul(amplitudeSIncident_2_3, t_s_2_3);
        complexFloat amplitudePIncident_2_3    = complexMul(amplitudePTransmitted_1_2, propagationPhaseChange);
        complexFloat amplitudePReflected_2_3   = complexMul(amplitudePIncident_2_3, r_p_2_3);
        complexFloat amplitudePTransmitted_2_3 = complexMul(amplitudePIncident_2_3, t_p_2_3);

        complexFloat r_s_2_1 = ReflectionSenkrecht(theta_2, theta_1, n_2, n_1);
        complexFloat t_s_2_1 = TransmittanceSenkrecht(theta_2, theta_1, n_2, n_1);
        complexFloat r_p_2_1 = ReflectionParallel(theta_2, theta_1, n_2, n_1);
        complexFloat t_p_2_1 = TransmittanceParallel(theta_2, theta_1, n_2, n_1);

        apparentAmplitudeSReflected   = amplitudeSReflected_1_2;
        apparentAmplitudeSTransmitted = amplitudeSTransmitted_2_3;
        apparentAmplitudePReflected   = amplitudePReflected_1_2;
        apparentAmplitudePTransmitted = amplitudePTransmitted_2_3;

        complexFloat coeff_s_reflected   = complexMul(complexMul(r_s_2_1, r_s_2_3), complexMul(propagationPhaseChange, propagationPhaseChange));
        complexFloat coeff_p_reflected   = complexMul(complexMul(r_p_2_1, r_p_2_3), complexMul(propagationPhaseChange, propagationPhaseChange));
        complexFloat coeff_s_transmitted = complexMul(complexMul(r_s_2_3, r_s_2_1), complexMul(propagationPhaseChange, propagationPhaseChange));
        complexFloat coeff_p_transmitted = complexMul(complexMul(r_p_2_3, r_p_2_1), complexMul(propagationPhaseChange, propagationPhaseChange));

        complexFloat tmp_s_reflected = complexDiv(complexMul(amplitudeSReflected_2_3, propagationPhaseChange), complexSub(1.0, coeff_s_reflected));
        complexFloat tmp_p_reflected = complexDiv(complexMul(amplitudePReflected_2_3, propagationPhaseChange), complexSub(1.0, coeff_p_reflected));
        complexFloat tmp_s_transmitted = complexDiv(complexMul(complexMul(complexMul(amplitudeSReflected_2_3, propagationPhaseChange), r_s_2_1), propagationPhaseChange), complexSub(1.0, coeff_s_transmitted));
        complexFloat tmp_p_transmitted = complexDiv(complexMul(complexMul(complexMul(amplitudePReflected_2_3, propagationPhaseChange), r_p_2_1), propagationPhaseChange), complexSub(1.0, coeff_p_transmitted));

        apparentAmplitudeSReflected   = complexAdd(apparentAmplitudeSReflected,   complexMul(tmp_s_reflected,   t_s_2_1));
        apparentAmplitudePReflected   = complexAdd(apparentAmplitudePReflected,   complexMul(tmp_p_reflected,   t_p_2_1));
        apparentAmplitudeSTransmitted = complexAdd(apparentAmplitudeSTransmitted, complexMul(tmp_s_transmitted, t_s_2_3));
        apparentAmplitudePTransmitted = complexAdd(apparentAmplitudePTransmitted, complexMul(tmp_p_transmitted, t_p_2_3));
    }

    float FresnelThinFilmInterferenceReflected(in float cosTheta, in float thickness, in float lambda, in complexFloat n0, in complexFloat n1, in complexFloat n2) {
        complexFloat apparentAmplitudeSReflected;
        complexFloat apparentAmplitudePReflected;
        complexFloat apparentAmplitudeSTransmitted;
        complexFloat apparentAmplitudePTransmitted;
        ThinFilmAmplitudes(acos(cosTheta), lambda, thickness, n0, n1, n2, apparentAmplitudeSReflected, apparentAmplitudePReflected, apparentAmplitudeSTransmitted, apparentAmplitudePTransmitted);
        return saturate((square(complexAbs(apparentAmplitudeSReflected)) + square(complexAbs(apparentAmplitudePReflected))) / 2.0);
    }

    float FresnelThinFilmInterferenceTransmitted(in float cosTheta, in float thickness, in float lambda, in complexFloat n0, in complexFloat n1, in complexFloat n2) {
        complexFloat apparentAmplitudeSReflected;
        complexFloat apparentAmplitudePReflected;
        complexFloat apparentAmplitudeSTransmitted;
        complexFloat apparentAmplitudePTransmitted;
        ThinFilmAmplitudes(acos(cosTheta), lambda, thickness, n0, n1, n2, apparentAmplitudeSReflected, apparentAmplitudePReflected, apparentAmplitudeSTransmitted, apparentAmplitudePTransmitted);

        complexFloat theta_0 = complexFloat(acos(cosTheta), 0.0);
        complexFloat theta_1 = SnellsLaw(theta_0, n0, n1);
        complexFloat theta_2 = SnellsLaw(theta_0, n0, n2);

        if(complexAbs(complexSin(theta_1)) > 1.0 || complexAbs(complexSin(theta_2)) > 1.0) {
            return 0.0;
        }

        complexFloat cosThetaI = complexCos(theta_0);

        complexFloat cosThetaT_1 = complexCos(theta_1);
        complexFloat cosThetaT_2 = complexCos(theta_2);

        float beamRatio = complexAbs(complexDiv(complexMul(n2, cosThetaT_2), complexMul(n0, cosThetaI)));

        return saturate(beamRatio * (square(complexAbs(apparentAmplitudeSTransmitted)) + square(complexAbs(apparentAmplitudePTransmitted))) / 2.0);
    }
#endif