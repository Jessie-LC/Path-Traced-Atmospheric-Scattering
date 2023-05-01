#if CAMERA_RESPONSE == 0
#include "responseCurveAGFAColorFuture100CD.glsl"
#elif CAMERA_RESPONSE == 1
#include "responseCurveDSCS315.glsl"
#elif CAMERA_RESPONSE == 2
#include "responseCurveFP2900ZG.glsl"
#elif CAMERA_RESPONSE == 3
#include "responseCurveFCICD.glsl"
#elif CAMERA_RESPONSE == 4
#include "responseCurveAGFAColorOptima2.glsl"
#endif