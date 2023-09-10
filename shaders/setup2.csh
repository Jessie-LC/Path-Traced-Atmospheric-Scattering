#version 450

layout (local_size_x = 16, local_size_y = 16) in;
const ivec3 workGroups = ivec3(32, 32, 1);

layout (rgba32f) uniform image2D usStandardAtmosphere_image;

uniform sampler3D noise3D;

#include "/lib/universal/universal.glsl"
#include "/lib/rng/pcg.glsl"

#include "/lib/atmosphere/constants.glsl"
#include "/lib/atmosphere/usStandardAtmosphere1976.glsl"

void US_StandardAtmosphereLookupUVReverse(vec2 coord, out float R) {
	float uvR  = RemoveUvMargin(coord.y, 512);

	const float H = sqrt(atmosphereRadiusSquared - atmosphereLowerLimitSquared);

	float rho = H * uvR;
	R = sqrt(rho * rho + atmosphereLowerLimitSquared);
}

void main() {
    ivec2 fragPos = ivec2(gl_GlobalInvocationID.xy);

    float R;
    US_StandardAtmosphereLookupUVReverse(vec2(fragPos) / 512.0, R);

    usStandardAtmosphere1976LU us1976LU = usStandardAtmosphere1976Lookup(R - planetRadius);
    usStandardAtmosphere1976LU seaLevel = usStandardAtmosphere1976Lookup(0.0);

    if((R - planetRadius) > 85e3) {
        us1976LU.layerLocalDensity = exp(-R * inverseScaleHeights.x + scaledPlanetRadius.x) / 3.0;
    }

    imageStore(usStandardAtmosphere_image, fragPos, vec4(us1976LU.layerLocalDensity / seaLevel.layerLocalDensity, max(us1976LU.layerLocalPressure, 0.0), max(us1976LU.layerLocalTemp, 0.0), 0.0));
}