#version 450

#extension GL_ARB_gpu_shader_int64 : enable
#extension AMD_gpu_shader_int64 : enable

#define CAMERA_RESPONSE 4 //[0 1 2 3 4]

layout (local_size_x = 32, local_size_y = 32) in;
const ivec3 workGroups = ivec3(32, 32, 1);

layout (rgba32f) uniform image2D cameraResponseLUT_image;

uniform sampler2D cameraResponseLUT;

#include "/lib/tonemap/response.glsl"

void main() {
    ivec2 fragPos = ivec2(gl_GlobalInvocationID.xy);

    float intensityRLUT = texelFetch(cameraResponseLUT, fragPos, 0).r + intensityColorR[fragPos.y];
    float intensityGLUT = texelFetch(cameraResponseLUT, fragPos, 0).g + intensityColorG[fragPos.y];
    float intensityBLUT = texelFetch(cameraResponseLUT, fragPos, 0).b + intensityColorB[fragPos.y];

    imageStore(cameraResponseLUT_image, fragPos, vec4(intensityRLUT, intensityGLUT, intensityBLUT, 0.0));
}