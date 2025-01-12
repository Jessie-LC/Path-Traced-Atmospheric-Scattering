#version 450

/*
    Copyright (C) 2023-2025  Jessie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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