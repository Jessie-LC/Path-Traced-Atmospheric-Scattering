#if !defined LIB_PCG_RNG
#define LIB_PCG_RNG
    // http://www.jcgt.org/published/0009/03/02/
    uvec3 pcg3d(uvec3 v) {
        v = v*1664525u + 1013904223u;
        v.x += v.y*v.z; v.y += v.z*v.x; v.z += v.x*v.y;
        v ^= v >> 16u;
        v.x += v.y*v.z; v.y += v.z*v.x; v.z += v.x*v.y;
        return v;
    }
    uvec4 pcg4d(uvec4 v) {
        v = v*1664525u + 1013904223u;
        v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
        v ^= v >> 16u;
        v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
        return v;
    }

    // https://nullprogram.com/blog/2018/07/31/
    uint lowbias32(uint x) {
        x ^= x >> 16;
        x *= 0x7feb352du;
        x ^= x >> 15;
        x *= 0x846ca68bu;
        x ^= x >> 16;
        return x;
    }

    uint randState;
    void InitRand(uint seed) { randState = lowbias32(seed); }
    uint RandNext() { return randState = lowbias32(randState); }
    #define RandNext2() uvec2(RandNext(), RandNext())
    #define RandNext3() uvec3(RandNext2(), RandNext())
    #define RandNext4() uvec4(RandNext3(), RandNext())
    #define RandNextF() (float(RandNext()) / float(0xffffffffu))
    #define RandNext2F() (vec2(RandNext2()) / float(0xffffffffu))
    #define RandNext3F() (vec3(RandNext3()) / float(0xffffffffu))
    #define RandNext4F() (vec4(RandNext4()) / float(0xffffffffu))
#endif