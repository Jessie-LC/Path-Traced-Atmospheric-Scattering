#if !defined RNG_MSWS
#define RNG_MSWS
    uint64_t x, w1, s;

    uint msws() {
        x *= x; 
        x += (w1 += s); 
        return uint(x = (x >> 32u) | (x << 32u));
    }

    void InitMSWS(uint64_t seed) {
        x = 0u; w1 = 0u;
        s = (((uint64_t(1890726812u) << 32u) | seed) << 1u) | uint64_t(1u);

        msws(); msws();
    }

    #define RandNext_MSWS() msws()
    #define RandNext2_MSWS() uvec2(msws(), msws())
    #define RandNext3_MSWS() uvec3(RandNext2(), msws())
    #define RandNext4_MSWS() uvec4(RandNext3(), msws())

    #define RandNextF_MSWS() (float(RandNext_MSWS() & 0x00ffffffu) / float(0x00ffffff))
    #define RandNextF2_MSWS() square(float(RandNext_MSWS() & 0x00ffffffu) / float(0x00ffffff))
    #define RandNext2F_MSWS() (vec2(RandNext2_MSWS() & 0x00ffffffu) / float(0x00ffffff))
    #define RandNext3F_MSWS() (vec3(RandNext3_MSWS() & 0x00ffffffu) / float(0x00ffffff))
#endif