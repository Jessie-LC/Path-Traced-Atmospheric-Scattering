#define _square(x) (x*x)
#define _cube(x) (x*x*x)

#define _saturate(x) clamp(x, 0.0, 1.0)
#define _saturateInt(x) clamp(x, 0, 1)

#define _rcp(x) (1.0 / x)

#define _log10(x, y) (log2(x) / log2(y))

#define _pow5(x) (x*x*x*x*x)

#define landmask(x) (x < 1.0)

#define ConeAngleToSolidAngle(x) (tau*(1.0 - cos(x)))