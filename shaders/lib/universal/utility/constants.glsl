const float pi = radians(180.0);
const float tau = radians(360.0);
const float rpi = 1./acos(-1.);
const float hpi = acos(0.);
const float pi4 = 12.5663706;
const float pidiv2 = pi * 0.5;
const float piSquared = pi*pi;

const vec3 lumacoeff_rec709 = vec3(0.2125, 0.7154, 0.0721);

const float ringRotation = 45.0;

const int noiseTextureResolution = 128;

const mat3 xyzToRGBMatrix_D50 = mat3(
     3.1338561, -1.6168667, -0.4906146,
    -0.9787684,  1.9161415,  0.0334540,
     0.0719453, -0.2289914,  1.4052427
);

const mat3 xyzToRGBMatrix_D65 = mat3(
	 3.2409699419, -1.5373831776, -0.4986107603,
	-0.9692436363,  1.8759675015,  0.0415550574,
	 0.0556300797, -0.2039769589,  1.0569715142
);

const mat3 rgbToXYZMatrix_D65 = mat3(
    0.4124564,  0.3575761,  0.1804375,
    0.2126729,  0.7151522,  0.0721750,
    0.0193339,  0.1191920,  0.9503041
);

const mat3 rgbToXYZMatrix_D50 = mat3(
    0.4360747,  0.3850649,  0.1430804,
    0.2225045,  0.7168786,  0.0606169,
    0.0139322,  0.0971045,  0.7141733
);

const mat3 xyzToDisplayP3Matrix = mat3(
     2.4933963, -0.9313459, -0.4026945,
    -0.8294868,  1.7626597,  0.0236246,
     0.0358507, -0.0761827,  0.9570140
);

const mat3 D50_2_D65 = mat3(
     0.9555766, -0.0230393,  0.0631636,
    -0.0282895,  1.0099416,  0.0210077,
     0.0122982, -0.0204830,  1.3299098
);

const int SpectraEntries = 8;

const float SpectraWavelengths[SpectraEntries] = float[SpectraEntries](
    6000.0,
    7800.0,
    8250.0,
    9000.0,
    9700.0,
    14500.0,
    17800.0,
    18000.0
);
const float SpectraValues[SpectraEntries] = float[SpectraEntries](
    0.0,
	0.0,
	0.87,
	0.7,
	1.0,
	0.25,
    0.0,
	0.0
);

const mat3 sRGB_D50_2_sRGB_D65 = (rgbToXYZMatrix_D50 * D50_2_D65) * xyzToRGBMatrix_D65;

const float preethamWavelengths[33] = float[](
    450.0,
    460.0,
    470.0,
    480.0,
    490.0,
    500.0,
    510.0,
    520.0,
    530.0,
    540.0,
    550.0,
    560.0,
    570.0,
    580.0,
    590.0,
    600.0,
    610.0,
    620.0,
    630.0,
    640.0,
    650.0,
    660.0,
    670.0,
    680.0,
    690.0,
    700.0,
    710.0,
    720.0,
    730.0,
    740.0,
    750.0,
    760.0,
    770.0
);

const float preethamOzone[33] = float[](
    0.3,
    0.6,
    0.9,
    1.4,
    2.1,
    3.0,
    4.0,
    4.8,
    6.3,
    7.5,
    8.5,
    10.3,
    12.0,
    12.0,
    11.5,
    12.5,
    12.0,
    10.5,
    9.0,
    7.9,
    6.7,
    5.7,
    4.8,
    3.6,
    2.8,
    2.3,
    1.8,
    1.4,
    1.1,
    1.0,
    0.9,
    0.7,
    0.4
);

const float ozoneCrossSection[441] = float[441](
    6.80778e-24,
    6.72106e-24,
    6.66971e-24,
    6.87827e-24,
    7.63950e-24,
    9.04948e-24,
    1.02622e-23,
    1.05505e-23,
    1.00303e-23,
    9.66106e-24,
    9.92189e-24,
    1.09556e-23,
    1.25580e-23,
    1.41026e-23,
    1.47637e-23,
    1.47593e-23,
    1.47278e-23,
    1.58624e-23,
    1.85887e-23,
    2.24518e-23,
    2.55393e-23,
    2.69799e-23,
    2.71700e-23,
    2.59416e-23,
    2.48635e-23,
    2.54214e-23,
    2.82383e-23,
    3.24288e-23,
    3.57929e-23,
    3.71260e-23,
    3.68606e-23,
    3.68453e-23,
    3.94202e-23,
    4.53838e-23,
    5.34144e-23,
    6.21470e-23,
    6.89342e-23,
    7.15131e-23,
    7.02279e-23,
    6.68451e-23,
    6.40985e-23,
    6.51215e-23,
    7.09484e-23,
    7.97582e-23,
    8.64796e-23,
    8.84751e-23,
    8.70576e-23,
    8.63018e-23,
    9.02596e-23,
    1.00978e-22,
    1.17515e-22,
    1.36812e-22,
    1.55564e-22,
    1.70593e-22,
    1.77413e-22,
    1.74113e-22,
    1.63969e-22,
    1.53825e-22,
    1.50061e-22,
    1.57091e-22,
    1.73006e-22,
    1.89872e-22,
    1.99180e-22,
    1.99402e-22,
    1.95992e-22,
    1.95885e-22,
    2.05664e-22,
    2.28073e-22,
    2.60481e-22,
    2.97688e-22,
    3.36293e-22,
    3.73195e-22,
    4.00920e-22,
    4.10428e-22,
    4.00191e-22,
    3.77591e-22,
    3.55367e-22,
    3.43550e-22,
    3.50045e-22,
    3.72443e-22,
    3.99079e-22,
    4.17388e-22,
    4.24576e-22,
    4.24739e-22,
    4.25983e-22,
    4.36706e-22,
    4.63007e-22,
    5.06215e-22,
    5.61756e-22,
    6.25345e-22,
    6.92671e-22,
    7.60101e-22,
    8.17582e-22,
    8.53087e-22,
    8.59583e-22,
    8.41161e-22,
    8.09704e-22,
    7.77762e-22,
    7.58661e-22,
    7.61105e-22,
    7.82768e-22,
    8.13525e-22,
    8.41416e-22,
    8.60281e-22,
    8.69574e-22,
    8.77739e-22,
    8.90289e-22,
    9.18185e-22,
    9.63101e-22,
    1.02541e-21,
    1.10497e-21,
    1.19583e-21,
    1.29472e-21,
    1.39640e-21,
    1.49041e-21,
    1.57014e-21,
    1.62239e-21,
    1.64414e-21,
    1.63511e-21,
    1.60943e-21,
    1.57830e-21,
    1.55493e-21,
    1.54503e-21,
    1.55300e-21,
    1.57805e-21,
    1.61238e-21,
    1.64978e-21,
    1.68423e-21,
    1.71542e-21,
    1.74504e-21,
    1.77787e-21,
    1.81470e-21,
    1.86234e-21,
    1.92426e-21,
    1.99836e-21,
    2.08321e-21,
    2.17570e-21,
    2.27551e-21,
    2.37767e-21,
    2.48026e-21,
    2.57787e-21,
    2.66735e-21,
    2.74553e-21,
    2.80416e-21,
    2.84156e-21,
    2.86077e-21,
    2.86533e-21,
    2.85907e-21,
    2.85266e-21,
    2.86095e-21,
    2.87845e-21,
    2.92588e-21,
    2.97008e-21,
    3.02468e-21,
    3.08141e-21,
    3.13490e-21,
    3.18141e-21,
    3.22207e-21,
    3.26213e-21,
    3.29445e-21,
    3.32516e-21,
    3.35579e-21,
    3.38847e-21,
    3.41886e-21,
    3.45674e-21,
    3.50070e-21,
    3.55041e-21,
    3.61007e-21,
    3.67904e-21,
    3.76616e-21,
    3.85792e-21,
    3.95625e-21,
    4.05115e-21,
    4.14698e-21,
    4.23308e-21,
    4.31213e-21,
    4.37493e-21,
    4.44152e-21,
    4.49554e-21,
    4.54212e-21,
    4.59922e-21,
    4.65627e-21,
    4.70549e-21,
    4.75188e-21,
    4.78362e-21,
    4.79933e-21,
    4.79812e-21,
    4.78287e-21,
    4.74991e-21,
    4.70931e-21,
    4.65747e-21,
    4.61692e-21,
    4.57024e-21,
    4.52700e-21,
    4.48817e-21,
    4.45931e-21,
    4.43458e-21,
    4.41148e-21,
    4.40927e-21,
    4.40508e-21,
    4.41249e-21,
    4.43419e-21,
    4.46445e-21,
    4.50560e-21,
    4.56963e-21,
    4.64735e-21,
    4.73301e-21,
    4.82020e-21,
    4.91050e-21,
    4.99163e-21,
    5.06017e-21,
    5.11838e-21,
    5.16436e-21,
    5.18613e-21,
    5.19008e-21,
    5.17248e-21,
    5.13839e-21,
    5.07516e-21,
    5.00213e-21,
    4.92632e-21,
    4.84196e-21,
    4.75813e-21,
    4.66949e-21,
    4.58682e-21,
    4.50504e-21,
    4.42659e-21,
    4.34938e-21,
    4.27621e-21,
    4.20827e-21,
    4.14570e-21,
    4.08986e-21,
    4.03221e-21,
    3.99139e-21,
    3.94294e-21,
    3.90294e-21,
    3.85486e-21,
    3.80352e-21,
    3.75269e-21,
    3.69724e-21,
    3.64581e-21,
    3.59756e-21,
    3.53981e-21,
    3.48189e-21,
    3.42639e-21,
    3.36507e-21,
    3.30716e-21,
    3.24798e-21,
    3.19212e-21,
    3.13235e-21,
    3.07385e-21,
    3.01187e-21,
    2.94933e-21,
    2.88675e-21,
    2.83154e-21,
    2.77990e-21,
    2.73430e-21,
    2.69151e-21,
    2.64926e-21,
    2.60694e-21,
    2.56838e-21,
    2.52929e-21,
    2.49407e-21,
    2.45557e-21,
    2.41588e-21,
    2.37737e-21,
    2.33497e-21,
    2.29460e-21,
    2.25198e-21,
    2.21134e-21,
    2.16653e-21,
    2.12952e-21,
    2.09231e-21,
    2.05119e-21,
    2.01199e-21,
    1.96873e-21,
    1.93030e-21,
    1.89301e-21,
    1.85458e-21,
    1.80984e-21,
    1.76722e-21,
    1.72459e-21,
    1.68500e-21,
    1.64647e-21,
    1.60911e-21,
    1.57194e-21,
    1.53783e-21,
    1.50400e-21,
    1.47295e-21,
    1.44342e-21,
    1.41512e-21,
    1.38809e-21,
    1.36429e-21,
    1.34049e-21,
    1.31934e-21,
    1.30100e-21,
    1.28154e-21,
    1.26035e-21,
    1.23594e-21,
    1.20922e-21,
    1.18024e-21,
    1.14995e-21,
    1.11892e-21,
    1.09140e-21,
    1.06392e-21,
    1.03712e-21,
    1.01065e-21,
    9.84534e-22,
    9.58011e-22,
    9.31230e-22,
    9.06905e-22,
    8.83424e-22,
    8.61809e-22,
    8.41371e-22,
    8.23199e-22,
    8.07479e-22,
    7.92359e-22,
    7.78960e-22,
    7.66918e-22,
    7.56724e-22,
    7.45938e-22,
    7.36321e-22,
    7.26761e-22,
    7.17708e-22,
    7.10170e-22,
    7.04603e-22,
    7.00521e-22,
    6.95807e-22,
    6.87983e-22,
    6.75690e-22,
    6.59167e-22,
    6.38658e-22,
    6.17401e-22,
    5.97986e-22,
    5.79980e-22,
    5.64879e-22,
    5.52304e-22,
    5.40930e-22,
    5.28950e-22,
    5.14905e-22,
    5.00676e-22,
    4.86900e-22,
    4.74324e-22,
    4.63744e-22,
    4.54117e-22,
    4.47413e-22,
    4.42084e-22,
    4.38598e-22,
    4.35751e-22,
    4.32496e-22,
    4.30002e-22,
    4.28472e-22,
    4.27365e-22,
    4.29043e-22,
    4.31385e-22,
    4.35345e-22,
    4.40512e-22,
    4.46268e-22,
    4.50925e-22,
    4.51983e-22,
    4.49671e-22,
    4.41359e-22,
    4.27561e-22,
    4.09127e-22,
    3.88901e-22,
    3.68851e-22,
    3.50462e-22,
    3.34368e-22,
    3.20386e-22,
    3.08569e-22,
    2.99026e-22,
    2.90708e-22,
    2.83838e-22,
    2.77892e-22,
    2.72682e-22,
    2.67864e-22,
    2.63381e-22,
    2.60147e-22,
    2.57597e-22,
    2.55903e-22,
    2.54995e-22,
    2.55263e-22,
    2.56910e-22,
    2.59848e-22,
    2.64943e-22,
    2.72251e-22,
    2.81519e-22,
    2.92565e-22,
    3.03612e-22,
    3.13434e-22,
    3.20710e-22,
    3.23925e-22,
    3.21425e-22,
    3.14522e-22,
    3.03211e-22,
    2.89017e-22,
    2.73981e-22,
    2.59406e-22,
    2.46085e-22,
    2.34234e-22,
    2.23936e-22,
    2.14826e-22,
    2.06425e-22,
    1.98427e-22,
    1.90789e-22,
    1.83692e-22,
    1.77111e-22,
    1.71523e-22,
    1.66604e-22,
    1.63367e-22,
    1.60371e-22,
    1.57834e-22,
    1.55432e-22,
    1.53961e-22,
    1.52632e-22,
    1.51695e-22,
    1.51650e-22,
    1.53341e-22,
    1.56550e-22,
    1.61557e-22,
    1.68582e-22,
    1.76205e-22,
    1.84627e-22,
    1.93246e-22,
    2.01741e-22,
    2.09583e-22,
    2.16778e-22,
    2.22566e-22,
    2.25770e-22,
    2.25611e-22,
    2.22491e-22,
    2.16317e-22,
    2.07365e-22,
    1.96101e-22,
    1.82575e-22,
    1.69093e-22,
    1.55152e-22,
    1.42655e-22,
    1.31245e-22,
    1.21519e-22,
    1.12924e-22,
    1.05472e-22
);

const float cloudAbsorption[441] = {
    7.44369 / 100.0,
    7.44426 / 100.0,
    7.44237 / 100.0,
    7.44366 / 100.0,
    7.43659 / 100.0,
    7.44048 / 100.0,
    7.44024 / 100.0,
    7.43617 / 100.0,
    7.43727 / 100.0,
    7.438 / 100.0,
    7.43854 / 100.0,
    7.43895 / 100.0,
    7.4397 / 100.0,
    7.43496 / 100.0,
    7.43371 / 100.0,
    7.42992 / 100.0,
    7.4302 / 100.0,
    7.42976 / 100.0,
    7.42477 / 100.0,
    7.42714 / 100.0,
    7.43138 / 100.0,
    7.42491 / 100.0,
    7.42602 / 100.0,
    7.42303 / 100.0,
    7.42749 / 100.0,
    7.41861 / 100.0,
    7.42385 / 100.0,
    7.42102 / 100.0,
    7.41913 / 100.0,
    7.42937 / 100.0,
    7.41909 / 100.0,
    7.41842 / 100.0,
    7.41615 / 100.0,
    7.41964 / 100.0,
    7.42028 / 100.0,
    7.41438 / 100.0,
    7.4127 / 100.0,
    7.40987 / 100.0,
    7.40793 / 100.0,
    7.41323 / 100.0,
    7.41197 / 100.0,
    7.40603 / 100.0,
    7.41231 / 100.0,
    7.41469 / 100.0,
    7.41554 / 100.0,
    7.41148 / 100.0,
    7.40935 / 100.0,
    7.40483 / 100.0,
    7.4056 / 100.0,
    7.41553 / 100.0,
    7.40949 / 100.0,
    7.40617 / 100.0,
    7.40402 / 100.0,
    7.41178 / 100.0,
    7.40851 / 100.0,
    7.39922 / 100.0,
    7.41048 / 100.0,
    7.40713 / 100.0,
    7.40338 / 100.0,
    7.40337 / 100.0,
    7.40289 / 100.0,
    7.4053 / 100.0,
    7.40893 / 100.0,
    7.40206 / 100.0,
    7.40191 / 100.0,
    7.39898 / 100.0,
    7.40057 / 100.0,
    7.40244 / 100.0,
    7.40001 / 100.0,
    7.39979 / 100.0,
    7.42757 / 100.0,
    7.37902 / 100.0,
    7.40995 / 100.0,
    7.4009 / 100.0,
    7.39471 / 100.0,
    7.39483 / 100.0,
    7.39576 / 100.0,
    7.40352 / 100.0,
    7.40168 / 100.0,
    7.39522 / 100.0,
    7.40305 / 100.0,
    7.39564 / 100.0,
    7.40137 / 100.0,
    7.39787 / 100.0,
    7.3976 / 100.0,
    7.40342 / 100.0,
    7.39832 / 100.0,
    7.4012 / 100.0,
    7.39519 / 100.0,
    7.3978 / 100.0,
    7.39494 / 100.0,
    7.40038 / 100.0,
    7.39864 / 100.0,
    7.39472 / 100.0,
    7.39884 / 100.0,
    7.38913 / 100.0,
    7.39315 / 100.0,
    7.39335 / 100.0,
    7.39168 / 100.0,
    7.39241 / 100.0,
    7.40136 / 100.0,
    7.39309 / 100.0,
    7.39771 / 100.0,
    7.40095 / 100.0,
    7.39643 / 100.0,
    7.39608 / 100.0,
    7.39205 / 100.0,
    7.3932 / 100.0,
    7.39594 / 100.0,
    7.39562 / 100.0,
    7.39805 / 100.0,
    7.39671 / 100.0,
    7.39611 / 100.0,
    7.39008 / 100.0,
    7.4027 / 100.0,
    7.38009 / 100.0,
    7.38824 / 100.0,
    7.40125 / 100.0,
    7.39234 / 100.0,
    7.39779 / 100.0,
    7.40083 / 100.0,
    7.39202 / 100.0,
    7.40509 / 100.0,
    7.40316 / 100.0,
    7.35708 / 100.0,
    7.3944 / 100.0,
    7.41114 / 100.0,
    7.39429 / 100.0,
    7.39545 / 100.0,
    7.3836 / 100.0,
    7.395 / 100.0,
    7.39458 / 100.0,
    7.38669 / 100.0,
    7.38021 / 100.0,
    7.39716 / 100.0,
    7.40051 / 100.0,
    7.38904 / 100.0,
    7.39028 / 100.0,
    7.38837 / 100.0,
    7.39004 / 100.0,
    7.3951 / 100.0,
    7.39345 / 100.0,
    7.39258 / 100.0,
    7.38934 / 100.0,
    7.39572 / 100.0,
    7.39238 / 100.0,
    7.39593 / 100.0,
    7.4001 / 100.0,
    7.39993 / 100.0,
    7.39302 / 100.0,
    7.39317 / 100.0,
    7.39492 / 100.0,
    7.39096 / 100.0,
    7.38996 / 100.0,
    7.38615 / 100.0,
    7.3936 / 100.0,
    7.39401 / 100.0,
    7.39131 / 100.0,
    7.39627 / 100.0,
    7.39108 / 100.0,
    7.39523 / 100.0,
    7.39286 / 100.0,
    7.40267 / 100.0,
    7.39808 / 100.0,
    7.38849 / 100.0,
    7.396 / 100.0,
    7.39183 / 100.0,
    7.40558 / 100.0,
    7.3985 / 100.0,
    7.3936 / 100.0,
    7.39533 / 100.0,
    7.39579 / 100.0,
    7.39362 / 100.0,
    7.3895 / 100.0,
    7.3963 / 100.0,
    7.40148 / 100.0,
    7.39548 / 100.0,
    7.39092 / 100.0,
    7.4 / 100.0,
    7.39153 / 100.0,
    7.39496 / 100.0,
    7.39388 / 100.0,
    7.40492 / 100.0,
    7.40778 / 100.0,
    7.37187 / 100.0,
    7.43124 / 100.0,
    7.38747 / 100.0,
    7.39022 / 100.0,
    7.39501 / 100.0,
    7.40354 / 100.0,
    7.39401 / 100.0,
    7.39009 / 100.0,
    7.4205 / 100.0,
    7.39296 / 100.0,
    7.40742 / 100.0,
    7.41192 / 100.0,
    7.403 / 100.0,
    7.40521 / 100.0,
    7.41501 / 100.0,
    7.40956 / 100.0,
    7.40306 / 100.0,
    7.40933 / 100.0,
    7.40741 / 100.0,
    7.40444 / 100.0,
    7.40599 / 100.0,
    7.40644 / 100.0,
    7.4097 / 100.0,
    7.412 / 100.0,
    7.41288 / 100.0,
    7.41278 / 100.0,
    7.41837 / 100.0,
    7.41316 / 100.0,
    7.40908 / 100.0,
    7.41697 / 100.0,
    7.42829 / 100.0,
    7.41044 / 100.0,
    7.406 / 100.0,
    7.41235 / 100.0,
    7.41266 / 100.0,
    7.42464 / 100.0,
    7.41415 / 100.0,
    7.41256 / 100.0,
    7.416 / 100.0,
    7.42248 / 100.0,
    7.41641 / 100.0,
    7.41091 / 100.0,
    7.41435 / 100.0,
    7.42029 / 100.0,
    7.41867 / 100.0,
    7.41802 / 100.0,
    7.42992 / 100.0,
    7.40679 / 100.0,
    7.4308 / 100.0,
    7.42164 / 100.0,
    7.41831 / 100.0,
    7.41508 / 100.0,
    7.41346 / 100.0,
    7.43364 / 100.0,
    7.41665 / 100.0,
    7.41976 / 100.0,
    7.41635 / 100.0,
    7.43723 / 100.0,
    7.4273 / 100.0,
    7.42517 / 100.0,
    7.42444 / 100.0,
    7.42259 / 100.0,
    7.42394 / 100.0,
    7.41983 / 100.0,
    7.42762 / 100.0,
    7.41976 / 100.0,
    7.43244 / 100.0,
    7.42106 / 100.0,
    7.42936 / 100.0,
    7.42495 / 100.0,
    7.44414 / 100.0,
    7.42078 / 100.0,
    7.44676 / 100.0,
    7.42098 / 100.0,
    7.44445 / 100.0,
    7.41864 / 100.0,
    7.42941 / 100.0,
    7.4251 / 100.0,
    7.43182 / 100.0,
    7.42954 / 100.0,
    7.4314 / 100.0,
    7.43721 / 100.0,
    7.44456 / 100.0,
    7.43109 / 100.0,
    7.42906 / 100.0,
    7.4298 / 100.0,
    7.4508 / 100.0,
    7.4342 / 100.0,
    7.43586 / 100.0,
    7.43306 / 100.0,
    7.43021 / 100.0,
    7.43121 / 100.0,
    7.44066 / 100.0,
    7.44295 / 100.0,
    7.44108 / 100.0,
    7.43997 / 100.0,
    7.43272 / 100.0,
    7.44028 / 100.0,
    7.43747 / 100.0,
    7.43427 / 100.0,
    7.44622 / 100.0,
    7.4391 / 100.0,
    7.4376 / 100.0,
    7.45194 / 100.0,
    7.44479 / 100.0,
    7.45874 / 100.0,
    7.45051 / 100.0,
    7.44496 / 100.0,
    7.53381 / 100.0,
    7.45718 / 100.0,
    7.45751 / 100.0,
    7.4896 / 100.0,
    7.45706 / 100.0,
    7.46034 / 100.0,
    7.45948 / 100.0,
    7.44555 / 100.0,
    7.45287 / 100.0,
    7.45549 / 100.0,
    7.45493 / 100.0,
    7.46184 / 100.0,
    7.46189 / 100.0,
    7.46408 / 100.0,
    7.45588 / 100.0,
    7.45753 / 100.0,
    7.46198 / 100.0,
    7.46271 / 100.0,
    7.46544 / 100.0,
    7.4632 / 100.0,
    7.47342 / 100.0,
    7.47482 / 100.0,
    7.48452 / 100.0,
    7.48221 / 100.0,
    7.48353 / 100.0,
    7.4997 / 100.0,
    7.49911 / 100.0,
    7.4993 / 100.0,
    7.5522 / 100.0,
    7.50267 / 100.0,
    7.51144 / 100.0,
    7.50578 / 100.0,
    7.53346 / 100.0,
    7.51552 / 100.0,
    7.52735 / 100.0,
    7.52743 / 100.0,
    7.57831 / 100.0,
    7.53197 / 100.0,
    7.52978 / 100.0,
    7.54173 / 100.0,
    7.56614 / 100.0,
    7.55705 / 100.0,
    7.59301 / 100.0,
    7.54522 / 100.0,
    7.56463 / 100.0,
    7.55138 / 100.0,
    7.55786 / 100.0,
    7.57221 / 100.0,
    7.60814 / 100.0,
    7.5882 / 100.0,
    7.57242 / 100.0,
    7.59993 / 100.0,
    7.6082 / 100.0,
    7.58127 / 100.0,
    7.5884 / 100.0,
    7.58572 / 100.0,
    7.59131 / 100.0,
    7.60518 / 100.0,
    7.61574 / 100.0,
    7.61563 / 100.0,
    7.61516 / 100.0,
    7.60157 / 100.0,
    7.60026 / 100.0,
    7.60652 / 100.0,
    7.67741 / 100.0,
    7.64687 / 100.0,
    7.67279 / 100.0,
    7.70159 / 100.0,
    7.72103 / 100.0,
    7.7161 / 100.0,
    7.68539 / 100.0,
    7.64107 / 100.0,
    7.57223 / 100.0,
    7.53758 / 100.0,
    7.53856 / 100.0,
    7.53314 / 100.0,
    7.56183 / 100.0,
    7.58084 / 100.0,
    7.62992 / 100.0,
    7.65618 / 100.0,
    7.68096 / 100.0,
    7.70527 / 100.0,
    7.73576 / 100.0,
    7.64485 / 100.0,
    7.60282 / 100.0,
    7.62212 / 100.0,
    7.77431 / 100.0,
    7.75376 / 100.0,
    7.67093 / 100.0,
    7.61626 / 100.0,
    7.57092 / 100.0,
    7.66357 / 100.0,
    7.58137 / 100.0,
    7.57466 / 100.0,
    7.65937 / 100.0,
    7.67859 / 100.0,
    7.60688 / 100.0,
    7.61323 / 100.0,
    8.10502 / 100.0,
    7.7119 / 100.0,
    7.65393 / 100.0,
    7.64227 / 100.0,
    7.61787 / 100.0,
    7.60777 / 100.0,
    7.60467 / 100.0,
    7.59035 / 100.0,
    7.58669 / 100.0,
    7.5884 / 100.0,
    7.7001 / 100.0,
    7.62653 / 100.0,
    7.61129 / 100.0,
    7.60498 / 100.0,
    7.60185 / 100.0,
    7.59896 / 100.0,
    7.59434 / 100.0,
    7.59264 / 100.0,
    7.92075 / 100.0,
    8.19364 / 100.0,
    7.6099 / 100.0,
    7.61504 / 100.0,
    7.60997 / 100.0,
    7.61078 / 100.0,
    7.62613 / 100.0,
    7.65093 / 100.0,
    7.6267 / 100.0,
    7.64118 / 100.0,
    7.62631 / 100.0,
    7.62729 / 100.0,
    7.63407 / 100.0,
    7.63392 / 100.0,
    7.64808 / 100.0,
    7.64236 / 100.0,
    7.65502 / 100.0,
    7.65061 / 100.0,
    7.6647 / 100.0,
    7.6494 / 100.0,
    7.66776 / 100.0,
    7.66824 / 100.0,
    7.66392 / 100.0,
    7.68011 / 100.0,
    7.68083 / 100.0,
    7.66624 / 100.0,
    7.67259 / 100.0,
    7.69495 / 100.0,
    7.68608 / 100.0,
    7.68834 / 100.0,
    7.69685 / 100.0,
    7.72064 / 100.0,
    7.71863 / 100.0
};