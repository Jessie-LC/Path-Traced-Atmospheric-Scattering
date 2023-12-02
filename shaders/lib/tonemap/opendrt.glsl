// Not sure if these should be natural or base 10
// Seems it'll only affect hlg though, which will be unused
#define _expf(x) exp(x)
#define _logf(x) log(x)

/*  OpenDRT -------------------------------------------------/
      v0.2.8
      Written by Jed Smith
      Ported to GLSL by Jacob Eriksson
      https://github.com/jedypod/open-display-transform

      License: GPL v3
-------------------------------------------------*/

#define Lp 100.0 // Lp [100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 750 760 770 780 790 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000]
#define gb 0.12 // Grey boost [0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0]

#define i_xyz 0
#define i_ap0 1
#define i_ap1 2
#define i_p3d65 3
#define i_rec2020 4
#define i_rec709 5
#define i_awg3 6
#define i_awg4 7
#define i_rwg 8
#define i_sgamut3 9
#define i_sgamut3cine 10
#define i_bmdwg 11
#define i_egamut 12
#define i_davinciwg 13
#define in_gamut i_rec709 // Input gamut [i_xyz i_ap0 i_ap1 i_p3d65 i_rec2020 i_rec709 i_awg3 i_awg4 i_rwg i_sgamut3 i_sgamut3cine i_bmdwg i_egamut i_davinciwg]
// XYZ, ACES 2065-1, ACEScg, P3D65, Rec.2020, Rec.709, Arri Wide Gamut 3, Arri Wide Gamut 4, Red Wide Gamut RGB, Sony SGamut3, Sony SGamut3Cine, Blackmagic Wide Gamut, Filmlight E - Gamut, DaVinci Wide Gamut

#define Rec709 0
#define P3D65 1
#define Rec2020 2
#define display_gamut Rec709 // Display gamut [Rec709 P3D65 Rec2020]
// Rec.709, P3 D65, Rec.2020

#define lin 0
#define srgb 1
#define rec1886 2
#define dci 3
#define pq 4
#define hlg 5
#define EOTF srgb // Display EOTF [lin srgb rec1886 dci pq hlg]
// Linear, 2.2 Power sRGB Display, 2.4 Power Rec .1886, 2.6 Power DCI, ST 2084 PQ, HLG

// Gamut Conversion Matrices
#define matrix_ap0_to_xyz mat3(vec3(0.93863094875f, -0.00574192055f, 0.017566898852f), vec3(0.338093594922f, 0.727213902811f, -0.065307497733f), vec3(0.000723121511f, 0.000818441849f, 1.0875161874f))
#define matrix_ap1_to_xyz mat3(vec3(0.652418717672f, 0.127179925538f, 0.170857283842f), vec3(0.268064059194f, 0.672464478993f, 0.059471461813f), vec3(-0.00546992851f, 0.005182799977f, 1.08934487929f))
#define matrix_rec709_to_xyz mat3(vec3(0.412390917540f, 0.357584357262f, 0.180480793118f), vec3(0.212639078498f, 0.715168714523f, 0.072192311287f), vec3(0.019330825657f, 0.119194783270f, 0.950532138348f))
#define matrix_p3d65_to_xyz mat3(vec3(0.486571133137f, 0.265667706728f, 0.198217317462f), vec3(0.228974640369f, 0.691738605499f, 0.079286918044f), vec3(-0.000000000000f, 0.045113388449, 1.043944478035f))
#define matrix_rec2020_to_xyz mat3(vec3(0.636958122253f, 0.144616916776f, 0.168880969286f), vec3(0.262700229883f, 0.677998125553f, 0.059301715344f), vec3(0.000000000000f, 0.028072696179, 1.060985088348f))
#define matrix_arriwg3_to_xyz mat3(vec3(0.638007619284f, 0.214703856337f, 0.097744451431f), vec3(0.291953779f, 0.823841041511f, -0.11579482051f), vec3(0.002798279032f, -0.067034235689f, 1.15329370742f))
#define matrix_arriwg4_to_xyz mat3(vec3(0.704858320407f, 0.12976029517f, 0.115837311474f), vec3(0.254524176404f, 0.781477732712f, -0.036001909116f), vec3(0.0f, 0.0f, 1.08905775076f))
#define matrix_redwg_to_xyz mat3(vec3(0.735275208950f, 0.068609409034f, 0.146571278572f), vec3(0.286694079638f, 0.842979073524f, -0.129673242569f), vec3(-0.079680845141f, -0.347343206406, 1.516081929207f))
#define matrix_sonysgamut3_to_xyz mat3(vec3(0.706482713192f, 0.128801049791f, 0.115172164069f), vec3(0.270979670813f, 0.786606411221f, -0.057586082034f), vec3(-0.009677845386f, 0.004600037493f, 1.09413555865f))
#define matrix_sonysgamut3cine_to_xyz mat3(vec3(0.599083920758f, 0.248925516115f, 0.102446490178f), vec3(0.215075820116f, 0.885068501744f, -0.100144321859f), vec3(-0.032065849545f, -0.027658390679f, 1.14878199098f))
#define matrix_bmdwg_to_xyz mat3(vec3(0.606538414955f, 0.220412746072f, 0.123504832387f), vec3(0.267992943525f, 0.832748472691f, -0.100741356611f), vec3(-0.029442556202f, -0.086612440646, 1.205112814903f))
#define matrix_egamut_to_xyz mat3(vec3(0.705396831036f, 0.164041340351f, 0.081017754972f), vec3(0.280130714178f, 0.820206701756f, -0.100337378681f), vec3(-0.103781513870f, -0.072907261550, 1.265746593475f))
#define matrix_davinciwg_to_xyz mat3(vec3(0.700622320175f, 0.148774802685f, 0.101058728993f), vec3(0.274118483067f, 0.873631775379f, -0.147750422359f), vec3(-0.098962903023f, -0.137895315886, 1.325916051865f))

#define matrix_xyz_to_rec709 mat3(vec3(3.2409699419f, -1.53738317757f, -0.498610760293f), vec3(-0.969243636281f, 1.87596750151f, 0.041555057407f), vec3(0.055630079697f, -0.203976958889f, 1.05697151424f))
#define matrix_xyz_to_p3d65 mat3(vec3(2.49349691194f, -0.931383617919f, -0.402710784451f), vec3(-0.829488969562f, 1.76266406032f, 0.023624685842f), vec3(0.035845830244f, -0.076172389268f, 0.956884524008f))
#define matrix_xyz_to_rec2020 mat3(vec3(1.71665118797f, -0.355670783776f, -0.253366281374f), vec3(-0.666684351832f, 1.61648123664f, 0.015768545814f), vec3(0.017639857445f, -0.042770613258f, 0.942103121235f))

/* Math helper functions ----------------------------*/

// Return identity 3x3 matrix
mat3 identity() {
  return mat3(vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f));
}

// Multiply 3x3 matrix m and vec3 vector v
vec3 vdot(mat3 m, vec3 v) {
  return vec3(m[0].x*v.x + m[0].y*v.y + m[0].z*v.z, m[1].x*v.x + m[1].y*v.y + m[1].z*v.z, m[2].x*v.x + m[2].y*v.y + m[2].z*v.z);
}

// Safe division of float a by float b
float sdivf(float a, float b) {
  if (b == 0.0f) return 0.0f;
  else return a/b;
}

// Safe division of vec3 a by float b
vec3 sdivf3f(vec3 a, float b) {
  return vec3(sdivf(a.x, b), sdivf(a.y, b), sdivf(a.z, b));
}

// Safe division of vec3 a by vec3 b
vec3 sdivf3f3(vec3 a, vec3 b) {
  return vec3(sdivf(a.x, b.x), sdivf(a.y, b.y), sdivf(a.z, b.z));
}

// Safe power function raising float a to power float b
float spowf(float a, float b) {
  if (a <= 0.0f) return a;
  else return pow(a, b);
}

// Safe power function raising vec3 a to power float b
vec3 spowf3(vec3 a, float b) {
  return vec3(pow(a.x, b), pow(a.y, b), pow(a.z, b));
}

// Return max of vec3 a and float mn
vec3 maxf3(float mn, vec3 a) { return vec3(max(a.x, mn), max(a.y, mn), max(a.z, mn)); }

// Return min of vec3 a and float mx
vec3 minf3(float mx, vec3 a) { return vec3(min(a.x, mx), min(a.y, mx), min(a.z, mx)); }

// Return the hypot or length of vec3 a
float hypotf3(vec3 a) { return  sqrt(spowf(a.x, 2.0f) + spowf(a.y, 2.0f) + spowf(a.z, 2.0f)); }

vec3 eotf_hlg(vec3 rgb, int inverse) {
  // Aply the HLG Forward or Inverse EOTF. Implements the full ambient surround illumination model
  // ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
  // ITU-R Rep BT.2390-8: https://www.itu.int/pub/R-REP-BT.2390
  // Perceptual Quantiser (PQ) to Hybrid Log-Gamma (HLG) Transcoding: https://www.bbc.co.uk/rd/sites/50335ff370b5c262af000004/assets/592eea8006d63e5e5200f90d/BBC_HDRTV_PQ_HLG_Transcode_v2.pdf

  const float HLG_Lw = 1000.0f;
  // const float HLG_Lb = 0.0f;
  const float HLG_Ls = 5.0f;
  const float h_a = 0.17883277f;
  const float h_b = 1.0f - 4.0f*0.17883277f;
  const float h_c = 0.5f - h_a*_logf(4.0f*h_a);
  const float h_g = 1.2f*pow(1.111f, log2(HLG_Lw/1000.0f))*pow(0.98f, log2(max(1e-6f, HLG_Ls)/5.0f));
  if (inverse == 1) {
    float Yd = 0.2627f*rgb.x + 0.6780f*rgb.y + 0.0593f*rgb.z;
    // HLG Inverse OOTF
    rgb = rgb*pow(Yd, (1.0f - h_g)/h_g);
    // HLG OETF
    rgb.x = rgb.x <= 1.0f/12.0f ? sqrt(3.0f*rgb.x) : h_a*_logf(12.0f*rgb.x - h_b) + h_c;
    rgb.y = rgb.y <= 1.0f/12.0f ? sqrt(3.0f*rgb.y) : h_a*_logf(12.0f*rgb.y - h_b) + h_c;
    rgb.z = rgb.z <= 1.0f/12.0f ? sqrt(3.0f*rgb.z) : h_a*_logf(12.0f*rgb.z - h_b) + h_c;
  } else {
    // HLG Inverse OETF
    rgb.x = rgb.x <= 0.5f ? rgb.x*rgb.x/3.0f : (_expf((rgb.x - h_c)/h_a) + h_b)/12.0f;
    rgb.y = rgb.y <= 0.5f ? rgb.y*rgb.y/3.0f : (_expf((rgb.y - h_c)/h_a) + h_b)/12.0f;
    rgb.z = rgb.z <= 0.5f ? rgb.z*rgb.z/3.0f : (_expf((rgb.z - h_c)/h_a) + h_b)/12.0f;
    // HLG OOTF
    float Ys = 0.2627f*rgb.x + 0.6780f*rgb.y + 0.0593f*rgb.z;
    rgb = rgb*pow(Ys, h_g - 1.0f);
  }
  return rgb;
}


vec3 eotf_pq(vec3 rgb, int inverse) {
  /* Apply the ST-2084 PQ Forward or Inverse EOTF
      ITU-R Rec BT.2100-2 https://www.itu.int/rec/R-REC-BT.2100
      ITU-R Rep BT.2390-9 https://www.itu.int/pub/R-REP-BT.2390
      Note: in the spec there is a normalization for peak display luminance. 
      For this function we assume the input is already normalized such that 1.0 = 10,000 nits
  */
  
  // const float Lp = 1.0f;
  const float m1 = 2610.0f/16384.0f;
  const float m2 = 2523.0f/32.0f;
  const float c1 = 107.0f/128.0f;
  const float c2 = 2413.0f/128.0f;
  const float c3 = 2392.0f/128.0f;

  if (inverse == 1) {
    // rgb /= Lp;
    rgb = spowf3(rgb, m1);
    rgb = spowf3((c1 + c2*rgb)/(1.0f + c3*rgb), m2);
  } else {
    rgb = spowf3(rgb, 1.0f/m2);
    rgb = spowf3((rgb - c1)/(c2 - c3*rgb), 1.0f/m1);
    // rgb *= Lp;
  }
  return rgb;
}


vec3 narrow_hue_angles(vec3 v) {
  return vec3(
    min(2.0f, max(0.0f, v.x - (v.y + v.z))),
    min(2.0f, max(0.0f, v.y - (v.x + v.z))),
    min(2.0f, max(0.0f, v.z - (v.x + v.y))));
}

float tonescale(float x, float m, float s, float c, int invert) {
  if (invert == 0) {
    return spowf(m*x/(x + s), c);
  } else {
    float ip = 1.0f/c;
    return spowf(s*x, ip)/(m - spowf(x, ip));
  }
}

float flare(float x, float fl, int invert) {
  if (invert == 0) {
    return spowf(x, 2.0f)/(x+fl);
  } else {
    return (x + sqrt(x*(4.0f*fl + x)))/2.0f;
  }
}

// https://www.desmos.com/calculator/gfubm2kvlu
float powerp(float x, float p, float m) {
  float y = x <= 0.0f ? x : x*spowf(spowf(x/m, 1.0f/p) + 1.0f, -p);
  return y;
}

// https://www.desmos.com/calculator/jrff9lrztn
float powerptoe(float x, float p, float m, float t0) {
  float y = x > t0 ? x : (x-t0)*spowf(spowf((t0-x)/(t0-m), 1.0f/p) + 1.0f, -p) + t0;
  return y;
}


vec3 transform(float p_R, float p_G, float p_B)
{
  // **************************************************
  // Parameter Setup
  // --------------------------------------------------

  // Dechroma
  float dch = 0.4f;

  // Chroma contrast
  float chc_p = 1.2f; // amount of contrast
  float chc_m = 0.5f; // pivot of contrast curve

  // Tonescale parameters
  float c = 1.1f; // contrast
  float fl = 0.1f; // flare/glare compensation

  // Weights: controls the "vibrancy" of each channel, and influences all other aspects of the display-rendering.
  vec3 weights = vec3(0.25f, 0.45f, 0.3f) / 2.0;

  // Hue Shift RGB controls
  vec3 hs = vec3(0.3f, -0.1f, -0.5f);

  // Input gamut conversion matrix (CAT02 chromatic adaptation to D65)
  mat3 in_to_xyz;
  if (in_gamut == i_xyz) in_to_xyz = identity();
  else if (in_gamut == i_ap0) in_to_xyz = matrix_ap0_to_xyz;
  else if (in_gamut == i_ap1) in_to_xyz = matrix_ap1_to_xyz;
  else if (in_gamut == i_p3d65) in_to_xyz = matrix_p3d65_to_xyz;
  else if (in_gamut == i_rec2020) in_to_xyz = matrix_rec2020_to_xyz;
  else if (in_gamut == i_rec709) in_to_xyz = matrix_rec709_to_xyz;
  else if (in_gamut == i_awg3) in_to_xyz = matrix_arriwg3_to_xyz;
  else if (in_gamut == i_awg4) in_to_xyz = matrix_arriwg4_to_xyz;
  else if (in_gamut == i_rwg) in_to_xyz = matrix_redwg_to_xyz;
  else if (in_gamut == i_sgamut3) in_to_xyz = matrix_sonysgamut3_to_xyz;
  else if (in_gamut == i_sgamut3cine) in_to_xyz = matrix_sonysgamut3cine_to_xyz;
  else if (in_gamut == i_bmdwg) in_to_xyz = matrix_bmdwg_to_xyz;
  else if (in_gamut == i_egamut) in_to_xyz = matrix_egamut_to_xyz;
  else if (in_gamut == i_davinciwg) in_to_xyz = matrix_davinciwg_to_xyz;

  mat3 xyz_to_display;
  if (display_gamut == Rec709) xyz_to_display = matrix_xyz_to_rec709;
  else if (display_gamut == P3D65) xyz_to_display = matrix_xyz_to_p3d65;
  else if (display_gamut == Rec2020) xyz_to_display = matrix_xyz_to_rec2020;

  int eotf;
  if (EOTF == lin)          eotf = 0;
  else if (EOTF == srgb)    eotf = 1;
  else if (EOTF == rec1886) eotf = 2;
  else if (EOTF == dci)     eotf = 3;
  else if (EOTF == pq)      eotf = 4;
  else if (EOTF == hlg)     eotf = 5;
  
  /* Display Scale ---------------*
      Remap peak white in display linear depending on the selected inverse EOTF.
      In our tonescale model, 1.0 is 100 nits, and as we scale up peak display luminance (Lp),
      we multiply up by the same amount. So if Lp=1,000, peak output of the tonescale model
      will be 10.0.

      So in ST2084 PQ, 1.0 is 10,000 nits, so we need to divide by 100 to fit out output into the 
      container.

      Similarly in HLG, 1.0 is 1,000 nits, so we need to divide by 10.

      If we are in an SDR mode, instead we just scale the peak so it hits display 1.0.
  */
  const float ds = eotf == 4 ? Lp/10000.0f : eotf == 5 ? Lp/1000.0f : 1.0f;
  
  /* Tonescale Parameters 
      ----------------------
    For the tonescale compression function, we use one inspired by the wisdom shared by Daniele Siragusano
    on the tonescale thread on acescentral: https://community.acescentral.com/t/output-transform-tone-scale/3498/224

    This is a variation which puts the power function _after_ the display-linear scale, which allows a simpler and exact
    solution for the intersection constraints. The resulting function is pretty much identical to Daniele's but simpler.
    Here is a desmos graph with the math. https://www.desmos.com/calculator/hglnae2ame

    And for more info on the derivation, see the "Michaelis-Menten Constrained" Tonescale Function here:
    https://colab.research.google.com/drive/1aEjQDPlPveWPvhNoEfK4vGH5Tet8y1EB#scrollTo=Fb_8dwycyhlQ

    For the user parameter space, we include the following creative controls:
    - Lp: display peak luminance. This sets the display device peak luminance and allows rendering for HDR.
    - contrast: This is a pivoted power function applied after the hyperbolic compress function, 
        which keeps middle grey and peak white the same but increases contrast in between.
    - flare: Applies a parabolic toe compression function after the hyperbolic compression function. 
        This compresses values near zero without clipping. Used for flare or glare compensation.
    - gb: Grey Boost. This parameter controls how many stops to boost middle grey per stop of peak luminance increase.

    Notes on the other non user-facing parameters:
    - (px, py): This is the peak luminance intersection constraint for the compression function.
        px is the input scene-linear x-intersection constraint. That is, the scene-linear input value 
        which is mapped to py through the compression function. By default this is set to 128 at Lp=100, and 256 at Lp=1000.
        Here is the regression calculation using a logarithmic function to match: https://www.desmos.com/calculator/chdqwettsj
    - (gx, gy): This is the middle grey intersection constraint for the compression function.
        Scene-linear input value gx is mapped to display-linear output gy through the function.
        Why is gy set to 0.11696 at Lp=100? This matches the position of middle grey through the Rec709 system.
        We use this value for consistency with the Arri and TCAM Rec.1886 display rendering transforms.
  */

  // input scene-linear peak x intercept
  float px = 256.0*_logf(Lp)/_logf(100.0) - 128.0f;
  // output display-linear peak y intercept
  float py = Lp/100.0f;
  // input scene-linear middle grey x intercept
  float gx = 0.18f;
  // output display-linear middle grey y intercept
  float gy = 11.696f/100.0f*(1.0f + gb*_logf(py)/_logf(2.0f));
  // s0 and s are input x scale for middle grey intersection constraint
  // m0 and m are output y scale for peak white intersection constraint
  float s0 = flare(gy, fl, 1);
  float m0 = flare(py, fl, 1);
  float ip = 1.0f/c;
  float s = (px*gx*(pow(m0, ip) - pow(s0, ip)))/(px*pow(s0, ip) - gx*pow(m0, ip));
  float m = pow(m0, ip)*(s + px)/px;



  /* Rendering Code ------------------------------------------ */

  vec3 rgb = vec3(p_R, p_G, p_B);

  // Convert into display gamut
  rgb = vdot(in_to_xyz, rgb);
  rgb = vdot(xyz_to_display, rgb);

  /* Take the the weighted sum of RGB. The weights
      scale the vector of each color channel, controlling the "vibrancy".
      We use this as a vector norm for separating color and intensity.
  */ 
  weights *= rgb; // multiply rgb by weights
  float lum = max(1e-8f, weights.x + weights.y + weights.z); // take the norm

  // RGB Ratios
  vec3 rats = sdivf3f(rgb, lum);

  // Apply tonescale function to lum
  float ts;
  ts = tonescale(lum, m, s, c, 0);
  ts = flare(ts, fl, 0);

  // Normalize so peak luminance is at 1.0
  ts *= 100.0f/Lp;

  // Clamp ts to display peak
  ts = min(1.0f, ts);

  /* Gamut Compress ------------------------------------------ *
    Most of our data is now inside of the display gamut cube, but there may still be some gradient disruptions
    due to highly chromatic colors going outside of the display cube on the lower end and then being clipped
    whether implicitly or explicitly. To combat this, our last step is to do a soft clip or gamut compression.
    In RGB Ratios, 0,0,0 is the gamut boundary, and anything outside of gamut will have one or more negative 
    components. So to compress the gamut we use lift these negative values and compress them into a small range
    near 0. We use the "PowerP" hyperbolic compression function but it could just as well be anything.
  */
  rats.x = powerptoe(rats.x, 0.05f, -0.05f, 1.0f);
  rats.y = powerptoe(rats.y, 0.05f, -0.05f, 1.0f);
  rats.z = powerptoe(rats.z, 0.05f, -0.05f, 1.0f);

  /* Calculate RGB CMY hue angles from the input RGB.
    The classical way of calculating hue angle from RGB is something like this
    mx = max(r,g,b)
    mn = min(r,g,b)
    c = mx - mn
    hue = (c==0?0:r==mx?((g-b)/c+6)%6:g==mx?(b-r)/c+2:b==mx?(r-g)/c+4:0)
    With normalized chroma (distance from achromatic), being calculated like this
    chroma = (mx - mn)/mx
    chroma can also be calculated as 1 - mn/mx

    Here we split apart the calculation for hue and chroma so that we have access to RGB CMY
    individually without having to linear step extract the result again.

    To do this, we first calculate the "wide" hue angle: 
      wide hue RGB = (RGB - mn)/mx
      wide hue CMY = (mx - RGB)/mx
    and then "narrow down" the hue angle for each with channel subtraction (see narrow_hue_angles() function).
  */
  
  float mx = max(rats.x, max(rats.y, rats.z));
  float mn = min(rats.x, min(rats.y, rats.z));

  vec3 rats_h = sdivf3f(rats - mn, mx);
  rats_h = narrow_hue_angles(rats_h);

  // Calculate "Chroma" (the normalized distance from achromatic).
  float rats_ch = 1.0f - sdivf(mn, mx);


  /* Chroma Value Compression ------------------------------------------ *
      RGB ratios may be greater than 1.0, which can result in discontinuities in highlight gradients.
      We compensate for this by normalizing the RGB Ratios so that max(r,g,b) does not exceed 1, and then mix
      the result. The factor for the mix is derived from tonescale * chroma, then taking only the top end of
      this with a compression function, so that we normalize only bright and saturated pixels.
  */

  // Normalization mix factor based on ccf * rgb chroma, smoothing transitions between r->g hue gradients
  float chf = ts*max(spowf(rats_h.x, 2.0f), max(spowf(rats_h.y, 2.0f), spowf(rats_h.z, 2.0f)));
  
  float chf_m = 0.25f;
  float chf_p = 0.65f;
  chf = 1.0f - spowf(spowf(chf/chf_m, 1.0f/chf_p)+1.0f, -chf_p);

  // Max of rgb ratios
  float rats_mx = max(rats.x, max(rats.y, rats.z));

  // Normalized rgb ratios
  vec3 rats_n = sdivf3f(rats, rats_mx);

  // Mix based on chf
  rats = rats_n*chf + rats*(1.0f - chf);


  /* Chroma Compression ------------------------------------------ *
      Here we set up the chroma compression factor, used to lerp towards 1.0
      in RGB Ratios, thereby compressing color towards display peak.
      This factor is driven by ts, biased by a power function to control chroma compression amount `dch`.
  */
  // float ccf = 1.0f - pow(ts, 1.0f/dch);
  float ccf = 1.0f - (pow(ts, 1.0f/dch)*(1.0f-ts) + ts*ts);

  // Apply chroma compression to RGB Ratios
  rats = rats*ccf + 1.0f - ccf;


  /* Chroma Compression Hue Shift ------------------------------------------ *
      Since we compress chroma by lerping in a straight line towards 1.0 in rgb ratios, this can result in perceptual hue shifts
      due to the Abney effect. For example, pure blue compressed in a straight line towards achromatic appears to shift in hue towards purple.

      To combat this, and to add another important user control for image appearance, we add controls to curve the hue paths 
      as they move towards achromatic. We include only controls for primary colors: RGB. In my testing, it was of limited use to
      control hue paths for CMY.

      To accomplish this, we use the inverse of the chroma compression factor multiplied by the RGB hue angles as a factor
      for a lerp between the various rgb components.

      We don't include the toe chroma compression for this hue shift. It is mostly important for highlights.
  */
  vec3 hsf = ccf*rats_h;
  
  // Apply hue shift to RGB Ratios
  vec3 rats_hs = vec3(rats.x + hsf.z*hs.z - hsf.y*hs.y, rats.y + hsf.x*hs.x - hsf.z*hs.z, rats.z + hsf.y*hs.y - hsf.x*hs.x);

  // Mix hue shifted RGB ratios by ts, so that we shift where highlights were chroma compressed plus a bit.
  rats = rats_hs*ts + rats*(1.0f - ts);


  /* Chroma Contrast
      Without this step, mid-range chroma in shadows and midtones looks too grey and dead.
      This is common with chromaticity-linear view transforms.
      In order to improve skin-tone rendering and overal "vibrance" of the image, which we
      are used to seeing with per-channel style view transforms, we boost mid-range chroma
      in shadows and midtones using a "chroma contrast" setup.
      
      Basically we take classical chroma (distance from achromatic), we take the compressed tonescale curve, 
      and we apply a contrast to the tonescale curve mixed by a parabolic center extraction of chroma, 
      so that we do not boost saturation at grey (increases noise), nor do we boost saturation of highly
      saturated colors which might already be near the edge of the gamut volume.
  */
  float chc_f = 4.0f*rats_ch*(1.0f - rats_ch);
  float chc_sa = min(2.0f, sdivf(lum, chc_m*spowf(sdivf(lum, chc_m), chc_p)*chc_f + lum*(1.0f - chc_f)));
  float chc_L = 0.23f*rats.x + 0.69f*rats.y + 0.08f*rats.z; // Roughly P3 weights, doesn't matter
  
  // Apply mid-range chroma contrast saturation boost
  rats = chc_L*(1.0f - chc_sa) + rats*chc_sa;

  // Apply tonescale to RGB Ratios
  rgb = rats*ts;

  // Apply display scale
  rgb *= ds;

  // Clamp
  rgb = maxf3(0.0f, rgb);
  rgb = minf3(ds, rgb);

  // Apply inverse Display EOTF
  float eotf_p = 2.0f + eotf * 0.2f;
  if ((eotf > 0) && (eotf < 4)) {
    rgb = spowf3(rgb, 1.0f/eotf_p);
  } else if (eotf == 4) {
    rgb = eotf_pq(rgb, 1);
  } else if (eotf == 5) {
    rgb = eotf_hlg(rgb, 1);
  }
  
  return rgb;
}