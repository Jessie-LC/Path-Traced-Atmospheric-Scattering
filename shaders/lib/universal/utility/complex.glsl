#define complexFloat vec2

complexFloat complexAdd(complexFloat a, complexFloat b) {
	complexFloat c;
	c.x=a.x+b.x;
	c.y=a.y+b.y;
	return c;
}

complexFloat complexAdd(float a, complexFloat b) {
	complexFloat c;
	c.x=a+b.x;
	c.y= +b.y;
	return c;
}

complexFloat complexAdd(complexFloat a, float b) {
	complexFloat c;
	c.x=a.x+b;
	c.y=  a.y;
	return c;
}

complexFloat complexSub(complexFloat a, complexFloat b) {
	complexFloat c;
	c.x=a.x-b.x;
	c.y=a.y-b.y;
	return c;
}

complexFloat complexSub(complexFloat a, float b) {
	complexFloat c;
	c.x=a.x-b;
	c.y=  a.y;
	return c;
}

complexFloat complexSub(float a, complexFloat b) {
	complexFloat c;
	c.x=a-b.x;
	c.y= -b.y;
	return c;
}

complexFloat complexMul(complexFloat a, complexFloat b) {
	complexFloat c;
	c.x=a.x*b.x-a.y*b.y;
	c.y=a.y*b.x+a.x*b.y;
	return c;
}

complexFloat complexMul(float x, complexFloat a) {
	complexFloat c;
	c.x=x*a.x;
	c.y=x*a.y;
	return c;
}

complexFloat complexMul(complexFloat x, float a) {
	complexFloat c;
	c.x=x.x*a;
	c.y=x.y*a;
	return c;
}

complexFloat complexSquare(complexFloat x) {
    return complexMul(x, x);
}

complexFloat complexConjugate(complexFloat z) {
	complexFloat c;
	c.x=z.x;
	c.y = -z.y;
	return c;
}

complexFloat complexDiv(complexFloat a, complexFloat b) {
	complexFloat c;
	float r,den;
	if (abs(b.x) >= abs(b.y)) {
		r=b.y/b.x;
		den=b.x+r*b.y;
		c.x=(a.x+r*a.y)/(den);
		c.y=(a.y-r*a.x)/(den);
	} else {
		r=b.x/b.y;
		den=b.y+r*b.x;
		c.x=(a.x*r+a.y)/(den);
		c.y=(a.y*r-a.x)/(den);
	}
	return c;
}

complexFloat complexDiv(complexFloat a, float b) {
    return complexDiv(a, complexFloat(b, 0.0));
}

complexFloat complexDiv(float a, complexFloat b) {
    return complexDiv(complexFloat(a, 0.0), b);
}

float complexAbs(complexFloat z) {
	return sqrt(z.x*z.x + z.y*z.y);
}

complexFloat complexSqrt(complexFloat z) {
	complexFloat c;
	float x,y,w,r;
	if ((z.x == 0.0) && (z.y == 0.0)) {
		c.x=0.0;
		c.y=0.0;
		return c;
	} else {
		x=abs(z.x);
		y=abs(z.y);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.x >= 0.0) {
			c.x=w;
			c.y=z.y/(2.0*w);
		} else {
			c.y=(z.y >= 0.0) ? w : -w;
			c.x=z.y/(2.0*c.y);
		}
		return c;
	}
}

complexFloat complexExp(complexFloat z) {
	return complexMul(exp(z.x), complexFloat(cos(z.y), sin(z.y)));
}

complexFloat complexLog(complexFloat z) {
    return complexFloat(0.5 * log(z.x * z.x + z.y * z.y), atan(z.y, z.x));
}

complexFloat complexPow(complexFloat z, float power) {
    float r = complexAbs(z);
    float theta = atan(z.y, z.x);
    theta *= power;
    r = pow(r, power);
    return complexFloat(cos(theta) * r, sin(theta) * r);
}

complexFloat complexCbrt(complexFloat z) {
    complexFloat root = complexPow(z, 1.0/3.0);
    return root;
}

complexFloat complexSinh(complexFloat z) {
    return complexFloat(sinh(z.x) * cos(z.y), cosh(z.x) * sin(z.y));
}
complexFloat complexCosh(complexFloat z) {
    return complexFloat(cosh(z.x) * cos(z.y), sinh(z.x) * sin(z.y));
}

complexFloat complexSin(complexFloat z) {
	z = complexDiv(complexSub(complexExp(complexMul(complexFloat(0.0, 1.0), z)), complexExp(complexMul(complexFloat(0.0, -1.0), z))), complexFloat(0.0, 2.0));
    return z;
}
complexFloat complexCos(complexFloat z) {
	z = complexDiv(complexAdd(complexExp(complexMul(complexFloat(0.0, 1.0), z)), complexExp(complexMul(complexFloat(0.0, -1.0), z))), complexFloat(2.0, 0.0));
    return z;
}
complexFloat complexArgument(complexFloat z) {
  return complexFloat(atan(z.x, z.y), 0.0);
}
complexFloat complexArcsin(complexFloat z) {
  return complexDiv(
	  complexLog(
		  complexAdd(
			  complexMul(complexFloat(0.0, 1.0), z), 
			  complexMul(
				  sqrt(
					  complexAbs(
						  complexSub(1.0, complexSquare(z))
						)
					),
			  complexExp(
				  complexMul(
					  complexFloat(0.0, 0.5), 
					  complexArgument(
						  complexSub(1.0, complexSquare(z))
					  	)
					)
				)
				)
			)
		), 
		complexFloat(0.0, 1.0)
	);
}

complexFloat makeNeg(complexFloat a) {
    return complexFloat(-a.x, -a.y);
}

complexFloat vectorToComplex(vec2 z) {
	return complexFloat(z.x, z.y);
}

vec2 complexToVector(complexFloat z) {
	return vec2(z.x, z.y);
}