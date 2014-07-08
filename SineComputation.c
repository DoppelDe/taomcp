/*
 ============================================================================
 Name        : SineComputation.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdint.h>

double taylorSine(double value, int accuracy) {
	double i = 1;
	double q = value * value;
	double sine = value;
	double nextTerm = value;
	double sign = -1;
	double accVal = 1 / pow(10, accuracy);

	while (fabs(nextTerm) > accVal) {
		nextTerm = nextTerm * q / (2 * i * (2 * i + 1));
		sine = sine + sign * nextTerm;
		sign = -sign;
		i++;
	}
	return sine;
}

double arcTans[] = { 0.7853981633974482789994908671360462903976,
		0.4636476090008060935154787784995278343558,
		0.2449786631268641434733268624768243171275,
		0.1243549945467614381566789916178095154464,
		0.0624188099959573500230547438150097150356,
		0.0312398334302682774421544564802388777025,
		0.0156237286204768312941615349132007395383,
		0.0078123410601011111439873069173245312413,
		0.0039062301319669717573901390750279460917,
		0.0019531225164788187584341550007138721412,
		0.0009765621895593194594364927496599193546,
		0.0004882812111948982899262139412144279049,
		0.0002441406201493617712447448120371973346,
		0.0001220703118936702078530659454358442417,
		0.0000610351561742087725935014541622791739,
		0.0000305175781155260957271547345159845577,
		0.0000152587890613157615423778681873479002,
		0.0000076293945311019699810389967098434028,
		0.0000038146972656064961417507561819428830,
		0.0000019073486328101869647792853193490892,
		0.0000009536743164059608441276310632217506,
		0.0000004768371582030888422810640820542760,
		0.0000002384185791015579736676881098325631,
		0.0000001192092895507806808997385635169597,
		0.0000000596046447753905522081060953335646,
		0.0000000298023223876953025738326493636679,
		0.0000000149011611938476545956387748939447,
		0.0000000074505805969238281250000000000000,
		0.0000000037252902984619140625000000000000,
		0.0000000018626451492309570312500000000000,
		0.0000000009313225746154785156250000000000,
		0.0000000004656612873077392578125000000000,
		0.0000000002328306436538696289062500000000,
		0.0000000001164153218269348144531250000000,
		0.0000000000582076609134674072265625000000,
		0.0000000000291038304567337036132812500000,
		0.0000000000145519152283668518066406250000,
		0.0000000000072759576141834259033203125000,
		0.0000000000036379788070917129516601562500,
		0.0000000000018189894035458564758300781250,
		0.0000000000009094947017729282379150390625,
		0.0000000000004547473508864641189575195312,
		0.0000000000002273736754432320594787597656,
		0.0000000000001136868377216160297393798828,
		0.0000000000000568434188608080148696899414,
		0.0000000000000284217094304040074348449707,
		0.0000000000000142108547152020037174224854,
		0.0000000000000071054273576010018587112427,
		0.0000000000000035527136788005009293556213,
		0.0000000000000017763568394002504646778107,
		0.0000000000000008881784197001252323389053,
		0.0000000000000004440892098500626161694527,
		0.0000000000000002220446049250313080847263,
		0.0000000000000001110223024625156540423632,
		0.0000000000000000555111512312578270211816,
		0.0000000000000000277555756156289135105908,
		0.0000000000000000138777878078144567552954,
		0.0000000000000000069388939039072283776477,
		0.0000000000000000034694469519536141888238,
		0.0000000000000000017347234759768070944119,
		0.0000000000000000008673617379884035472060,
		0.0000000000000000004336808689942017736030,
		0.0000000000000000002168404344971008868015,
		0.0000000000000000001084202172485504434007,
		0.0000000000000000000542101086242752217004,
		0.0000000000000000000271050543121376108502,
		0.0000000000000000000135525271560688054251,
		0.0000000000000000000067762635780344027125,
		0.0000000000000000000033881317890172013563,
		0.0000000000000000000016940658945086006781,
		0.0000000000000000000008470329472543003391,
		0.0000000000000000000004235164736271501695,
		0.0000000000000000000002117582368135750848,
		0.0000000000000000000001058791184067875424,
		0.0000000000000000000000529395592033937712,
		0.0000000000000000000000264697796016968856,
		0.0000000000000000000000132348898008484428,
		0.0000000000000000000000066174449004242214,
		0.0000000000000000000000033087224502121107,
		0.0000000000000000000000016543612251060553,
		0.0000000000000000000000008271806125530277,
		0.0000000000000000000000004135903062765138,
		0.0000000000000000000000002067951531382569,
		0.0000000000000000000000001033975765691285,
		0.0000000000000000000000000516987882845642,
		0.0000000000000000000000000258493941422821,
		0.0000000000000000000000000129246970711411,
		0.0000000000000000000000000064623485355705,
		0.0000000000000000000000000032311742677853,
		0.0000000000000000000000000016155871338926,
		0.0000000000000000000000000008077935669463,
		0.0000000000000000000000000004038967834732,
		0.0000000000000000000000000002019483917366,
		0.0000000000000000000000000001009741958683,
		0.0000000000000000000000000000504870979341,
		0.0000000000000000000000000000252435489671,
		0.0000000000000000000000000000126217744835,
		0.0000000000000000000000000000063108872418,
		0.0000000000000000000000000000031554436209,
		0.0000000000000000000000000000015777218104 };
int arcTansInt[] = { 843314856, 497837829, 263043836, 133525158, 67021686,
		33543515, 16775850, 8388437, 4194282, 2097149, 1048575, 524287, 262143,
		131071, 65535, 32767, 16383, 8191, 4095, 2047, 1023, 511, 255, 127, 63,
		31, 15, 7, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0 };
int64_t arcTansLongInt[] = { 3622009730124208800LL, 2138197196257991700LL,
		1129764674086562800LL, 573486187337333400LL, 287855949760867800LL,
		144068300052059500LL, 72051727179705000LL, 36028062613110100LL,
		18014304036978600LL, 9007186378227700LL, 4503595336597500LL,
		2251795520815100LL, 1125895612923900LL, 562945658978300LL,
		281470682005500LL, 140733193519100LL, 70364449275900LL,
		35180077154300LL, 17587891093500LL, 8791798063100LL, 4393751547900LL,
		2194728290300LL, 1095216661500LL, 545460847100LL, 270582939900LL,
		133143986300LL, 64424509500LL, 30064771100LL, 12884901900LL,
		4294967300LL, 0LL, 0, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL,
		0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0, 0LL,
		0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL,
		0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL,
		0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL, 0LL };
int kInversInt = 652032873;
int64_t kInversLongInt = 2800459870029452800LL;
int HALF_INTMAX = INT32_MAX / 2;
int64_t HALF_LONGINTMAX = INT64_MAX / 2;
double kInvers = 0.6072529350088812561694;
double pi = 3.1415926535897932384626433832795;

#if defined(i386) || defined(i486) || \
	defined(intel) || defined(x86) || defined(i86pc) || \
	defined(__alpha) || defined(__osf__)
#define __LITTLE_ENDIAN
#endif

#ifdef __LITTLE_ENDIAN
#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x
#else
#define __HI(x) *(int*)&x
#define __LO(x) *(1+(int*)&x)
#define __HIp(x) *(int*)x
#define __LOp(x) *(1+(int*)x)
#endif

static int two_over_pi[] = { 0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD,
		0xC0DB62, 0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
		0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129, 0xA73EE8,
		0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41, 0x3991D6, 0x398353,
		0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8, 0x97FFDE, 0x05980F, 0xEF2F11,
		0x8B5A0A, 0x6D1F6D, 0x367ECF, 0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D,
		0x7527BA, 0xC7EBE5, 0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F,
		0x8D5D08, 0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
		0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880, 0x4D7327,
		0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B, };
static int npio2_hw[] = { 0x3FF921FB, 0x400921FB, 0x4012D97C, 0x401921FB,
		0x401F6A7A, 0x4022D97C, 0x4025FDBB, 0x402921FB, 0x402C463A, 0x402F6A7A,
		0x4031475C, 0x4032D97C, 0x40346B9C, 0x4035FDBB, 0x40378FDB, 0x403921FB,
		0x403AB41B, 0x403C463A, 0x403DD85A, 0x403F6A7A, 0x40407E4C, 0x4041475C,
		0x4042106C, 0x4042D97C, 0x4043A28C, 0x40446B9C, 0x404534AC, 0x4045FDBB,
		0x4046C6CB, 0x40478FDB, 0x404858EB, 0x404921FB, };
static double zero = 0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
half = 5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
two24 = 1.67772160000000000000e+07, /* 0x41700000, 0x00000000 */
invpio2 = 6.36619772367581382433e-01, /* 0x3FE45F30, 0x6DC9C883 */
pio2_1 = 1.57079632673412561417e+00, /* 0x3FF921FB, 0x54400000 */
pio2_1t = 6.07710050650619224932e-11, /* 0x3DD0B461, 0x1A626331 */
pio2_2 = 6.07710050630396597660e-11, /* 0x3DD0B461, 0x1A600000 */
pio2_2t = 2.02226624879595063154e-21, /* 0x3BA3198A, 0x2E037073 */
pio2_3 = 2.02226624871116645580e-21, /* 0x3BA3198A, 0x2E000000 */
pio2_3t = 8.47842766036889956997e-32; /* 0x397B839A, 0x252049C1 */

static int init_jk[] = { 2, 3, 4, 6 };
static double PIo2[] = { 1.57079625129699707031e+00, /* 0x3FF921FB, 0x40000000 */
7.54978941586159635335e-08, /* 0x3E74442D, 0x00000000 */
5.39030252995776476554e-15, /* 0x3CF84698, 0x80000000 */
3.28200341580791294123e-22, /* 0x3B78CC51, 0x60000000 */
1.27065575308067607349e-29, /* 0x39F01B83, 0x80000000 */
1.22933308981111328932e-36, /* 0x387A2520, 0x40000000 */
2.73370053816464559624e-44, /* 0x36E38222, 0x80000000 */
2.16741683877804819444e-51, /* 0x3569F31D, 0x00000000 */
};

static double one = 1.0, twon24 = 5.96046447753906250000e-08; /* 0x3E700000, 0x00000000 */

/*
Copied from FDLIBM:

FDLIBM (Freely Distributable LIBM) is a C math library 
for machines that support IEEE 754 floating-point arithmetic. 
In this release, only double precision is supported.

FDLIBM is intended to provide a reasonably portable (see 
assumptions below), reference quality (below one ulp for
major functions like sin,cos,exp,log) math library 
(libm.a).  For a copy of FDLIBM, please see
	http://www.netlib.org/fdlibm/
or
	http://www.validlab.com/software/
*/
int __kernel_rem_pio2(double *x, double *y, int e0, int nx, int prec,
		const int *ipio2) {
	int jz, jx, jv, jp, jk, carry, n, iq[20], i, j, k, m, q0, ih;
	double z, fw, f[20], fq[20], q[20];

	/* initialize jk*/
	jk = init_jk[prec];
	jp = jk;

	/* determine jx,jv,q0, note that 3>q0 */
	jx = nx - 1;
	jv = (e0 - 3) / 24;
	if (jv < 0)
		jv = 0;
	q0 = e0 - 24 * (jv + 1);

	/* set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk] */
	j = jv - jx;
	m = jx + jk;
	for (i = 0; i <= m; i++, j++)
		f[i] = (j < 0) ? zero : (double) ipio2[j];

	/* compute q[0],q[1],...q[jk] */
	for (i = 0; i <= jk; i++) {
		for (j = 0, fw = 0.0; j <= jx; j++)
			fw += x[j] * f[jx + i - j];
		q[i] = fw;
	}

	jz = jk;
	recompute:
	/* distill q[] into iq[] reversingly */
	for (i = 0, j = jz, z = q[jz]; j > 0; i++, j--) {
		fw = (double) ((int) (twon24 * z));
		iq[i] = (int) (z - two24 * fw);
		z = q[j - 1] + fw;
	}

	/* compute n */
	z = scalbn(z, q0); /* actual value of z */
	z -= 8.0 * floor(z * 0.125); /* trim off integer >= 8 */
	n = (int) z;
	z -= (double) n;
	ih = 0;
	if (q0 > 0) { /* need iq[jz-1] to determine n */
		i = (iq[jz - 1] >> (24 - q0));
		n += i;
		iq[jz - 1] -= i << (24 - q0);
		ih = iq[jz - 1] >> (23 - q0);
	} else if (q0 == 0)
		ih = iq[jz - 1] >> 23;
	else if (z >= 0.5)
		ih = 2;

	if (ih > 0) { /* q > 0.5 */
		n += 1;
		carry = 0;
		for (i = 0; i < jz; i++) { /* compute 1-q */
			j = iq[i];
			if (carry == 0) {
				if (j != 0) {
					carry = 1;
					iq[i] = 0x1000000 - j;
				}
			} else
				iq[i] = 0xffffff - j;
		}
		if (q0 > 0) { /* rare case: chance is 1 in 12 */
			switch (q0) {
			case 1:
				iq[jz - 1] &= 0x7fffff;
				break;
			case 2:
				iq[jz - 1] &= 0x3fffff;
				break;
			}
		}
		if (ih == 2) {
			z = one - z;
			if (carry != 0)
				z -= scalbn(one, q0);
		}
	}

	/* check if recomputation is needed */
	if (z == zero) {
		j = 0;
		for (i = jz - 1; i >= jk; i--)
			j |= iq[i];
		if (j == 0) { /* need recomputation */
			for (k = 1; iq[jk - k] == 0; k++)
				; /* k = no. of terms needed */

			for (i = jz + 1; i <= jz + k; i++) { /* add q[jz+1] to q[jz+k] */
				f[jx + i] = (double) ipio2[jv + i];
				for (j = 0, fw = 0.0; j <= jx; j++)
					fw += x[j] * f[jx + i - j];
				q[i] = fw;
			}
			jz += k;
			goto recompute;
		}
	}

	/* chop off zero terms */
	if (z == 0.0) {
		jz -= 1;
		q0 -= 24;
		while (iq[jz] == 0) {
			jz--;
			q0 -= 24;
		}
	} else { /* break z into 24-bit if necessary */
		z = scalbn(z, -q0);
		if (z >= two24) {
			fw = (double) ((int) (twon24 * z));
			iq[jz] = (int) (z - two24 * fw);
			jz += 1;
			q0 += 24;
			iq[jz] = (int) fw;
		} else
			iq[jz] = (int) z;
	}

	/* convert integer "bit" chunk to floating-point value */
	fw = scalbn(one, q0);
	for (i = jz; i >= 0; i--) {
		q[i] = fw * (double) iq[i];
		fw *= twon24;
	}

	/* compute PIo2[0,...,jp]*q[jz,...,0] */
	for (i = jz; i >= 0; i--) {
		for (fw = 0.0, k = 0; k <= jp && k <= jz - i; k++)
			fw += PIo2[k] * q[i + k];
		fq[jz - i] = fw;
	}

	/* compress fq[] into y[] */
	switch (prec) {
	case 0:
		fw = 0.0;
		for (i = jz; i >= 0; i--)
			fw += fq[i];
		y[0] = (ih == 0) ? fw : -fw;
		break;
	case 1:
	case 2:
		fw = 0.0;
		for (i = jz; i >= 0; i--)
			fw += fq[i];
		y[0] = (ih == 0) ? fw : -fw;
		fw = fq[0] - fw;
		for (i = 1; i <= jz; i++)
			fw += fq[i];
		y[1] = (ih == 0) ? fw : -fw;
		break;
	case 3: /* painful */
		for (i = jz; i > 0; i--) {
			fw = fq[i - 1] + fq[i];
			fq[i] += fq[i - 1] - fw;
			fq[i - 1] = fw;
		}
		for (i = jz; i > 1; i--) {
			fw = fq[i - 1] + fq[i];
			fq[i] += fq[i - 1] - fw;
			fq[i - 1] = fw;
		}
		for (fw = 0.0, i = jz; i >= 2; i--)
			fw += fq[i];
		if (ih == 0) {
			y[0] = fq[0];
			y[1] = fq[1];
			y[2] = fw;
		} else {
			y[0] = -fq[0];
			y[1] = -fq[1];
			y[2] = -fw;
		}
	}
	return n & 7;
}

int __ieee754_rem_pio2(double x, double *y) {
	double z, w, t, r, fn;
	double tx[3];
	int e0, i, j, nx, n, ix, hx;

	hx = __HI(x); /* high word of x */
	ix = hx & 0x7fffffff;
	if (ix <= 0x3fe921fb) /* |x| ~<= pi/4 , no need for reduction */
	{
		y[0] = x;
		y[1] = 0;
		return 0;
	}
	if (ix < 0x4002d97c) { /* |x| < 3pi/4, special case with n=+-1 */
		if (hx > 0) {
			z = x - pio2_1;
			if (ix != 0x3ff921fb) { /* 33+53 bit pi is good enough */
				y[0] = z - pio2_1t;
				y[1] = (z - y[0]) - pio2_1t;
			} else { /* near pi/2, use 33+33+53 bit pi */
				z -= pio2_2;
				y[0] = z - pio2_2t;
				y[1] = (z - y[0]) - pio2_2t;
			}
			return 1;
		} else { /* negative x */
			z = x + pio2_1;
			if (ix != 0x3ff921fb) { /* 33+53 bit pi is good enough */
				y[0] = z + pio2_1t;
				y[1] = (z - y[0]) + pio2_1t;
			} else { /* near pi/2, use 33+33+53 bit pi */
				z += pio2_2;
				y[0] = z + pio2_2t;
				y[1] = (z - y[0]) + pio2_2t;
			}
			return -1;
		}
	}
	if (ix <= 0x413921fb) { /* |x| ~<= 2^19*(pi/2), medium size */
		t = fabs(x);
		n = (int) (t * invpio2 + half);
		fn = (double) n;
		r = t - fn * pio2_1;
		w = fn * pio2_1t; /* 1st round good to 85 bit */
		if (n < 32 && ix != npio2_hw[n - 1]) {
			y[0] = r - w; /* quick check no cancellation */
		} else {
			j = ix >> 20;
			y[0] = r - w;
			i = j - (((__HI(y[0])) >> 20) & 0x7ff);
			if (i > 16) { /* 2nd iteration needed, good to 118 */
				t = r;
				w = fn * pio2_2;
				r = t - w;
				w = fn * pio2_2t - ((t - r) - w);
				y[0] = r - w;
				i = j - (((__HI(y[0])) >> 20) & 0x7ff);
				if (i > 49) { /* 3rd iteration need, 151 bits acc */
					t = r; /* will cover all possible cases */
					w = fn * pio2_3;
					r = t - w;
					w = fn * pio2_3t - ((t - r) - w);
					y[0] = r - w;
				}
			}
		}
		y[1] = (r - y[0]) - w;
		if (hx < 0) {
			y[0] = -y[0];
			y[1] = -y[1];
			return -n;
		} else
			return n;
	}
	/*
	 * all other (large) arguments
	 */
	if (ix >= 0x7ff00000) { /* x is inf or NaN */
		y[0] = y[1] = x - x;
		return 0;
	}
	/* set z = scalbn(|x|,ilogb(x)-23) */
	__LO(z) = __LO(x);
	e0 = (ix >> 20) - 1046; /* e0 = ilogb(z)-23; */
	__HI(z) = ix - (e0 << 20);
	for (i = 0; i < 2; i++) {
		tx[i] = (double) ((int) (z));
		z = (z - tx[i]) * two24;
	}
	tx[2] = z;
	nx = 3;
	while (tx[nx - 1] == zero)
		nx--; /* skip zero term */
	n = __kernel_rem_pio2(tx, y, e0, nx, 2, two_over_pi);
	if (hx < 0) {
		y[0] = -y[0];
		y[1] = -y[1];
		return -n;
	}
	return n;
}

int reduceRange(double value, double* reducedValue) {
	double result[2];
	int n = __ieee754_rem_pio2(value, result);
	*reducedValue = result[0];
	return n;
}

double cordicSineOnDouble(double value) {
	int n = 0;
	double xn = kInvers;
	double yn = 0;
	double zn;
	int calc = reduceRange(value, &zn);
	double lastDiff = 1;
	double acc = 1 / pow(10, 10);
	double currentPow = 2;
	double factor = 0.5;
	while (fabs(lastDiff) >= acc && n < 50) {
		double d;
		if (zn >= 0) {
			d = 1;
		} else {
			d = -1;
		}
		currentPow *= factor;
		lastDiff = d * xn * currentPow;
		double xn1 = xn - d * yn * currentPow;
		yn += lastDiff;
		zn = zn - d * arcTans[n];
		xn = xn1;
		n++;
	}
	switch (calc & 3) {
	case 0:
		return yn;
	case 1:
		return xn;
	case 2:
		return -yn;
	default:
		return -xn;
	}
}

int cordicSineOnInt(int value, int calc) {
	int n = 0;
	int xn = kInversInt;
	int yn = 0;
	int zn = value;
	int lastDiff = HALF_INTMAX;
	int acc = 1;
	while (abs(lastDiff) >= acc && n < 50) {
		int d;
		if (zn >= 0) {
			d = 1;
		} else {
			d = -1;
		}
		lastDiff = d * (xn >> n);
		int xn1 = xn - d * (yn >> n);
		yn += lastDiff;
		zn = zn - d * arcTansInt[n];
		xn = xn1;
		n++;
	}
	switch (calc & 3) {
	case 0:
		return yn;
	case 1:
		return xn;
	case 2:
		return -yn;
	default:
		return -xn;
	}
}

double cordicBasedOnInt(double value) {

	double reducedValue;
	int calc = reduceRange(value, &reducedValue);
	int result = cordicSineOnInt(reducedValue * HALF_INTMAX, calc);
	return (double) result / HALF_INTMAX;
}

int64_t cordicSineOnLongInt(int64_t value) {
	int64_t n = 0;
	int64_t xn = kInversLongInt;
	int64_t yn = 0;
	int64_t zn = value;
	int64_t lastDiff = HALF_LONGINTMAX;
	int64_t acc = 1;
	while (abs(lastDiff) >= acc && n < 50) {
		int64_t d;
		if (zn >= 0) {
			d = 1;
		} else {
			d = -1;
		}
		lastDiff = d * (xn >> n);
		int64_t xn1 = xn - d * (yn >> n);
		yn += lastDiff;
		zn = zn - d * arcTansLongInt[n];
		xn = xn1;
		n++;
	}
	return yn;
}

double cordicBasedOnLongInt(double value) {

	int64_t result = cordicSineOnLongInt(1 * HALF_LONGINTMAX);
	return (double) result / HALF_LONGINTMAX;
}

double calculatePi() {

	int max = 10000;
	double sum = 0;
	int i;
	for (i = 1; i < max; ++i) {
		sum += 6.0 / (i * i);
	}
	return sqrt(sum);
}

int main(void) {
	double value = 10000000;
	int i;
	clock_t begin, between, end;
	begin = clock();
	double r1, r2;
	for (i = 0; i < 10000000; i++)
		r1 = cordicSineOnDouble(value);
	between = clock();
	for (i = 0; i < 10000000; i++)
		r2 = cordicBasedOnInt(value);
	end = clock();
	float t1 = between - begin;
	float t2 = end - between;
	printf("%.20F %.20F\n", t1, t2);
	printf("%.10F %.10F %.10F", r1, r2, sin(value));
	return 0;
}
