#include <cstdio>
#include <string>
#include <gsl\gsl_sf.h>
#include <gsl\gsl_const.h>
#include <gsl\gsl_randist.h>
#include <gsl\gsl_math.h>
#include <iostream>
#include <gsl/gsl_const_mksa.h>
#include <math.h>
#include <gsl\gsl_complex.h>
#include <gsl\gsl_complex_math.h>
#include <gsl/gsl_linalg.h>



#define TAG_TRANSMIT_GAIN 1
#define TAG_RECEIVED_GAIN 1
#define TAG_REFLECTION_COEFFICIENT 1
#define READER_TRANSMIT_GAIN 1
#define READER_RECEIVED_GAIN 1
#define FREQUENCY 920
#define LOSSCONST -27.55
#define PATH_LOSS_EXPONENT 1.8
#define GAUSSIAN_NOISE 5.2
#define TRANSMITTER_POWER 2

typedef struct {
	double x;
	double y;
}coordinate;

coordinate antennas[3];
coordinate object;

int hex2dec(std::string hex)
{
	int dec = hex[1] * 16;
	dec += hex[0] * 16 * 16;
	return dec;
}

int tansRSSI(std::string rssiHex)
{

	int rssiDec = hex2dec(rssiHex);
	int rssi = 130 - rssiDec;
	return rssi;
}

int setdefultAntennaCoordinate()
{
	antennas[0].x = 0;
	antennas[0].y = 0;
	antennas[1].x = 0;
	antennas[1].y = 2;
	antennas[2].x = 2;
	antennas[2].y = 2;
	return 0;
}

double CalcTransmitterdBm()
{
	double result;
	result = 10 * gsl_sf_log(1000 * TRANSMITTER_POWER);
	return result;
}

double CalcWaveLength()
{
	double result;
	result = GSL_CONST_MKSA_SPEED_OF_LIGHT / (FREQUENCY * 1000000);
	return result;
}

double CalcFreeSpacePathLoss(double distance)
{

	double fspl;
	double wavelength = CalcWaveLength();
	double exponent = wavelength / (4 * M_PI * distance);
	exponent = gsl_pow_4(exponent);
	fspl = READER_RECEIVED_GAIN*READER_TRANSMIT_GAIN*(TAG_RECEIVED_GAIN*TAG_REFLECTION_COEFFICIENT*TAG_TRANSMIT_GAIN)*exponent;
	//fspl = 10 * gsl_sf_log(fspl);
	return fspl;
}

double CalcDownLinkPathLos(double distance)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	double pl;
	double fspl = CalcFreeSpacePathLoss(1);
	pl = fspl + 10 * PATH_LOSS_EXPONENT * gsl_sf_log(distance) + gsl_ran_gaussian(r,GAUSSIAN_NOISE*GAUSSIAN_NOISE);
	return pl;
}

double CalcUpLinkPathLoss(double distance)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	double pl;
	double fspl = CalcFreeSpacePathLoss(2);
	pl = fspl + 10 * 2 * PATH_LOSS_EXPONENT * gsl_sf_log(distance) + gsl_ran_gaussian(r, GAUSSIAN_NOISE*GAUSSIAN_NOISE);
	return pl;
}

double powReal(double power)
{
	gsl_complex z;
	z.dat[0] = M_E;
	z.dat[1] = 0;
	gsl_complex result = gsl_complex_pow_real(z, power);
	return result.dat[0];
}

double CalcDistanceByRSSIsimple(double RSSI) 
{
	double distance;
	double exponent;
	double fspl = CalcFreeSpacePathLoss(2);
	const gsl_rng_type * T;
	gsl_rng * r;
	double transmitdBm = CalcTransmitterdBm();
	exponent = (transmitdBm - RSSI - fspl - gsl_ran_gaussian(r, GAUSSIAN_NOISE*GAUSSIAN_NOISE));
	exponent = exponent / 36;
	distance = powReal(exponent);
	return distance;
}

void CalcCoordinateByDistance()
{
	double a11 = 2 * (antennas[1].x - antennas[0].x);
	double a12 = 2 * (antennas[1].y - antennas[0].y);
	double a21 = 2 * (antennas[2].x - antennas[0].x);
	double a22 = 2 * (antennas[2].y - antennas[0].y);
	double a_data[] = { a11,a12,
						a21,a22 };
	double distance[3];
	double b11 = gsl_pow_2(distance[0]) - gsl_pow_2(distance[1]) - gsl_pow_2(antennas[0].x) + gsl_pow_2(antennas[1].x) - gsl_pow_2(antennas[0].y) + gsl_pow_2(antennas[1].y);
	double b12 = gsl_pow_2(distance[0]) - gsl_pow_2(distance[2]) - gsl_pow_2(antennas[0].x) + gsl_pow_2(antennas[2].x) - gsl_pow_2(antennas[0].y) + gsl_pow_2(antennas[2].y);
	double b_data[] = { b11,b12 };
	gsl_matrix_view m = gsl_matrix_view_array(a_data, 2, 2);
	gsl_vector_view b = gsl_vector_view_array(b_data, 2);
	gsl_vector *ans = gsl_vector_alloc(2);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(2);
	gsl_linalg_LU_decomp(&m.matrix, p, &s);
	gsl_linalg_LU_solve(&m.matrix, p, &b.vector, ans);
	object.x = ans->data[0];
	object.y = ans->data[1];
	gsl_permutation_free(p);
	gsl_vector_free(ans);
	return;
}

int main()
{
    return 0;
}