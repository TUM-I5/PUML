#ifndef SEISMIC_VELOCITY_H_
#define SEISMIC_VELOCITY_H_

double AsagiDerived(int, double x, double y, double z, int nx, int ny, int nz, double *xNc, double *yNc,double *zNc, double *VpNc);
double landers61(int, double, double, double z);
double sumatra1223_high(int group, double, double, double z);
double sumatra1223_low(int group, double, double, double z);
double sumatra1224(int group, double, double, double z);

#endif
