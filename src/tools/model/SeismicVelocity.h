#ifndef SEISMIC_VELOCITY_H_
#define SEISMIC_VELOCITY_H_

extern "C" {
enum MaterialParameter {
  Vp      = (1u << 0),
  Vs      = (1u << 1),
  Density = (1u << 2),
  Qp      = (1u << 3),
  Qs      = (1u << 4)
};

void get_material_parameters(double* lat, double* lon, double* height, unsigned numberOfPoints, unsigned mask, double* values);

}

#endif
