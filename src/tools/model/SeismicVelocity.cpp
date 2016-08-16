#include "SeismicVelocity.h"

#include <cassert>
#include <cstdio>
#include <cmath>

static double const Layers[] = { -100., -300., -1000., -3000., -6000., -31000. };
// vp, vs, density
static double const LayerLameParams[][3] = {
  {1925.9102873142, 876.2568612893, 2300.0000000000},
  {3473.1266454957, 1736.5633227478, 2300.0000000000},
  {4224.2887430448, 2296.9042711374, 2600.0000000000},
  {5100.0000000000, 2800.0000000000, 2700.0000000000},
  {5922.6430205131, 3288.7200310316, 2870.0000000000},
  {5933.4162648030, 3374.1453944454, 3500.0000000000},
  {7800.0000000000, 4500.0000000000, 3200.0000000000}
};

/*static double const Layers[] = { -2000. };
// vp, vs, density
static double const LayerLameParams[][3] = {
  {4000.0, 2000.0, 2600.0},
  {6000.0, 3464.0, 2700.0}
};*/

void get_material_parameters(double* lat, double* lon, double* height, unsigned numberOfPoints, unsigned mask, double* values)
{
  assert((mask & ~(Vp | Vs | Density)) == 0);
  
  unsigned params[3];
  unsigned nParams = 0;
  if (mask & Vp) params[nParams++] = 0;
  if (mask & Vs) params[nParams++] = 1;
  if (mask & Density) params[nParams++] = 2;

  unsigned numberOfLayers = sizeof(Layers) / sizeof(double);
  for (unsigned c = 0; c < numberOfPoints; ++c) {
    unsigned l = 0;
    while (l < numberOfLayers && height[c] < Layers[l]) {
      ++l;
    }
    for (unsigned i = 0; i < nParams; ++i) {
      values[c * nParams + i] = LayerLameParams[l][ params[i] ];
    }
  }  
}
