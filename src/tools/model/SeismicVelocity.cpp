#include "SeismicVelocity.h"

double landers61(int, double, double, double z)
{
  if (z > -0.1) {
    return 1925.9102873142;
  }
  if (z >= -0.3) {
    return 3473.1266454957;
  }
  if (z >= -1.0) {
    return 4224.2887430448;
  }
  if (z >= -3.0) {
    return 5100.0000000000;
  }
  if (z >= -6.0) {
    return 5922.6430205131;
  }
  if (z >= -31.0) {
    return 5933.4162648030;
  } else {
    return 7800.0000000000;
  }
  return -1.0;
}

double sumatra1224(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 1:
      cp = 6500.0;
      break;
    case 2:
      cp = 5000.0;
      break;
    case 3:
      cp = 7100.0;
      break;
    case 4:
    case 5:
      cp = 8000.0;
      break;
    case 6:
      if (z > -6000.0) {
        cp = 6000.0;
      } else if (z >= -12000.0) {
        cp = 6600.0;
      } else if (z > -23000.0) {
        cp = 7100.0;
      } else {
        cp = 8000.0;
      }
      break;
    default:
      break;
  }
  return cp;
}
