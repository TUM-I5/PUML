#include "SeismicVelocity.h"
#include <iostream>

double AsagiDerived(int, double x, double y, double z, int nx, int ny, int nz, double *xNc, double *yNc,double *zNc, double *VpNc)
{
   int i1,j1,k1;
   double ax,ay,az;
   for (i1 = 0; i1 < nx; i1++)
   {
      if (xNc[i1]>=x) 
         break;
   }
   for (j1 = 0; j1 < ny; j1++)
   {
      if (yNc[j1]>=y) 
         break;
   }
   for (k1 = 0; k1 < nz; k1++)
   {
      if (zNc[k1]>=z) 
         break;
   }
   

   if (i1==0) {
       std::cerr<<"i1 =0 "<<x<<" "<<xNc[0];
       i1=1;
       ax=0.0;
   } else if (x>xNc[nx-1]) {
       std::cerr<<"i1 =nx-1 "<<x<<" "<<xNc[nx-1];
       i1=nx-1;
       ax = 1.0;
   } else {
       ax = (x-xNc[i1-1])/(xNc[i1]-xNc[i1-1]);
   }

   if (j1==0) {
       std::cerr<<"j1 =0 "<<y<<" "<<yNc[0];
       j1=1;
       ay=0.0;
   } else if (y>yNc[ny-1]) {
       std::cerr<<"j1 =ny-1 "<<y<<" "<<yNc[ny-1];
       j1=ny-1;
       ay = 1.0;
   } else {
       ay = (y-yNc[j1-1])/(yNc[j1]-yNc[j1-1]);
   }

   if (k1==0) {
       std::cerr<<"k1 =0 "<<z<<" "<<zNc[0];
       k1=1;
       az=0.0;
   } else if (z>zNc[nz-1]) {
       std::cerr<<"k1 =nz-1 "<<z<<" "<<zNc[nz-1];
       k1=nz-1;
       az = 1.0;
   } else {
       az = (z-zNc[k1-1])/(zNc[k1]-zNc[k1-1]);
   }
   if ((std::min(ax,std::min(ay,az))<0)|(std::max(ax,std::max(ay,az))>1))
      std::cerr<<ax<<" "<<ay<<" "<<az;
   int na=ny*nz;
   double vp  = ax*ay*az*VpNc[i1*na+j1*nz+k1] + ax*ay*(1.0-az)*VpNc[i1*na+j1*nz+k1-1] + 
                                   ax*(1.0-ay)*az*VpNc[i1*na+(j1-1)*nz+k1] + ax*(1.0-ay)*(1.0-az)*VpNc[i1*na+(j1-1)*nz+k1-1] + 
                                   (1.0-ax)*ay*az*VpNc[(i1-1)*na+j1*nz+k1] + (1.0-ax)*ay*(1.0-az)*VpNc[(i1-1)*na+j1*nz+k1-1] + 
                                   (1.0-ax)*(1.0-ay)*az*VpNc[(i1-1)*na+(j1-1)*nz+k1] + (1.0-ax)*(1.0-ay)*(1.0-az)*VpNc[(i1-1)*na+(j1-1)*nz+k1-1];
  
  return vp;
}

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

/*
BedrockVelModel(1,:) = (/  -6d3, 2550d0,18589500000d0,26571000000d0/) => 5000 m/s
BedrockVelModel(2,:) = (/  -8d3, 2850d0,39016500000d0,42379500000d0/) => 6500 m/s
BedrockVelModel(3,:) = (/ -12d3, 3050d0,50027625000d0,53695250000d0/) => 7100 m/s
BedrockVelModel(4,:) = (/-6d3,2720d0,33320000000d0,31280000000d0/)    => 6000 m/s
BedrockVelModel(5,:) = (/-12d3,2860d0,41298400000d0,41984800000d0/)   => 6600 m/s
BedrockVelModel(6,:) = (/-23d3,3050d0,46390500000d0,60969500000d0/)   => 7100 m/s
BedrockVelModel(7,:) = (/ -5d10, 3330d0,65942325000d0,81235350000d0/) => 8000 m/s  */
double sumatra1223_high(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 4: // LVZ
    case 2: // layer 1
      cp = 5000.0;
      break;
    case 3: // layer 2
      cp = 6500.0;
      break;
    case 7: // layer 3
      cp = 7100.0;
      break;
    case 6: // layer 4
      cp = 8000.0;
      break;
    case 1: // continental
    case 5: // box
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

double sumatra1223_low(int group, double, double, double z)
{
  double cp = -1.0;
  switch (group) {
    case 2: // LVZ
    case 7: // layer 1
      cp = 5000.0;
      break;
    case 4: // layer 2
      cp = 6500.0;
      break;
    case 3: // layer 3
      cp = 7100.0;
      break;
    case 6: // layer 4
      cp = 8000.0;
      break;
    case 1: // continental
    case 5: // box
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
