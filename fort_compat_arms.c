#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arms.h"
#include "fort_compat_arms.h"

double VS(double x, void *VS_data)
{
  double y;
  double e, b;
  e = 3.0/2.0;
  b = ((struct VS_param*) VS_data)->b;

  y = -(b/pow(x,e)) - e*log(x);

  return y;
}

int arms_VS_setup(int ninit, double b, ARMS_STATE **state)
{
  int err;
  double xl = 0.0, xr = 1.0, xprev = 0.9;
  struct VS_param *VS_data;
  VS_data = (struct VS_param*)malloc(sizeof(struct VS_param));
  VS_data->b = b;
  int dometrop;
  if(b >= 0.4) {
      dometrop = 0;
  } else {
      dometrop = 1;
  }

  err = arms_simple_setup(ninit, &xl, &xr, VS, VS_data, dometrop, &xprev, state);
  if(*state == NULL) printf("Null state, bro.\n");
  return err;
}

int arms_vs_setup_(int *ninit, double *b, ARMS_STATE **state)
{
  int ninit_ = *ninit;
  double b_ = *b;

  return arms_VS_setup(ninit_, b_, state);
}

int arms_sample_(double *xsamp, int *nsamp, ARMS_STATE **state)
{
    int nsamp_ = *nsamp;
    return arms_sample(xsamp, nsamp_, *state);
}

void arms_cleanup_(ARMS_STATE **state)
{
    arms_cleanup(*state);
}
