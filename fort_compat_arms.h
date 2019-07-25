#ifndef FORT_COMPAT_ARMS
#define FORT_COMPAT_ARMS
int arms_VS_setup(int ninit, double b, ARMS_STATE **state);
int arms_vs_setup_(int *ninit, double *b, ARMS_STATE **state);
int arms_sample_(double *xsamp, int *nsamp, ARMS_STATE **state);
void arms_cleanup_();
double VS(double x, void *VS_data);
struct VS_param {
  double b;
};
#endif
