#include <stdio.h>
#include <qedfvcoef.h>

int main(void)
{
  qedfv_context *ctx = qedfv_create_context();

  dvec3 v = {0.999, 0., 0.};
  ctx->debug = true;
  ctx->nmax = 425;
  double s = qedfv_coef(0., v, 0.3, ctx);
  printf("%f\n", s);

  qedfv_destroy_context(ctx);

  return 0;
}
