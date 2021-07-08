#ifndef _LUTEXP_CU_H_
#define _LUTEXP_CU_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// overwrite weak attribute linkage!
void __config_lutexp_cuda(double *xrand, int nrand);
void __compute_lutexp_cuda(uint64_t iters);

#ifdef __cplusplus
}
#endif
#endif /* #ifndef _LUTEXP_CU_H_ */
