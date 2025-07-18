/* nausha.h version 1.0; B D McKay, January 2025.
 * See credits and documentation in nausha.c. */

#ifndef _NAUSHA_H_   /* only process this file once */
#define _NAUSHA_H_

#include "gtools.h"

#if HAVE_STDINT_H
typedef unsigned char nsword8;
typedef uint32_t nsword32;
typedef uint64_t nsword64;
#else
typedef unsigned char nsword8;
typedef unsigned int nsword32;  /* Must be 32 bits */
typedef unsigned long long nsword64;  /* Must be 64 bits */
#endif

typedef struct {
    nsword8 data[64];
    unsigned long long bitlen;
    nsword32 datalen;
    nsword32 state[8];
} SHA256_CTX;

#ifdef __cplusplus
extern "C" {
#endif

void sha256_init(SHA256_CTX *ctx);
void sha256_update(SHA256_CTX *ctx, const nsword8 data[], size_t len);
void sha256_update_small(SHA256_CTX *ctx, const int data[], size_t len);
void sha256_final(SHA256_CTX *ctx, nsword8 hash[]);
void sha256(nsword8 hash[], nsword8 data[], size_t len);

void shahash(graph *g, int m, int n, nsword8 hash[]);
void shahash_sg(sparsegraph *sg, nsword8 hash[]);

#ifdef __cplusplus
}
#endif

#endif  /* _NAUSHA_H_ */
