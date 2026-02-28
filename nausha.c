/* nausha.c version 1.1; B D McKay, August 2025 */

#include "nausha.h"
#if HAVE_ENDIAN_H && !defined(__STDC_ENDIAN_NATIVE__)
#include <endian.h>
#endif

/****************************************************************************
* This implementation of SHA256 is by Brad Conte, released to public domain.
* Author:     Brad Conte (brad AT bradconte.com) */

/* SHA256 is a 64-byte cryptographic hash code. The general method for
 * computing it is:
 *    sha256_init(SHA256_CTX *ctx);
 *    sha256_update(SHA256_CTX *ctx, const nsword8 data[], size_t len);
 *       any number of times, where data is an array of len bytes.
 *    sha256_final(SHA256_CTX *ctx, nsword8 hash[])
 * Then the result is in hash[0..63].
 */

#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

static const nsword32 k[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

#if defined(__STDC_ENDIAN_NATIVE__) && defined(__STDC_ENDIAN_LITTLE__)
#if __STDC_ENDIAN_NATIVE__ == __STDC_ENDIAN_LITTLE__
#define ISBIGENDIAN 0
#else
#define ISBIGENDIAN 1
#endif
#elif defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN)
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define ISBIGENDIAN 0
#else
#define ISBIGENDIAN 1
#endif
#elif defined(__DARWIN_BYTE_ORDER) && defined(__DARWIN_LITTLE_ENDIAN)
#if __DARWIN_BYTE_ORDER == __DARWIN_LITTLE_ENDIAN
#define ISBIGENDIAN 0
#else
#define ISBIGENDIAN 1
#endif
#else
static int
isbigendian(void)
{
    unsigned int x = 1;
    unsigned char *c = (char*)&x;
    return (*c == (char)1 ? 0 : 1);
}
#define ISBIGENDIAN isbigendian() 
#endif

static void
sha256_transform(SHA256_CTX *ctx, const nsword8 data[])
{
    nsword32 a, b, c, d, e, f, g, h, i, j, t1, t2, m[64];

    for (i = 0, j = 0; i < 16; ++i, j += 4)
            m[i] = (data[j] << 24) | (data[j + 1] << 16)
				| (data[j + 2] << 8) | (data[j + 3]);
    for ( ; i < 64; ++i)
            m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

    a = ctx->state[0];
    b = ctx->state[1];
    c = ctx->state[2];
    d = ctx->state[3];
    e = ctx->state[4];
    f = ctx->state[5];
    g = ctx->state[6];
    h = ctx->state[7];

    for (i = 0; i < 64; ++i) {
        t1 = h + EP1(e) + CH(e,f,g) + k[i] + m[i];
        t2 = EP0(a) + MAJ(a,b,c);
        h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
    }

    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}

void
sha256_init(SHA256_CTX *ctx)
{
    ctx->datalen = 0;
    ctx->bitlen = 0;
    ctx->state[0] = 0x6a09e667;
    ctx->state[1] = 0xbb67ae85;
    ctx->state[2] = 0x3c6ef372;
    ctx->state[3] = 0xa54ff53a;
    ctx->state[4] = 0x510e527f;
    ctx->state[5] = 0x9b05688c;
    ctx->state[6] = 0x1f83d9ab;
    ctx->state[7] = 0x5be0cd19;
}

void
sha256_update(SHA256_CTX *ctx, const nsword8 data[], size_t len)
{
    size_t i;

    for (i = 0; i < len; ++i) {
        ctx->data[ctx->datalen] = data[i];
        ctx->datalen++;
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

void
sha256_update_small(SHA256_CTX *ctx, const int data[], size_t len)
/* Version for values < 2^16, ignoring all but the lower order two bytes
 * of each int. */
{
    size_t i;
    nsword8 x;

    for (i = 0; i < len; ++i) {
        x = data[i] & 0xFF;
        ctx->data[ctx->datalen] = x;
        ctx->datalen++;
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
        x = ((nsword32)data[i] >> 8) & 0xFF;
        ctx->data[ctx->datalen] = x;
        ctx->datalen++;
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

void
sha256_final(SHA256_CTX *ctx, nsword8 hash[])
{
    nsword32 i;

    i = ctx->datalen;

    // Pad whatever data is left in the buffer.
    if (ctx->datalen < 56) {
        ctx->data[i++] = 0x80;
        while (i < 56)
            ctx->data[i++] = 0x00;
    }
    else {
        ctx->data[i++] = 0x80;
        while (i < 64)
                ctx->data[i++] = 0x00;
        sha256_transform(ctx, ctx->data);
        memset(ctx->data, 0, 56);
    }

    // Append to the padding the total message's length in bits and transform.
    ctx->bitlen += ctx->datalen * 8;
    ctx->data[63] = ctx->bitlen;
    ctx->data[62] = ctx->bitlen >> 8;
    ctx->data[61] = ctx->bitlen >> 16;
    ctx->data[60] = ctx->bitlen >> 24;
    ctx->data[59] = ctx->bitlen >> 32;
    ctx->data[58] = ctx->bitlen >> 40;
    ctx->data[57] = ctx->bitlen >> 48;
    ctx->data[56] = ctx->bitlen >> 56;
    sha256_transform(ctx, ctx->data);

    // Since SHA uses big endian, we need the endianness of the current
    // runtime to copy the result into the final hash.
 
    if (ISBIGENDIAN)
        for (i = 0; i < 4; ++i) {
            hash[i]      = (ctx->state[1] >> (i * 8)) & 0x000000ff;
            hash[i + 4]  = (ctx->state[0] >> (i * 8)) & 0x000000ff;
            hash[i + 8]  = (ctx->state[3] >> (i * 8)) & 0x000000ff;
            hash[i + 12] = (ctx->state[2] >> (i * 8)) & 0x000000ff;
            hash[i + 16] = (ctx->state[5] >> (i * 8)) & 0x000000ff;
            hash[i + 20] = (ctx->state[4] >> (i * 8)) & 0x000000ff;
            hash[i + 24] = (ctx->state[7] >> (i * 8)) & 0x000000ff;
            hash[i + 28] = (ctx->state[6] >> (i * 8)) & 0x000000ff;
        }
    else
        for (i = 0; i < 4; ++i) {
            hash[i]      = (ctx->state[0] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 4]  = (ctx->state[1] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 8]  = (ctx->state[2] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 20] = (ctx->state[5] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 24] = (ctx->state[6] >> (24 - i * 8)) & 0x000000ff;
            hash[i + 28] = (ctx->state[7] >> (24 - i * 8)) & 0x000000ff;
        }
}

void
sha256(nsword8 hash[], nsword8 data[], size_t len)
/* SHA256 for a array of len bytes */
{
    SHA256_CTX ctx;

    sha256_init(&ctx);
    sha256_update(&ctx, data, len);
    sha256_final(&ctx,hash);
}

/**************************************************************************/

void
shahash(graph *g, int m, int n, nsword8 hash[])
/* Make SHA256 for a dense graph */
{
    SHA256_CTX ctx;
    graph *gi;
    int i,j,nb,sh;
    setword w;
    nsword8 x;

    nb = (n + 7)/8;  /* Bytes per row */

    sha256_init(&ctx);
    for (i = 0; i < n; ++i)
    {
        gi = g + m*i;
        w = *gi;
        sh = WORDSIZE - 8;
        for (j = 0; j < nb; ++j)
        {
            if (sh < 0) { w = *(++gi); sh = WORDSIZE - 8; }
            x = (w >> sh) & 0xFFU;

            ctx.data[ctx.datalen] = x;
            ctx.datalen++;
            if (ctx.datalen == 64)
            {
                sha256_transform(&ctx, ctx.data);
                ctx.bitlen += 512;
                ctx.datalen = 0;
            }
            sh -= 8;
        }
    }
    sha256_final(&ctx,hash);
}

static void
sha256_int2(SHA256_CTX *ctx, int k)
/* Insert 2-byte integer */
{
    nsword8 x;
    int w;
    
    w = k;

    x = w & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }

    x = (w >> 8) & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }
}

static void
sha256_int4(SHA256_CTX *ctx, int k)
/* Insert 4-byte integer */
{
    nsword8 x;
    int w;
    
    w = k;

    x = w & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }

    x = (w >> 8) & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }

    x = (w >> 16) & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }

    x = (w >> 24) & 0xFFU;
    ctx->data[ctx->datalen] = x;
    ctx->datalen++;
    if (ctx->datalen == 64) {
        sha256_transform(ctx, ctx->data);
        ctx->bitlen += 512;
        ctx->datalen = 0;
    }
}

void
shahash_sg(sparsegraph *sg, nsword8 hash[])
/* Make SHA256 for a sparse graph. You need to sort the adjacency lists
 * first (sortlists_sg()) if you want this to be consistent. */
{
    SHA256_CTX ctx;
    int n,i,j,*e,*d,*ei;
    size_t *v;

    SG_VDE(sg,v,d,e);
    n = sg->nv;
    sha256_init(&ctx);

    if (n < 65535)
    {
	for (i = 0; i < n; ++i)
	{
            sha256_int2(&ctx,d[i]);
            ei = e + v[i];
            for (j = 0; j < d[i]; ++j) sha256_int2(&ctx,ei[j]);
	}
        sha256_final(&ctx,hash);
    }
    else
    {
        for (i = 0; i < n; ++i)
        {
            sha256_int4(&ctx,d[i]);
            ei = e + v[i];
            for (j = 0; j < d[i]; ++j) sha256_int4(&ctx,ei[j]);
        }
        sha256_final(&ctx,hash);
    }
}
