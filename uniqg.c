/* uniqg.c version 0.9; B D McKay, April 2024 */

#define USAGE "uniqg [-q] [-xFILE] [-Xfile] [-hFILE] [-fxxx] [-u|-S|-t] \n\
                      [-c] [-k] [-i# -I#:# -K#] [infile [outfile]]"

#define HELPTEXT \
" Remove duplicates from a file of graphs or digraphs.\n\
  The SHA256 cryptographic hash function is used for comparisons\n\
\n\
    -S  Use sparse representation internally.\n\
         Note that this changes the canonical labelling.\n\
         Multiple edges are not supported.  One loop per vertex is ok.\n\
    -t  Use Traces.\n\
\n\
    -u  No output, just count\n\
    -H  Write hash codes, not graphs (note: binary output)\n\
    -k  Write the input graph exactly, not a canonical graph\n\
    -c  Assume graphs from infile are canonically labelled already\n\
\n\
    -xFILE  Read a file of graphs and exclude them from the output\n\
    -XFILE  Like -xFILE but assume they are already canonically labelled\n\
    -hFILE  Read a file of hash codes and exclude them from the output\n\
    -F  Flush output for each new graph (expensive if there are many)\n\
\n\
    -fxxx  Specify a partition of the vertex set.  xxx is any\n\
       string of ASCII characters except nul.  This string is\n\
       considered extended to infinity on the right with the\n\
       character 'z'. The sequence 'x^N', where x is a character and N is\n\
       a number, is equivalent to writing 'x' N times.  One character is\n\
       associated with each vertex, in the order given.  The labelling\n\
       used obeys these rules:\n\
        (1) the new order of the vertices is such that the associated\n\
       characters are in ASCII ascending order\n\
        (2) if two graphs are labelled using the same string xxx,\n\
       the output graphs are identical iff there is an\n\
       associated-character-preserving isomorphism between them.\n\
       If a leading '-' is used, as in -f-xxx, the characters are\n\
       assigned to the vertices starting at the last vertex, and\n\
       the new order of the vertices respects decreasing ASCII order.\n\
\n\
    -i#  select an invariant (1 = twopaths, 2 = adjtriang(K), 3 = triples,\n\
        4 = quadruples, 5 = celltrips, 6 = cellquads, 7 = cellquins,\n\
        8 = distances(K), 9 = indsets(K), 10 = cliques(K), 11 = cellcliq(K),\n\
       12 = cellind(K), 13 = adjacencies, 14 = cellfano, 15 = cellfano2,\n\
       16 = refinvar(K))\n\
    -I#:#  select mininvarlevel and maxinvarlevel (default 1:1)\n\
    -K#   select invararg (default 3)\n\
\n\
    -q  suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h"
#include "nautinv.h"
#include "gutils.h"
#include "traces.h"

/****************************************************************************
* This implementation of SHA256 is by Brad Conte, released to public domain.
* Author:     Brad Conte (brad AT bradconte.com) */

#if HAVE_STDINT_H
typedef unsigned char word8;
typedef uint32_t word32;
typedef uint64_t word64;
#else
typedef unsigned char word8;
typedef unsigned int word32;  /* Must be 32 bits */
typedef unsigned long long word64;  /* Must be 64 bits */
#endif

typedef struct {
    word8 data[64];
    unsigned long long bitlen;
    word32 datalen;
    word32 state[8];
} SHA256_CTX;

void sha256_init(SHA256_CTX *ctx);
void sha256_update(SHA256_CTX *ctx, const word8 data[], size_t len);
void sha256_final(SHA256_CTX *ctx, word8 hash[]);
void sha256(word8 hash[], word8 data[], size_t len);

#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

static const word32 k[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

void sha256_transform(SHA256_CTX *ctx, const word8 data[])
{
    word32 a, b, c, d, e, f, g, h, i, j, t1, t2, m[64];

    for (i = 0, j = 0; i < 16; ++i, j += 4)
            m[i] = (data[j] << 24) | (data[j + 1] << 16) | (data[j + 2] << 8) | (data[j + 3]);
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

void sha256_init(SHA256_CTX *ctx)
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

void sha256_update(SHA256_CTX *ctx, const word8 data[], size_t len)
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

void sha256_update_small(SHA256_CTX *ctx, const int data[], size_t len)
/* Version for values < 2^16, ignoring high order zeros */
{
    size_t i;
    word8 x;

    for (i = 0; i < len; ++i) {
        x = data[i] & 0xFF;
        ctx->data[ctx->datalen] = x;
        ctx->datalen++;
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
        x = ((word32)data[i] >> 8) & 0xFF;
        ctx->data[ctx->datalen] = x;
        ctx->datalen++;
        if (ctx->datalen == 64) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

void sha256_final(SHA256_CTX *ctx, word8 hash[])
{
    word32 i;

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

    // Since this implementation uses little endian byte ordering and SHA uses big endian,
    // reverse all the bytes when copying the final state to the output hash.
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

void sha256(word8 hash[],word8 data[],size_t len)
{
    SHA256_CTX ctx;

    sha256_init(&ctx);
    sha256_update(&ctx, data, len);
    sha256_final(&ctx,hash);
}

static struct invarrec
{
    void (*entrypoint)(graph*,int*,int*,int,int,int,int*,
                      int,boolean,int,int);
    void (*entrypoint_sg)(graph*,int*,int*,int,int,int,int*,
                      int,boolean,int,int);
    char *name;
} invarproc[]
    = {{NULL, NULL, "none"},
       {twopaths,   NULL, "twopaths"},
       {adjtriang,  NULL, "adjtriang"},
       {triples,    NULL, "triples"},
       {quadruples, NULL, "quadruples"},
       {celltrips,  NULL, "celltrips"},
       {cellquads,  NULL, "cellquads"},
       {cellquins,  NULL, "cellquins"},
       {distances, distances_sg, "distances"},
       {indsets,    NULL, "indsets"},
       {cliques,    NULL, "cliques"},
       {cellcliq,   NULL, "cellcliq"},
       {cellind,    NULL, "cellind"},
       {adjacencies, adjacencies_sg, "adjacencies"},
       {cellfano,   NULL, "cellfano"},
       {cellfano2,  NULL, "cellfano2"},
       {refinvar,   NULL, "refinvar"}
      };

#define NUMINVARS ((int)(sizeof(invarproc)/sizeof(struct invarrec)))

static DEFAULTOPTIONS_GRAPH(dense_options);
static DEFAULTOPTIONS_SPARSEGRAPH(sparse_options);
static DEFAULTOPTIONS_DIGRAPH(digraph_options);
static DEFAULTOPTIONS_TRACES(traces_options);
static char *fmt;

typedef word64 hashcode[4];
typedef struct splaynode {
   struct splaynode *parent,*left,*right;
   hashcode hash;
} SPLAYNODE;

static SPLAYNODE *root=NULL;
static SPLAYNODE *freelist=NULL;

#define SPLAYALLOC 64  /* How many to allocate at once */
/* SPLAYNODEs in the free list are linked forwards using the parent field */

static SPLAYNODE*
allocsplaynode(size_t sz)
{
    SPLAYNODE *p;
    int i;

    if (freelist == NULL)
    {
        if ((p = malloc(SPLAYALLOC*sizeof(SPLAYNODE))) == NULL)
            gt_abort(">E uniqg ran out of memory\n");
        for (i = 0; i < SPLAYALLOC-1; ++i) p[i].parent = &p[i+1];
        p[SPLAYALLOC-1].parent = NULL;
        freelist = p;
    }

    p = freelist;
    freelist = p->parent;
    return p;
}

#define SPLAY_ALLOC allocsplaynode
#define SPLAY_FREE(pp) { pp->parent = freelist; freelist = pp; }

static int
hashcompare(hashcode hash1, hashcode hash2)
{
    if (hash1[0] < hash2[0]) return -1;
    if (hash1[0] > hash2[0]) return 1;
    if (hash1[1] < hash2[1]) return -1;
    if (hash1[1] > hash2[1]) return 1;
    if (hash1[2] < hash2[2]) return -1;
    if (hash1[2] > hash2[2]) return 1;
    if (hash1[3] < hash2[3]) return -1;
    if (hash1[3] > hash2[3]) return 1;
    return 0;
}

#define ACTION(p) {}
#define SCAN_ARGS
#define INSERT_ARGS , hashcode hash, boolean *isnew
#define PRESENT(p) *isnew = FALSE;
#define NOT_PRESENT(p) { *isnew = TRUE; p->hash[0] = hash[0]; \
   p->hash[1] = hash[1]; p->hash[2] = hash[2]; p->hash[3] = hash[3]; }
#define COMPARE(p) hashcompare(hash,p->hash)

#include "splay.c"

/**************************************************************************/

static void
shahash(graph *g, int m, int n, hashcode hash)
{
    sha256((unsigned char*)hash,(unsigned char*)g,m*(size_t)n*sizeof(graph));
}

static void
shahash_sg(sparsegraph *sg, hashcode hash)
{
    SHA256_CTX ctx;
    int n,i,*e,*d;
    size_t *v;

    SG_VDE(sg,v,d,e);
    n = sg->nv;
    sha256_init(&ctx);
    for (i = 0; i < n; ++i)
    {
        sha256_update(&ctx,(word8*)(d+i),sizeof(int));
        if (d[i] > 0) sha256_update(&ctx,(word8*)(e+v[i]),d[i]*sizeof(int));
        sha256_final(&ctx,(word8*)hash);
    }
}

static void
shahash_small_sg(sparsegraph *sg, hashcode hash)
{
    SHA256_CTX ctx;
    int n,i,*e,*d;
    size_t *v;

    SG_VDE(sg,v,d,e);
    n = sg->nv;
    sha256_init(&ctx);
    for (i = 0; i < n; ++i)
    {
        sha256_update_small(&ctx,d+i,1);
        if (d[i] > 0) sha256_update_small(&ctx,e+v[i],d[i]);
        sha256_final(&ctx,(word8*)hash);
    }
}

/**************************************************************************/

#define ALREADY   1     /* Already canonical, needs 2, 4 or 8 as well. */
#define USEDENSE  2
#define USESPARSE 4
#define USETRACES 8

static boolean
processone(FILE *f, int prog, graph **gcan, sparsegraph *sgcan,
   boolean *digraph, int *nv, boolean *isnew)
/* Read from f, return FALSE if EOF. Otherwise, set *digraph, *n.
   Label with prog and insert hash into tree. Set *isnew.
   If writecanon, set either gcan or sgcan to canonical form. */
{
    int m,n,loops;
    static SG_DECL(sg);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(graph,h,h_sz);
    DYNALLSTAT(setword,work,work_sz);
    DYNALLSTAT(set,active,active_sz);
    graph *g;
    TracesStats traces_stats;
    statsblk nauty_stats;
    hashcode hash;

    if ((prog & USEDENSE))
    {
        if ((g = readg_loops(f,NULL,0,&m,&n,&loops,digraph)) == NULL)
            return FALSE;
        *nv = n;
    }
    else
    {
        if (read_sgg_loops(f,&sg,&loops,digraph) == NULL) return FALSE;
        if ((*digraph || loops > 0) && (prog & USETRACES))
            gt_abort(">E Traces cannot handle digraphs or loops\n"); 
        *nv = n = sg.nv;
    }
    m = SETWORDSNEEDED(n);

    if ((prog & ALREADY))
    {
        if ((prog & USEDENSE)) *gcan = g;
        else 
        {
            SG_TRANSFER(*sgcan,sg);
            sortlists_sg(sgcan);
        }
    }    
    else
    {
        DYNALLOC1(int,lab,lab_sz,n,"uniqg malloc"); 
        DYNALLOC1(int,ptn,ptn_sz,n,"uniqg malloc"); 
        DYNALLOC1(int,orbits,orbits_sz,n,"uniqg malloc"); 
        if (!(prog & USETRACES))
        {
            DYNALLOC1(setword,work,work_sz,1000*m,"uniqg malloc");
            DYNALLOC1(set,active,active_sz,m,"uniqg malloc");
        }

        if ((prog & USEDENSE))
        {
            DYNALLOC2(graph,h,h_sz,m,n,"uniqg malloc");
            dense_options.digraph = (loops>0||*digraph);
            setlabptnfmt(fmt,lab,ptn,active,m,n);
            nauty(g,lab,ptn,active,orbits,&dense_options,&nauty_stats,
                        work,1000*m,m,n,h);
            *gcan = h;
        }
        else if ((prog & USESPARSE))
        {
            SG_ALLOC(*sgcan,n,sg.nde,"uniqg malloc");
            sparse_options.digraph = (loops>0||*digraph);
            setlabptnfmt(fmt,lab,ptn,active,m,n);
            nauty((graph*)(&sg),lab,ptn,active,orbits,&sparse_options,
                        &nauty_stats,work,1000*m,m,n,(graph*)(sgcan));
            sortlists_sg(sgcan);
        }
        else   /* Traces */
        {
            SG_ALLOC(*sgcan,n,sg.nde,"uniqg");
            sgcan->nv = n;
            sgcan->nde = sg.nde;
            setlabptnfmt(fmt,lab,ptn,NULL,0,n);
            Traces(&sg,lab,ptn,orbits,&traces_options,&traces_stats,sgcan);
            sortlists_sg(sgcan);
        }
    }

    if ((prog & USEDENSE))
        shahash(*gcan,m,n,hash);
    else if ((long)n < 65535)
        shahash_small_sg(sgcan,hash);
    else
        shahash_sg(sgcan,hash);

    splay_insert(&root,hash,isnew);
  
    if ((prog & USEDENSE) && !(prog & ALREADY))
    {
        if (*gcan == g) *gcan = NULL; 
        free(g);
    }

    return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    int m,n,codetype;
    int argnum,j;
    char *arg,sw;
    boolean quiet,badargs,digraph,isnew;
    boolean fswitch,hswitch,Hswitch,kswitch;
    boolean iswitch,Iswitch,Kswitch,Sswitch,cswitch;
    boolean uswitch,tswitch,xswitch,Xswitch,Fswitch;
    char *xarg,*Xarg,*harg;
    long minil,maxil;
    double t;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    unsigned long long nin,nout;
    int prog;
    hashcode hash;
    graph *h;
    size_t nr;
    int inv,mininvarlevel,maxinvarlevel,invararg;
    SG_DECL(sh);

    HELP; PUTVERSION;

    if (sizeof(word32) != 4) gt_abort(">E word32 has wrong length\n");
    if (sizeof(word64) != 8) gt_abort(">E word64 has wrong length\n");

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    quiet = badargs = Hswitch = kswitch = cswitch = FALSE;
    fswitch = xswitch = Xswitch = hswitch = FALSE;
    iswitch = Iswitch = Kswitch = Fswitch = FALSE;
    uswitch = Sswitch = tswitch = FALSE;
    infilename = outfilename = NULL;
    inv = 0;

    argnum = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('S',Sswitch)
                else SWBOOLEAN('t',tswitch)
                else SWBOOLEAN('k',kswitch)
                else SWBOOLEAN('c',cswitch)
                else SWBOOLEAN('H',Hswitch)
                else SWBOOLEAN('F',Fswitch)
                else SWINT('i',iswitch,inv,"uniqg -i")
                else SWINT('K',Kswitch,invararg,"uniqg -K")
                else SWRANGE('k',":-",Iswitch,minil,maxil,"uniqg -k")
                else SWRANGE('I',":-",Iswitch,minil,maxil,"uniqg -I")
                else if (sw == 'f')
                {
                    fswitch = TRUE;
                    fmt = arg;
                    break;
                }
                else if (sw == 'x')
                {
                    xswitch = TRUE;
                    xarg = arg;
                    break;
                }
                else if (sw == 'X')
                {
                    Xswitch = TRUE;
                    Xarg = arg;
                    break;
                }
                else if (sw == 'h')
                {
                    hswitch = TRUE;
                    harg = arg;
                    break;
                }
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1) infilename = arg;
            else if (argnum == 2) outfilename = arg;
            else                  badargs = TRUE;
        }
    }

    if (tswitch && Sswitch)
        gt_abort(">E uniqg: -t and -S are incompatible\n");

    if (iswitch && inv == 0) iswitch = FALSE;

    if (!Sswitch && iswitch && (inv > NUMINVARS))
        gt_abort(">E uniqg: -i value must be 0..16\n");
    if (tswitch && iswitch)
        gt_abort(">E uniqg: invariants are not available with -t\n");
    if (Sswitch && iswitch && invarproc[inv].entrypoint_sg == NULL)
        gt_abort(
            ">E uniqg: that invariant is not available in sparse mode\n");

    if (iswitch)
    {
        if (Iswitch)
        {
            mininvarlevel = minil;
            maxinvarlevel = maxil;
        }
        else
            mininvarlevel = maxinvarlevel = 1;
        if (!Kswitch) invararg = 3;
    }
    if (!fswitch) fmt = NULL;

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (Sswitch)      prog = USESPARSE;
    else if (tswitch) prog = USETRACES;
    else              prog = USEDENSE;

    if (!quiet)
    {
        fprintf(stderr,">A uniqg");
        if (fswitch || iswitch || tswitch || Sswitch || Hswitch
                 || kswitch || cswitch || Fswitch)
            fprintf(stderr," -");
        if (Sswitch) fprintf(stderr,"S");
        if (tswitch) fprintf(stderr,"t");
        if (Hswitch) fprintf(stderr,"H");
        if (kswitch) fprintf(stderr,"k");
        if (cswitch) fprintf(stderr,"c");
        if (Fswitch) fprintf(stderr,"F");
        if (iswitch)
            fprintf(stderr,"i=%s[%d:%d,%d]",invarproc[inv].name,
                    mininvarlevel,maxinvarlevel,invararg);
        if (fswitch) fprintf(stderr," -f%s",fmt);
        if (hswitch) fprintf(stderr," -h%s",harg);
        if (xswitch) fprintf(stderr," -x%s",xarg);
        if (Xswitch) fprintf(stderr," -X%s",Xarg);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    root = NULL;

    dense_options.getcanon = sparse_options.getcanon = TRUE;
    digraph_options.getcanon = traces_options.getcanon = TRUE;
    dense_options.defaultptn = sparse_options.defaultptn = FALSE;
    digraph_options.defaultptn = traces_options.defaultptn = FALSE;
    traces_options.verbosity = 0;
    if (iswitch)
    {
        dense_options.invarproc = invarproc[inv].entrypoint;
        dense_options.mininvarlevel = mininvarlevel;
        dense_options.maxinvarlevel = maxinvarlevel;
        dense_options.invararg = invararg;
        sparse_options.invarproc = invarproc[inv].entrypoint_sg;
        sparse_options.mininvarlevel = mininvarlevel;
        sparse_options.maxinvarlevel = maxinvarlevel;
        sparse_options.invararg = invararg;
        digraph_options.invarproc = invarproc[inv].entrypoint;
        digraph_options.mininvarlevel = mininvarlevel;
        digraph_options.maxinvarlevel = maxinvarlevel;
        digraph_options.invararg = invararg;
    }

    t = CPUTIME;

    if (xswitch)
    {
        if ((infile = opengraphfile(xarg,&codetype,FALSE,1)) == NULL)
            exit(1);
        nin = nout = 0;
        while (processone(infile,prog,&h,&sh,&digraph,&n,&isnew))
        {
            ++nin;
            if (isnew) ++nout;
        }
        fclose(infile);
        if (!quiet) fprintf(stderr,
                        ">x %llu exclusions (%llu unique) read from %s\n",
                        nin,nout,xarg);
    }
 
    if (Xswitch)
    {
        if ((infile = opengraphfile(Xarg,&codetype,FALSE,1)) == NULL)
            exit(1);
        nin = nout = 0;
        while (processone(infile,prog|ALREADY,&h,&sh,&digraph,&n,&isnew))
        {
            ++nin;
            if (isnew) ++nout;
        }
        fclose(infile);
        if (!quiet) fprintf(stderr,
                   ">X %llu labelled exclusions (%llu unique) read from %s\n",
                   nin,nout,Xarg);
    }
 
    if (hswitch)
    {
        if ((infile = fopen(harg,"rb")) == NULL)
            gt_abort_1(">E uniqg: Can't open %s for reading\n",harg);
        nin = nout = 0;
        while ((nr = fread(hash,1,32,infile)) != 0)
        {
            if (nr != 32) gt_abort_1(">E error reading %s\n",harg);
            ++nin;
            splay_insert(&root,hash,&isnew);
            if (isnew) ++nout;
        }
        fclose(infile);
        if (!quiet) fprintf(stderr,
                   ">H %llu hash codes (%llu unique) read from %s\n",
                   nin,nout,harg);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);

    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    if (!outfilename || outfilename[0] == '-')
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,"w")) == NULL)
        gt_abort_1(">E Can't open output file %s\n",outfilename);

    if (!uswitch && !Hswitch && (codetype&HAS_HEADER))
    {
        if      ((codetype&SPARSE6))  writeline(outfile,SPARSE6_HEADER);
        else if ((codetype&DIGRAPH6)) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    nin = nout = 0;

    while (processone(infile,(cswitch?(prog|ALREADY):prog),&h,&sh,
                      &digraph,&n,&isnew))
    {
        ++nin;
        if (isnew)
        {
            ++nout;
            if (uswitch)
            {
            }
            else if (Hswitch)
            {
                if (fwrite(root->hash,1,32,outfile) != 32)
                    gt_abort(">E error in writing hashcode\n");
            }
            else if (kswitch)
                writelast(outfile);
            else if ((prog & USEDENSE))
            {
                m = SETWORDSNEEDED(n);
                if (readg_code == SPARSE6)
                    writes6(outfile,h,m,n);
                else if (readg_code == DIGRAPH6)
                    writed6(outfile,h,m,n);
                else
                    writeg6(outfile,h,m,n);
            }
            else
            {
                if (readg_code == SPARSE6)
                    writes6_sg(outfile,&sh);
                else if (readg_code == DIGRAPH6)
                    writed6_sg(outfile,&sh);
                else
                    writeg6_sg(outfile,&sh);
            }
        }
        if (Fswitch) fflush(outfile);
        if ((prog & ALREADY) && (prog & USEDENSE)) free(h);
    }

    t = CPUTIME - t;

    if (!quiet)
    {
        if (uswitch)
            fprintf(stderr,
                ">Z %llu graphs read from %s, %llu unique; %.2f sec.\n",
                nin,infilename,nout,t);
        else if (Hswitch)
            fprintf(stderr,
                ">Z %llu graphs read from %s, "
                "%llu hashcodes written to %s; %.2f sec.\n",
                nin,infilename,nout,outfilename,t);
        else
            fprintf(stderr,
                ">Z %llu graphs read from %s, %llu written to %s; %.2f sec.\n",
                nin,infilename,nout,outfilename,t);
    }

    exit(0);
}
