/* distgraphg.c  version 1.2; B D McKay, Jan 2016. */

#define USAGE "distgraphg [-e|o|i] [-d#|-d#:#]... [-q] [infile [outfile]]"

#define HELPTEXT \
" Form graphs defined by distances in the input graphs. The -d,-e,-o-i\n\
  parameters define a set of distances. The output graph has an edge\n\
  if the distance in the input graph is in the set.\n\
\n\
  -d# -d#:# Specify a range of allowed distances.\n\
            Can be repeated up to 100 times.\n\
  -i Include infinity (unreachable vertices)\n\
  -e, -o  Only allow even or odd distances\n\
\n\
  If none of -d,-e,-o,-i is present, the square of the graph is taken.\n\
  If -e or -o appear without -d, all even or all odd distances are allowed.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

#define MAXDIST 100   /* Max number of distance ranges */
typedef struct {long lo,hi;} distrange;

/**************************************************************************/

static void
makedistances_sg(sparsegraph *sg, int u, int *dist)
/* Set dist[i] = distance from v to i.  dist[i]=n means infinity. */
{
    int i,j,n,head,tail,w;
    int *d,*e,*ew;
    size_t *v;
    DYNALLSTAT(int,queue,queue_sz);

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    DYNALLOC1(int,queue,queue_sz,n,"cagesearch");

    for (i = 0; i < n; ++i) dist[i] = n;
    queue[0] = u;
    dist[u] = 0;

    head = 0;
    tail = 1;
    while (tail < n && head < tail)
    {
        w = queue[head++];
        ew = e + v[w];
        for (j = 0; j < d[w]; ++j)
        {
            i = ew[j];
            if (dist[i] == n)
            {
                dist[i] = dist[w] + 1;
                queue[tail++] = i;
            }
        }
    }
}

/**************************************************************************/

static void
okdistances(distrange *dist, int numdist, int n, boolean infinite,
                boolean evenonly, boolean oddonly, int*distok)
/* Make vector of allowed distances, with n being infinity. */
{
    int i;
    long j,j0,j1;

    distok[0] = FALSE;
    if ((evenonly || oddonly) && numdist == 0 && !infinite)
        for (i = 1; i < n; ++i) distok[i] = TRUE;
    else
        for (i = 1; i < n; ++i) distok[i] = FALSE;
    distok[n] = infinite;

    for (i = 0; i < numdist; ++i)
    {
        j0 = (dist[i].lo < 1 ? 1 : dist[i].lo);
        j1 = (dist[i].hi == NOLIMIT || dist[i].hi >= n ? n-1 : dist[i].hi);
        for (j = j0; j <= j1; ++j) distok[j] = TRUE;
    }
    if (evenonly) for (i = 1; i < n; i += 2) distok[i] = FALSE;
    if (oddonly) for (i = 2; i < n; i += 2) distok[i] = FALSE;
}

/**************************************************************************/

static void
distgraph(sparsegraph *g, boolean *distok, sparsegraph *h)
/* h := distance graph of g */
{
    int *he,*hd;
    size_t *hv;
    int n,d,jj;
    size_t i,j;
    DYNALLSTAT(int,dist,dist_sz);

    n = g->nv;
    DYNALLOC1(int,dist,dist_sz,n,"distgraph");

    SG_ALLOC(*h,n,2L*n,"distgraph");
    SG_VDE(h,hv,hd,he);
    h->nv = n;
    j = 0;

    for (i = 0; i < n; ++i)
    {
        hv[i] = j;
        d = 0;
        if (h->elen < j + n)
        {
            DYNREALLOC(int,h->e,h->elen,j+2L*n,"distgraph");
            he = h->e;
        }
        makedistances_sg(g,i,dist);
        for (jj = 0; jj < n; ++jj)
            if (distok[dist[jj]]) { he[j++] = jj; ++d; }
        hd[i] = d;
    }
    h->nde = j;

    sortlists_sg(h);
}

/**************************************************************************/

static void
square(sparsegraph *g, sparsegraph *h)
/* h := square of g */
{
    DYNALLSTAT(int,mark,mark_sz);   /* vertex marker */
    int *ge,*gd,*he,*hd;
    size_t *gv,*hv;
    int n,jj,kk;
    size_t i,j,k,hnde,vi;

    SG_VDE(g,gv,gd,ge);
    n = g->nv;
    DYNALLOC1(int,mark,mark_sz,n,"distgraph");
    for (i = 0; i < n; ++i) mark[i] = -1;

    hnde = 0;
    for (i = 0; i < n; ++i) hnde += gd[i]*(size_t)gd[i];
    if (hnde > n*(size_t)(n-1)) hnde = n*(size_t)(n-1);
    SG_ALLOC(*h,n,hnde+1,"distgraph");
    SG_VDE(h,hv,hd,he);
    h->nv = n;

    vi = 0;
    for (i = 0; i < n; ++i)
    {
        mark[i] = i;
        hv[i] = vi;
        for (j = 0; j < gd[i]; ++j)
        {
            k = ge[gv[i]+j];
            if (mark[k] != i)
            {
                he[vi++] = k;
                mark[k] = i;
            }
            for (jj = 0; jj < gd[k]; ++jj)
            {
                kk = ge[gv[k]+jj];
                if (mark[kk] != i)
                {
                    he[vi++] = kk;
                    mark[kk] = i;
                }
            }
        }
        hd[i] = vi - hv[i];
    }
    h->nde = vi;

    sortlists_sg(h);
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet,dosquare;
    boolean infinite,evenonly,oddonly,dswitch;
    distrange dist[MAXDIST];
    int numdist,lastn,n;
    int i,j,argnum;
    int codetype,outcode,loops;
    boolean digraph;
    sparsegraph g,h;
    nauty_counter nin;
    char *arg,sw;
    DYNALLSTAT(boolean,distok,distok_sz);
    double t;

    HELP; PUTVERSION;

    SG_INIT(g);
    SG_INIT(h);

    infilename = outfilename = NULL;
    quiet = infinite = evenonly = oddonly = dswitch = FALSE;
    numdist = 0;

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
                     SWBOOLEAN('i',infinite)
                else SWBOOLEAN('e',evenonly)
                else SWBOOLEAN('o',oddonly)
                else SWBOOLEAN('q',quiet)
                else if (sw == 'd')
                {
                    if (numdist == MAXDIST) gt_abort_1(
                        ">E distgraphg: max distance ranges is %d\n",MAXDIST);
                    SWRANGE('d',":-",dswitch,dist[numdist].lo,
                                       dist[numdist].hi,"distgraphg -d");
                    ++numdist;
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

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (evenonly && oddonly)
        gt_abort(">E distgraphg: -e and -o are incompatible\n");

    dosquare = !evenonly && !oddonly && !infinite &&
       (numdist == 0 || (numdist == 1 && dist[0].lo == 1 && dist[0].hi == 2));

    if (!quiet)
    {
        fprintf(stderr,">A distgraphg");
        if (evenonly || oddonly || infinite)
            fprintf(stderr," -%s%s%s",(evenonly?"e":""),
                    (oddonly?"o":""),(infinite?"i":""));
        for (i = 0; i < numdist; ++i)
            fprintf(stderr," -d%ld:%ld",(dist[i].lo<=0?0:dist[i].lo),
                (dist[i].hi==NOLIMIT?0:dist[i].hi));
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
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

    if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    if (codetype&SPARSE6)  outcode = SPARSE6;
    else                   outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else if (outcode == GRAPH6) writeline(outfile,GRAPH6_HEADER);
        else  writeline(outfile,DIGRAPH6_HEADER);
    }

    nin = 0;
    lastn = -1;
    t = CPUTIME;
    while (read_sgg_loops(infile,&g,&loops,&digraph))
    {
        ++nin;

        if (dosquare)
            square(&g,&h);
        else
        {
            n = g.nv;
            if (n != lastn)
            {
                DYNALLOC1(boolean,distok,distok_sz,n+1,"distgraphg");
                okdistances(dist,numdist,n,infinite,
                                        evenonly,oddonly,distok);
                lastn = n;
            }
            distgraph(&g,distok,&h);
        }

        if (digraph) writed6_sg(outfile,&h);
        else if (outcode == SPARSE6) writes6_sg(outfile,&h);
        else writeg6_sg(outfile,&h);
    }
    t = CPUTIME - t;

    if (!quiet)
    {
        fprintf(stderr,">Z " COUNTER_FMT
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);
    }

    exit(0);
}
