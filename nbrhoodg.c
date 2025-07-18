/* nbrhoodg.c  version 2.0; B D McKay, Feb 2025. */

#define USAGE "nbrhoodg [-lq] [-c|-C] [-d#|d#:#] [-v#|-v#:#] [infile [outfile]]"

#define HELPTEXT \
" Extract neighbourhoods of vertices of a graph.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    No isomorph reduction is done. No null graphs are written.\n\
\n\
    -l  Canonically label outputs (default is same labelling as input)\n\
    -C  Extract closed neighbourhoods instead.\n\
    -c  Extract non-neighbourhoods instead.\n\
    -D# Extract neighbourhoods out to distance # (implies -C)\n\
    -d# -d#:# Only include vertices with original degree in the given range\n\
    -v# -v#:# Only include vertices with these vertex numbers (first is 0).\n\
        No empty graphs are output.\n\
        For digraphs, out-degree and out-neighbourhoods are used.\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 

static FILE *outfile;
static nauty_counter nout;
static int outcode;
static boolean digraph,dolabel;

static int
list_dist(graph *g, int m, int n, int v, int *perm,
                long mindist, long maxdist, boolean compl)
/* Put in perm[0..] the vertices at distance mindist..maxdist from v.
   If compl=TRUE, use the complement of that set.
   Assumes maxdist >= mindist. Returns the number of vertices. */
{
    int i,head,tail,w,dis,nperm;
    set *gw;
#if MAXN
    int queue[MAXN],dist[MAXN];
#else
    DYNALLSTAT(int,queue,queue_sz);
    DYNALLSTAT(int,dist,dist_sz);
#endif

    if (n == 0) return 0;

#if !MAXN
    DYNALLOC1(int,queue,queue_sz,n,"nbrhoodg");
    DYNALLOC1(int,dist,dist_sz,n,"nbrhoodg");
#endif

    for (i = 0; i < n; ++i) dist[i] = NAUTY_INFINITY;

    queue[0] = v;
    dist[v] = dis = 0;

    head = 0;
    tail = 1;
    while (dis <= maxdist && tail < n && head < tail)
    {
        w = queue[head++];
        gw = GRAPHROW(g,w,m);
        for (i = -1; (i = nextelement(gw,m,i)) >= 0;)
        {
            if (dist[i] == NAUTY_INFINITY)
            {
                dis = dist[i] = dist[w] + 1;
                queue[tail++] = i;
            }
        }
    }

    nperm = 0;
    if (compl)
    {
        for (i = 0; i < n; ++i)
            if (dist[i] < mindist || dist[i] > maxdist) perm[nperm++] = i;
    }
    else
    {
        for (i = 0; i < n; ++i)
            if (dist[i] >= mindist && dist[i] <= maxdist) perm[nperm++] = i;
    }

    return nperm;
}

/**************************************************************************/

static void
getsubgraph(graph *gin, int *perm, int nperm, graph *gout, int m, int n)
/* Makes a subgraph gin<perm[0..nperm-1]> */
{
    int i,j,k;
    int newm;
    set *gi,*wgi;

    newm = SETWORDSNEEDED(nperm);

    EMPTYGRAPH(gout,newm,nperm);

    for (i = 0, gi = (set*)gout; i < nperm; ++i, gi += newm)
    {
        wgi = GRAPHROW(gin,perm[i],m);
        for (j = 0; j < nperm; ++j)
        {
            k = perm[j];
            if (ISELEMENT(wgi,k)) ADDELEMENT(gi,j);
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile;
    boolean badargs,quiet,dswitch,vswitch,cswitch,Cswitch,Dswitch;
    int i,j,m,n,v,argnum;
    int codetype;
    graph *g,*gq;
    nauty_counter nin;
    char *arg,sw;
    setword *gv;
    long mindist,maxdist;
    long mindeg,maxdeg;
    long minvert,maxvert;
    int degv,msub,nsub;
    double t;
#if MAXN
    boolean degok[MAXN];
    boolean del[MAXN];
#else
    DYNALLSTAT(boolean,perm,perm_sz);
    DYNALLSTAT(graph,gsub,gsub_sz);
    DYNALLSTAT(graph,gx,gx_sz);
#endif

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
    Cswitch = cswitch = vswitch = dswitch = Dswitch = quiet = FALSE;
    mindist = maxdist = 0;

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
                     SWBOOLEAN('l',dolabel)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('c',cswitch)
                else SWBOOLEAN('C',Cswitch)
                else SWRANGE('D',":-",Dswitch,mindist,maxdist,">E nbrhoodg -D")
                else SWRANGE('v',":-",vswitch,minvert,maxvert,">E nbrhoodg -v")
                else SWRANGE('d',":-",dswitch,mindeg,maxdeg,">E nbrhoodg -d")
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

    if (Cswitch && (cswitch || Dswitch))
        gt_abort(">E nbrhoodg: -C is incompatible with -c and -D\n");
        
    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A nbrhoodg");
        if (dolabel) fprintf(stderr," -l");
        if (dswitch) fprintf(stderr," -d%ld:%ld",mindeg,maxdeg);
        if (vswitch) fprintf(stderr," -v%ld:%ld",minvert,maxvert);
        if (Dswitch) fprintf(stderr," -D%ld:%ld",mindist,maxdist);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (!dswitch)
    {
        mindeg = 0;
        maxdeg = NAUTY_INFINITY;
    }
    if (!vswitch)
    {
        minvert = 0;
        maxvert = NAUTY_INFINITY;
    }
    if (mindeg < 0) mindeg = 0;
    if (maxdeg == NOLIMIT) maxdeg = NAUTY_INFINITY;
    if (minvert < 0) minvert = 0;
    if (maxvert == NOLIMIT) maxvert = NAUTY_INFINITY;
    if (mindist < 0) mindist = 0;
    if (maxdist == NOLIMIT) maxdist = NAUTY_INFINITY;

    if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = nout = 0;
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

    if (codetype&SPARSE6)       outcode = SPARSE6;
    else if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    else                        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

#if !MAXN
        DYNALLOC2(graph,gsub,gsub_sz,m,n,"nbrhoodg");
        if (dolabel) DYNALLOC2(graph,gx,gx_sz,m,n,"nbrhoodg");
        DYNALLOC1(int,perm,perm_sz,n,"nbrhoodg");
#endif

        for (v = minvert, gv = GRAPHROW(g,minvert,m);
                              v < n && v <= maxvert; ++v, gv += m)
        {
            if (digraph) degv = 0;
            else degv = (ISELEMENT(gv,v) != 0);
            for (i = 0; i < m; ++i) degv += POPCOUNT(gv[i]);
            if (degv < mindeg || degv > maxdeg) continue;

            if (Dswitch)
            {
                nsub = list_dist(g,m,n,v,perm,mindist,maxdist,cswitch);
            }
            else if (Cswitch)
            {
                perm[0] = v;
                nsub = 1;
                for (i = 0; i < n; ++i)
                    if (ISELEMENT(gv,i) && i != v) perm[nsub++] = i;
            }
            else if (cswitch)
            {
                nsub = 0;
                for (i = 0; i < n; ++i)
                    if (!ISELEMENT(gv,i) && i != v) perm[nsub++] = i;
            }
            else
            {
                nsub = 0;
                for (i = 0; i < n; ++i)
                    if (ISELEMENT(gv,i) && i != v) perm[nsub++] = i;
            }

            if (nsub != 0)
            {
                getsubgraph(g,perm,nsub,gsub,m,n);
                msub = SETWORDSNEEDED(nsub);
                if (dolabel)
                {
                    fcanonise(gsub,msub,nsub,gx,NULL,digraph);
                    gq = gx;
                }
                else
                    gq = gsub;

                if (outcode == DIGRAPH6 || digraph) writed6(outfile,gq,msub,nsub);
                else if (outcode == SPARSE6)        writes6(outfile,gq,msub,nsub);
                else                                writeg6(outfile,gq,msub,nsub);
                ++nout;
            }
        }

        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,
             ">Z  " COUNTER_FMT " graphs read from %s, "
                    COUNTER_FMT " written to %s; %3.2f sec.\n",
             nin,infilename,nout,outfilename,t);

    exit(0);
}
