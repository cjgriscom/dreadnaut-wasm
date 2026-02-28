/* uniqg.c version 1.2; B D McKay, August 2025 */

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
    -H  Write hash codes (32 bytes), not graphs (note: binary output)\n\
          Note that the output depends on the endianness of the hardware,\n\
          so -H and -h can only be used together on compatible endianness.\n\
    -J  Write hash codes without newlines, followed by newlines.\n\
          This output cannot be read with -h.\n\
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
    -y  Write a 256-bit cryptographic hashcode to stderr. This depends on\n\
         the set of isomorphism classes but not their order. It also\n\
         depends on -i, -I, -K, -S, -t and -c.\n\
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
#include "nausha.h"

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

typedef nsword64 hashcode[4];
typedef struct record {
   struct record *left,*right;
   hashcode hash;
} Record;

#define NUMTREES 64
#define TREEMSK ((nsword64)(0x3f))  /* NUMTREES bits */
static Record *root[NUMTREES]={NULL};
static Record *freelist=NULL;

#define ALLOCBATCH 128  /* How many to allocate at once */
/* Records in the free list are linked forwards using the left field */

static Record*
allocrecord(void)
{
    Record *p;
    int i;

    if (freelist == NULL)
    {
        if ((p = malloc(ALLOCBATCH*sizeof(Record))) == NULL)
            gt_abort(">E uniqg ran out of memory\n");
        for (i = 0; i < ALLOCBATCH-1; ++i) p[i].left = &p[i+1];
        p[ALLOCBATCH-1].left = NULL;
        freelist = p;
    }

    p = freelist;
    freelist = p->left;
    return p;
}

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

static void
tree_insert(hashcode code, boolean *isnew)
{
    Record *p,*pp;
    boolean left;
    int comp,i;

    p = root[code[2]&TREEMSK];
    pp = NULL;

    while (p)
    {
        comp = hashcompare(code,p->hash);
        if (comp == 0)
        {
            *isnew = FALSE;
            return;
        }
        else if (comp < 0)
        {
            pp = p;
            p = p->left;
            left = TRUE;
        }
        else
        {
            pp = p;
            p = p->right;
            left = FALSE;
        }
    }

    p = allocrecord();
    for (i = 0; i < 4; ++i) p->hash[i] = code[i];
    p->left = p->right = NULL;
    if (pp == NULL) root[code[2]&TREEMSK] = p;
    else if (left)  pp->left = p;
    else            pp->right = p;
    *isnew = TRUE;
}

/**************************************************************************/

#define ALREADY   1     /* Already canonical, needs 2, 4 or 8 as well. */
#define USEDENSE  2
#define USESPARSE 4
#define USETRACES 8

static boolean
processone(FILE *f, int prog, graph **gcan, sparsegraph *sgcan,
   boolean *digraph, int *nv, boolean *isnew, hashcode hash)
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
    graph *g=NULL;
    TracesStats traces_stats;
    statsblk nauty_stats;

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
        shahash(*gcan,m,n,(nsword8*)hash);
    else
        shahash_sg(sgcan,(nsword8*)hash);

    tree_insert(hash,isnew);
  
    if ((prog & USEDENSE) && !(prog & ALREADY))
    {
        if (*gcan == g) *gcan = NULL; 
        free(g);
    }

    return TRUE;
}

/**************************************************************************/

static void
makeystring(hashcode hash, char *ystring)
{
    static char c[] = "0123456789" "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                   "abcdefghijklmnopqrstuvwxyz";
    int i,ii,j;
    nsword64 wk;

    j = 0;
    for (i = 0; i < 4; ++i)
    {
        if (i > 0) ystring[j++] = ' ';
        wk = hash[i];
        for (ii = 0; ii < 11; ++ii)
        {
            ystring[j++] = c[wk % 62];
            wk /= 62;
        }
    }
    ystring[j] = '\0';
}

/**************************************************************************/

static void
putlinecode(FILE *f, hashcode hash)
/* Write a 32-byte code without newline, followed by a newline */
{
    unsigned char c,*s,buff[33];
    int i;
    nsword64 w;

    s = buff;
    for (i = 0; i < 4; ++i)
    {
        w = hash[i];
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF; w >>= 8;
        *(s++) = (c == '\n' ? 'X' : c);
        c = w & 0xFF;
        *(s++) = (c == '\n' ? 'X' : c);
    }
    *s = '\n';

    if (fwrite(buff,1,33,f) != 33)
        gt_abort(">E uniqg: fwrite failed in putlinecode()\n");
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    int m,n,codetype;
    int argnum,j;
    char *arg,sw;
    boolean quiet,badargs,digraph,isnew,yswitch;
    boolean fswitch,hswitch,Hswitch,kswitch,Jswitch;
    boolean iswitch,Iswitch,Kswitch,Sswitch,cswitch;
    boolean uswitch,tswitch,xswitch,Xswitch,Fswitch;
    char *xarg,*Xarg,*harg;
    long minil,maxil;
    double t;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    unsigned long long nin,nout;
    int prog;
    hashcode hash,cumhash;
    char ystring[50];
    graph *h;
    size_t nr;
    int inv,mininvarlevel,maxinvarlevel,invararg;
    SG_DECL(sh);

    HELP; PUTVERSION;

    if (sizeof(nsword32) != 4) gt_abort(">E nsword32 has wrong length\n");
    if (sizeof(nsword64) != 8) gt_abort(">E nsword64 has wrong length\n");

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    quiet = badargs = Hswitch = kswitch = cswitch = FALSE;
    fswitch = xswitch = Xswitch = hswitch = Jswitch = FALSE;
    iswitch = Iswitch = Kswitch = Fswitch = FALSE;
    uswitch = Sswitch = tswitch = yswitch = FALSE;
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
                else SWBOOLEAN('J',Jswitch)
                else SWBOOLEAN('F',Fswitch)
                else SWBOOLEAN('y',yswitch)
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
    if ((Hswitch != 0) + (Jswitch != 0) + (kswitch != 0) > 1)
        gt_abort(">E uniqg: -k, -H and -J are incompatible\n");

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
                 || Jswitch || kswitch || cswitch || Fswitch)
            fprintf(stderr," -");
        if (Sswitch) fprintf(stderr,"S");
        if (tswitch) fprintf(stderr,"t");
        if (Hswitch) fprintf(stderr,"H");
        if (Jswitch) fprintf(stderr,"J");
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
        while (processone(infile,prog,&h,&sh,&digraph,&n,&isnew,hash))
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
        while (processone(infile,prog|ALREADY,&h,&sh,&digraph,&n,&isnew,hash))
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
            tree_insert(hash,&isnew);
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
    cumhash[0] = cumhash[1] = cumhash[2] = cumhash[3] = 0;

    while (processone(infile,(cswitch?(prog|ALREADY):prog),&h,&sh,
                      &digraph,&n,&isnew,hash))
    {
        ++nin;
        if (isnew)
        {
            ++nout;
            cumhash[0] += hash[0]; cumhash[1] += hash[1];
            cumhash[2] += hash[2]; cumhash[3] += hash[3];
            if (uswitch)
            {
            }
            else if (Hswitch)
            {
                if (fwrite(hash,1,32,outfile) != 32)
                    gt_abort(">E error in writing hashcode\n");
            }
            else if (Jswitch)
                putlinecode(outfile,hash);
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

    if (yswitch)
    {
        makeystring(cumhash,ystring);
        fprintf(stderr,">Y %s\n",ystring);
    }

    if (!quiet)
    {
        if (uswitch)
            fprintf(stderr,
                ">Z %llu graphs read from %s, %llu unique; %.2f sec.\n",
                nin,infilename,nout,t);
        else if (Hswitch || Jswitch)
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
