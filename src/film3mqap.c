/* FILM3MQAP - estimate model quality from FILM3 ensemble of protein structures - by David Jones, April 2012 */

/* V0.2 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQLEN 5000
#define MAXNSTRUC 10000

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#define BIG (1.0e32)
#define SMALL (1.0e-5)
#define ZERO (0.0)
#define PI (3.14159265)

#define NTOKENS 26
#define DEFR 1.80

/* A list of common PDB record types... */

#define HEADER             1
#define COMPND             2
#define SOURCE             3
#define AUTHOR             4
#define REVDAT             5
#define REMARK             6
#define SEQRES             7
#define FTNOTE             8
#define HET                9
#define FORMUL             10
#define HELIX              11
#define CRYST1             12
#define ORIGX1             13
#define ORIGX2             14
#define ORIGX3             15
#define SCALE1             16
#define SCALE2             17
#define SCALE3             18
#define ATOM               19
#define TER                20
#define HETATM             21
#define CONECT             22
#define ENDENT             23
#define JRNL               24
#define TURN               25
#define ENDMDL             26

#define MEMGRAIN	   100

#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define veczero(v) memset(v, 0, sizeof(v))
#define veccopy(a,b) (a[0]=b[0],a[1]=b[1],a[2]=b[2])
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) (a[0]=b[0]*c,a[1]=b[1]*c,a[2]=b[2]*c)
#define ran0() ((rand()&32767)/(double)32768.0)
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define SQR(x) ((x)*(x))
#define dist(a,b) sqrt(SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))

typedef double   Transform[4][4];
typedef double   Point[3];

typedef struct
{
    short           length;
    short           posn_a[MAXSEQLEN], posn_b[MAXSEQLEN];
}
ALNS;

struct pdbatm
{
    Point           pos;
    double          radius;	/* Radius of sphere */
    int             resnum, aac, nequiv;
    char            ispolar, isback, donors, acceptors;
    char            ndon, nacc, donflag, accflag;
    char            atmnam[5];
}              *atoms[MAXNSTRUC];

char            pdbfn[80], csdfn[80], logfn[80], keyword[40], buf[8192];

/* Record names for decoding record types */
char           *tokstr[] =
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL",
 "HELIX", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM", "TER",
 "HETATM", "CONECT", "END", "JRNL", "TURN", "ENDMDL"
};

/* Residue name to allow conversion of a.a. name into numeric code */
char           *rnames[] =
{
 "ALA", "ARG", "ASN", "ASP", "CYS",
 "GLN", "GLU", "GLY", "HIS", "ILE",
 "LEU", "LYS", "MET", "PHE", "PRO",
 "SER", "THR", "TRP", "TYR", "VAL",
 "ASX", "GLX", "UNK"
};

enum RESCODE
{
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys,
    Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val,
    Asx, Glx, Unk
};

int             nmodels, natoms[MAXNSTRUC];

int pat[MAXSEQLEN][MAXSEQLEN];


/***************************************************************************/

void
                fail(errstr)
    char           *errstr;
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}


/* Apply a Transform matrix to a point */
void
                transform_point(transform, p, tp)
    Transform       transform;	/* transform to apply to the point */
    Point           p;		/* the point to transform */
    Point           tp;		/* the returned point after transformation */
{
    int             i, j;
    Point           temp;

    temp[0] = p[0] + transform[0][3];
    temp[1] = p[1] + transform[1][3];
    temp[2] = p[2] + transform[2][3];
    tp[0] = dotprod(transform[0], temp);
    tp[1] = dotprod(transform[1], temp);
    tp[2] = dotprod(transform[2], temp);

#if 0
    printf("%g %g %g\n%g %g %g\n\n", p[0], p[1], p[2], tp[0], tp[1], tp[2]);
#endif
}

/* Eigen decomposition code for symmetric 6x6 matrix, taken/modified from the public
   domain Java Matrix library JAMA. */

#define NDIM 6

#define hypot2(x, y) (sqrt((x)*(x)+(y)*(y)))

/* Symmetric Householder reduction to tridiagonal form. */

void tred2(double V[NDIM][NDIM], double d[NDIM], double e[NDIM]) 
{
    int i, j, k;
    
    for (j = 0; j < NDIM; j++)
	d[j] = V[NDIM-1][j];

    /* Householder reduction to tridiagonal form. */

    for (i = NDIM-1; i > 0; i--)
    {
	/* Scale to avoid under/overflow. */

	double scale = 0.0;
	double h = 0.0;
	
	for (k = 0; k < i; k++)
	    scale = scale + fabs(d[k]);
	
	if (scale == 0.0)
	{
	    e[i] = d[i-1];
	    
	    for (j = 0; j < i; j++)
	    {
		d[j] = V[i-1][j];
		V[i][j] = 0.0;
		V[j][i] = 0.0;
	    }
	}
	else
	{
	    /* Generate Householder vector. */

	    for (k = 0; k < i; k++)
	    {
		d[k] /= scale;
		h += d[k] * d[k];
	    }
	    
	    double f = d[i-1];
	    double g = sqrt(h);
	    
	    if (f > 0)
		g = -g;
	    
	    e[i] = scale * g;
	    h = h - f * g;
	    d[i-1] = f - g;
	    
	    for (j = 0; j < i; j++)
		e[j] = 0.0;

	    /* Apply similarity transformation to remaining columns. */

	    for (j = 0; j < i; j++)
	    {
		f = d[j];
		V[j][i] = f;
		g = e[j] + V[j][j] * f;

		for (k = j+1; k <= i-1; k++)
		{
		    g += V[k][j] * d[k];
		    e[k] += V[k][j] * f;
		}
		
		e[j] = g;
	    }
	    
	    f = 0.0;
	    
	    for (j = 0; j < i; j++)
	    {
		e[j] /= h;
		f += e[j] * d[j];
	    }
	    
	    double hh = f / (h + h);
	    
	    for (j = 0; j < i; j++)
		e[j] -= hh * d[j];
	    
	    for (j = 0; j < i; j++)
	    {
		f = d[j];
		g = e[j];
		
		for (k = j; k <= i-1; k++)
		    V[k][j] -= (f * e[k] + g * d[k]);
		
		d[j] = V[i-1][j];
		
		V[i][j] = 0.0;
	    }
	}
	
	d[i] = h;
    }
    

    /* Accumulate transformations. */

    for (i = 0; i < NDIM-1; i++)
    {
	V[NDIM-1][i] = V[i][i];
	
	V[i][i] = 1.0;
	
	double h = d[i+1];
	
	if (h != 0.0)
	{
	    for (k = 0; k <= i; k++)
		d[k] = V[k][i+1] / h;
	    
	    for (j = 0; j <= i; j++)
	    {
		double g = 0.0;
		
		for (k = 0; k <= i; k++)
		    g += V[k][i+1] * V[k][j];
		
		for (k = 0; k <= i; k++)
		    V[k][j] -= g * d[k];
	    }
	}
	
	for (k = 0; k <= i; k++)
	    V[k][i+1] = 0.0;
    }
    
    for (j = 0; j < NDIM; j++)
    {
	d[j] = V[NDIM-1][j];
	V[NDIM-1][j] = 0.0;
    }
    
    V[NDIM-1][NDIM-1] = 1.0;
    
    e[0] = 0.0;
}
 

/* Symmetric tridiagonal QL algorithm. */

void tql2(double V[NDIM][NDIM], double d[NDIM], double e[NDIM]) 
{
    int i, j, k, l, m;
    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    
    for (i = 1; i < NDIM; i++)
	e[i-1] = e[i];
    
    e[NDIM-1] = 0.0;
    
    for (l = 0; l < NDIM; l++)
    {
	/* Find small subdiagonal element */

	tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
	
	int m = l;
	
	while (m < NDIM)
	{
	    if (fabs(e[m]) <= eps*tst1)
		break;
	    
	    m++;
	}
	

	/* If m == l, d[l] is an eigenvalue, otherwise, iterate. */

	if (m > l)
	{
	    int iter = 0;
	    
	    do
	    {
		iter = iter + 1;

		/* Compute implicit shift */

		double g = d[l];
		double p = (d[l+1] - g) / (2.0 * e[l]);
		double r = hypot2(p,1.0);
		
		if (p < 0)
		    r = -r;
		
		d[l] = e[l] / (p + r);
		d[l+1] = e[l] * (p + r);
		
		double dl1 = d[l+1];
		double h = g - d[l];
		
		for (i = l+2; i < NDIM; i++)
		    d[i] -= h;
		
		f = f + h;

		/* Implicit QL transformation. */

		p = d[m];
		
		double c = 1.0;
		double c2 = c;
		double c3 = c;
		double el1 = e[l+1];
		double s = 0.0;
		double s2 = 0.0;
		
		for (i = m-1; i >= l; i--)
		{
		    c3 = c2;
		    c2 = c;
		    s2 = s;
		    g = c * e[i];
		    h = c * p;
		    r = hypot2(p,e[i]);
		    e[i+1] = s * r;
		    s = e[i] / r;
		    c = p / r;
		    p = c * d[i] - s * g;
		    d[i+1] = h + s * (c * g + s * d[i]);

		    /* Accumulate transformation. */

		    for (k = 0; k < NDIM; k++)
		    {
			h = V[k][i+1];
			V[k][i+1] = s * V[k][i] + c * h;
			V[k][i] = c * V[k][i] - s * h;
		    }
		}
		
		p = -s * s2 * c3 * el1 * e[l] / dl1;
		e[l] = s * p;
		d[l] = c * p;
		
		/* Check for convergence. */
		
	    } while (fabs(e[l]) > eps*tst1);
	}
	
	d[l] = d[l] + f;
	e[l] = 0.0;
    }
    
    /* Sort eigenvalues and corresponding vectors. */

    for (i = 0; i < NDIM-1; i++)
    {
	double p = d[i];

	k = i;
	
	for (j = i+1; j < NDIM; j++)
	    if (d[j] > p)
	    {
		k = j;
		p = d[j];
	    }
	
	if (k != i)
	{
	    d[k] = d[i];
	    d[i] = p;
	    
	    for (j = 0; j < NDIM; j++)
	    {
		p = V[j][i];
		V[j][i] = V[j][k];
		V[j][k] = p;
	    }
	}
    }
}


void eigen_decomposition(double A[NDIM][NDIM], double d[NDIM], double V[NDIM][NDIM]) 
{
    int i, j;
    double e[NDIM];
    
    for (i = 0; i < NDIM; i++)
	for (j = 0; j < NDIM; j++)
	    V[i][j] = A[i][j];
    
    tred2(V, d, e);
    tql2(V, d, e);
}



int
                lsq_fit(double u[3][3], Transform R)
{
    double           du, omega[6][6], vom[6][6];
    double           dom[6], root2, h[3][3], k[3][3], sign;
    int             i, j, l, rot;

    /* Constant */
    root2 = sqrt(2.0);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = ZERO;

    for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    omega[i][j] = 0;

    /* Calculate determinant of U */
    du = u[0][0] * u[1][1] * u[2][2] - u[0][0] * u[1][2] * u[2][1]
	- u[0][1] * u[1][0] * u[2][2] + u[0][1] * u[1][2] * u[2][0]
	+ u[0][2] * u[1][0] * u[2][1] - u[0][2] * u[1][1] * u[2][0];

    /* If determinant is zero return */
    if (fabs(du) < SMALL)
	return TRUE;

    /* Make metric matrix omega */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    omega[i][j] = ZERO;
	    omega[i + 3][j + 3] = ZERO;
	    omega[i + 3][j] = u[j][i];
	    omega[i][j + 3] = u[i][j];
	}

    /* Diagonalise matrix */
    eigen_decomposition(omega, dom, vom);

    /* Check for degeneracy */
    if (du <= ZERO && fabs(dom[2] - dom[5]) < SMALL)
	return TRUE;

    /* Determine h and k */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    h[i][j] = root2 * vom[i][j];
	    k[i][j] = root2 * vom[i + 3][j];
	}

    sign = h[0][0] * h[1][1] * h[2][2] - h[0][0] * h[1][2] * h[2][1]
	- h[0][1] * h[1][0] * h[2][2] + h[0][1] * h[1][2] * h[2][0]
	+ h[0][2] * h[1][0] * h[2][1] - h[0][2] * h[1][1] * h[2][0];

    if (sign <= ZERO)
	for (i = 0; i < 3; i++)
	{
	    h[i][2] = -h[i][2];
	    k[i][2] = -k[i][2];
	}

    /* Determine rotation matrix */
    du /= fabs(du);
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = k[j][0] * h[i][0] + k[j][1] * h[i][1] + du * k[j][2] * h[i][2];

    R[0][3] = R[1][3] = R[2][3] = R[3][0] = R[3][1] = R[3][2] = ZERO;
    R[3][3] = 1.0;

    return FALSE;
}


/*
 * Trace back highest scoring path 
 */
void
                trace(short *posa, short *posb, int mati, int matj,
		      int pati, int patj, int lasti, int lastj, short *n)
{
    int             pij = pat[pati][patj], i, j;

    for (i = lasti + 1; i < mati; i++)
    {
	*(++posa) = i;
	*(++posb) = 0;
	(*n)++;
    }
    for (j = lastj + 1; j < matj; j++)
    {
	*(++posa) = 0;
	*(++posb) = j;
	(*n)++;
    }
    *(++posa) = mati;
    *(++posb) = matj;
    (*n)++;

    if (!pij)
	return;

    if (pij == 1)
	trace(posa, posb, mati + 1, matj + 1, pati + 1, patj + 1,
	      mati, matj, n);
    if (pij < 1)
	trace(posa, posb, mati + 1, matj - pij, pati + 1, patj - pij,
	      mati, matj, n);
    if (pij > 1)
	trace(posa, posb, mati + pij, matj + 1, pati + pij, patj + 1,
	      mati, matj, n);
}

int
                seqscore(const struct pdbatm *atom1, const struct pdbatm *atom2, const int gappen, const int gapext, ALNS * aln, const int seq1len, const int seq2len)
{
    short          *posa, *posb;
    int             trace_back = (aln != NULL);
    int             now = 0, last = 1;
    int             pati, patj, mati, matj, i, j;
    int             toprows[MAXSEQLEN + 1], topcol, toprow;
    int             maxrows[MAXSEQLEN + 1], maxscore, maxflg = FALSE, maxcol,
                    maxrow, diag, row, col;
    int             mat[2][MAXSEQLEN + 1];

    for (i = 1; i <= seq1len; i++)
    {
	mat[0][i] = mat[1][i] = maxrows[i] = -1000000;
	pat[i][seq2len] = toprows[i] = 0;
    }

    for (j = seq2len; j > 0; j--)
    {
	maxcol = -1000000;
	topcol = 0;

	if (trace_back)
	    pat[seq1len][j] = 0;

	for (i = seq1len; i > 0; i--)
	{
	    /* Get matrix element */
	    mat[now][i] = atom2[j - 1].aac == atom1[i - 1].aac ? 1 : -1;

	    if (j != seq2len && i != seq1len)
	    {
		diag = mat[last][i + 1];
		maxrow = maxrows[i];
		toprow = toprows[i];

		switch (gapext)
		{
		case 1:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) + 1;
			row = maxrow - gappen - (toprow - j) + 1;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		    break;
		case 2:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) - (topcol - i) + 1;
			row = maxrow - gappen - (toprow - j) - (toprow - j) + 1;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		    break;
		case 0:
		    col = maxcol - gappen;
		    row = maxrow - gappen;
		    break;
		default:
		    if (topcol)
		    {
			col = maxcol - gappen - (topcol - i) * gapext + gapext;
			row = maxrow - gappen - (toprow - j) * gapext + gapext;
		    }
		    else
		    {
			col = maxcol - gappen;
			row = maxrow - gappen;
		    }
		}

		if (diag >= col && diag >= row)
		{
		    mat[now][i] += diag;
		    if (trace_back)
			pat[i][j] = 1;
		}
		else
		{
		    if (row > col)
		    {
			mat[now][i] += row;
			if (trace_back)
			    pat[i][j] = -(toprow - j);
		    }
		    else
		    {
			mat[now][i] += col;
			if (trace_back)
			    pat[i][j] = topcol - i;
		    }
		}

		if (diag > maxrows[i])
		{
		    maxrows[i] = diag;
		    toprows[i] = j + 1;
		}
		if (diag > maxcol)
		{
		    maxcol = diag;
		    topcol = i + 1;
		}

		if ((i == 1 || j == 1) && (!maxflg || mat[now][i] > maxscore))
		{
		    maxflg = TRUE;
		    maxscore = mat[now][i];
		    if (trace_back)
		    {
			pati = i;
			patj = j;
			mati = matj = 1;
			if (i == 1)
			    matj = j;
			if (j == 1)
			    mati = i;
		    }
		}
	    }
	}
	now = !now;
	last = !last;
    }

    if (!trace_back)
	return (maxscore);

    posa = aln->posn_a;
    posb = aln->posn_b;
    aln->length = 0;
    trace(posa, posb, mati, matj, pati, patj, 0, 0, &aln->length);
    posa += aln->length;
    posb += aln->length;
    if (*posa == seq1len)
    {
	for (i = *(posb) + 1; i <= seq2len; i++)
	{
	    *(++posa) = 0;
	    *(++posb) = i;
	    (aln->length)++;
	}
	if (aln->length > MAXSEQLEN)
	    puts("score : max. align length exceeded!");
	return (maxscore);
    }
    if (*posb == seq2len)
	for (i = *(posa) + 1; i <= seq1len; i++)
	{
	    *(++posb) = 0;
	    *(++posa) = i;
	    (aln->length)++;
	}
    if (aln->length > MAXSEQLEN)
	puts("score : max. align length exceeded!");

    return (maxscore);
}

/* Get atomic coordinates from ATOM records */
void
                getcoord(x, y, z, chain, n, aacode, atmnam)
    double          *x, *y, *z;
    int            *n, *aacode;
    char           *atmnam, *chain;
{
    int             i, atomn;

    strncpy(atmnam, buf + 12, 4);
    atmnam[4] = '\0';
    if (atmnam[2] == ' ')
	atmnam[3] = ' ';
    sscanf(buf + 6, "%d", &atomn);
    *chain = buf[21];

    sscanf(buf + 22, "%d%lf%lf%lf", n, x, y, z);
    for (i = 0; i < 22; i++)
	if (!strncmp(rnames[i], buf + 17, 3))
	    break;
    *aacode = i;
}

int
                read_atoms(FILE *ifp, int *natoms, struct pdbatm **atmptr)
{
    int             i, token=ENDENT, namino, aac, resnum = 0, blksize, read_end = 0;
    double           x, y, z, u[3][3];
    char            atmnam[5], chain = '?';
    struct pdbatm  *atom = *atmptr;

    blksize = 10 * MEMGRAIN;
    if (!(atom = malloc(blksize * sizeof(struct pdbatm))))
	fail("read_atoms : Out of memory!");

    while (!feof(ifp))
    {
	if (!fgets(buf, 8192, ifp))
	    break;

	if (!isalpha(buf[0]))
	    continue;

	/* Read the record name */
	if (sscanf(buf, "%s", keyword) != 1)
	    continue;

	if (!keyword[0])	/* Odd - there isn't a record name! Skip. */
	    continue;

	token = 0;
	for (i = 1; i <= NTOKENS; i++)	/* Decode record type */
	    if (!strcmp(keyword, tokstr[i - 1]))
		token = i;

	if (token == ENDENT || token == ENDMDL)
	{
	    read_end = 1;
	    break;
	}
	
	switch (token)
	{
	case ATOM:
	    if (buf[16] != ' ' && buf[16] != 'A')
		continue;
	    buf[26] = ' ';
	    if (*natoms >= blksize)
	    {
		blksize += MEMGRAIN;
		atom = realloc(atom, blksize * sizeof(struct pdbatm));
		if (!atom)
		    fail("read_atoms : Out of Memory!");
	    }
	    getcoord(&x, &y, &z, &chain, &namino, &aac, atmnam);
	    if (strncmp(atmnam, " CA", 3))
		continue;
	    resnum++;
	    atom[*natoms].pos[0] = x;
	    atom[*natoms].pos[1] = y;
	    atom[*natoms].pos[2] = z;
	    atom[*natoms].resnum = resnum;
	    atom[*natoms].aac = aac;
	    atom[*natoms].isback = i > 0 && i < 4;
	    strcpy(atom[*natoms].atmnam, atmnam);
	    ++*natoms;
	    break;
	default:		/* Ignore all other types in this version */
	    break;
	}
    }

    *atmptr = atom;

    if (!read_end)
    {
	*natoms = 0;
	return 1;
    }
    

    return 0;
}

main(int argc, char **argv)
{
    int             i, j, k, ii, jj, kk, l, n, nmax = 0, at1, at2, nequiv, **neqmat, nclust[MAXNSTRUC], maxclusize=0, refpdb=FALSE;
    int             blksize, hashval, moda, modb, repnum, naln;
    double           x, y, z, d, r, rmsd, rmsdsum=0.0, u[3][3], **tmscmat, tmsc, matchsum, cutoff, tmscsum=0.0;
    int             mineqv = 56;
    double           eqvdist = 5.93, maxrms = 1000.0, d0sq, rsq, rvar[MAXSEQLEN];
    char            buf[1024], lastid[80], bval[10];
    FILE           *ifp, *ofp;
    Point           new, CG_a, CG_b;
    Transform       fr_xf;
    ALNS tplt;

    if (argc != 2 && argc != 4)
	fail("usage : film3mqap ensemble.pdb [input.pdb output.pdb]");

    ifp = fopen(argv[1], "r");	/* Open PDB file in TEXT mode */
    if (!ifp)
    {
	printf("PDB file error!\n");
	exit(-1);
    }

    /* Read models */
    for (nmodels=0; nmodels<MAXNSTRUC; nmodels++)
	if (read_atoms(ifp, natoms+nmodels, atoms+nmodels))
	    break;

    fclose(ifp);
    
/*    printf("%d models read from PDB file\n", nmodels); */

    /* Calculate all pairwise equivalences */

    for (moda=0; moda<nmodels; moda++)
    {
	nclust[moda] = 0;

	if (natoms[moda] > 0)
	    for (i=0; i<natoms[moda]; i++)
		atoms[moda][i].nequiv = 0;
    }

    for (i=0; i<natoms[moda]; i++)
	rvar[i] = 0.0;

    for (moda=0; moda<nmodels; moda++)
	for (modb=moda+1; modb<nmodels; modb++)
	    if (natoms[moda] > 10 && natoms[modb] > 10)
	    {
		if (seqscore(atoms[moda], atoms[modb], 4, 0, &tplt, natoms[moda], natoms[modb]) < 10)
		    continue;
		
		/* Calculate centroids */
		veczero(CG_a);
		veczero(CG_b);
		for (naln = 0, j = 1; j <= tplt.length; j++)
		    if (tplt.posn_a[j] > 0 && tplt.posn_b[j] > 0)
		    {
			vecadd(CG_a, CG_a, atoms[moda][tplt.posn_a[j]-1].pos);
			vecadd(CG_b, CG_b, atoms[modb][tplt.posn_b[j]-1].pos);
			naln++;
		    }
		
		vecscale(CG_a, CG_a, 1.0 / naln);
		vecscale(CG_b, CG_b, 1.0 / naln);
		
		for (i = 0; i < natoms[moda]; i++)
		    vecsub(atoms[moda][i].pos, atoms[moda][i].pos, CG_a);
		for (j = 0; j < natoms[modb]; j++)
		    vecsub(atoms[modb][j].pos, atoms[modb][j].pos, CG_b);
		
		/* Calculate U */
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++)
			u[i][j] = ZERO;
		
		for (j = 1; j <= tplt.length; j++)
		    if (tplt.posn_a[j] > 0 && tplt.posn_b[j] > 0)
			for (k = 0; k < 3; k++)
			    for (l = 0; l < 3; l++)
				u[k][l] += atoms[moda][tplt.posn_a[j]-1].pos[k] * atoms[modb][tplt.posn_b[j]-1].pos[l];
		
		if (lsq_fit(u, fr_xf))
		    continue; /* LSQ fit failed - all we can do is move to the next pair! */
		
		for (j = 0; j < natoms[modb]; j++)
		{
		    transform_point(fr_xf, atoms[modb][j].pos, new);
		    veccopy(atoms[modb][j].pos, new);
		}
		
		for (nequiv = 0, i = 1; i <= tplt.length; i++)
		    if (tplt.posn_a[i] > 0 && tplt.posn_b[i] > 0)
			nequiv++;
		
		d0sq = SQR(1.24 * pow(nequiv-15.0, 1.0/3.0) - 1.8);
		
		for (rmsd = tmsc = 0, i = 1; i <= tplt.length; i++)
		    if (tplt.posn_a[i] > 0 && tplt.posn_b[i] > 0)
		    {
			rsq = distsq(atoms[moda][tplt.posn_a[i]-1].pos, atoms[modb][tplt.posn_b[i]-1].pos);
			rmsd += rsq;
			tmsc += 1.0F / (1.0F + rsq/d0sq);
			rvar[tplt.posn_a[i]-1] += rsq;
		    }
		
		rmsd = sqrt(rmsd / nequiv);
		tmsc /= nequiv;
		
		if (tmsc > 0.4)
		{
		    nclust[moda]++;
		    nclust[modb]++;
		    
		    if (nclust[moda] > maxclusize)
			maxclusize = nclust[moda];
		    if (nclust[modb] > maxclusize)
			maxclusize = nclust[modb];
		}
		
		tmscsum += tmsc;
		rmsdsum += rmsd;
	    }
    
    tmscsum /= 0.5*nmodels*(nmodels-1);
    rmsdsum /= 0.5*nmodels*(nmodels-1);

    puts("Individual residue RMS values:");
    for (i=0; i<natoms[0]; i++)
	printf("%4d %f\n", i+1, sqrt(2.0*rvar[i]/nmodels/(nmodels-1)));
    
    printf("\nPercentage of structures in largest cluster = %f\nMean TM-Score = %lf\nMean RMSD = %lf\n", 100.0*maxclusize/nmodels, tmscsum, rmsdsum);
    
    printf("\nExpected TM-score for final model = %f\n", 0.661 * tmscsum + 0.289);
    puts("TM-score < 0.4 = Poor model with probably incorrect fold");
    puts("TM-score > 0.5 = Fair model with probably correct fold");
    puts("TM-score > 0.6 = Good model with probably correct fold\n");
    
    if (argc > 2)
    {
	ifp = fopen(argv[2], "r");	/* Open PDB file in TEXT mode */
	if (!ifp)
	    fail("Input PDB file error!\n");

	ofp = fopen(argv[3], "w");	/* Open PDB file in TEXT mode */
	if (!ofp)
	    fail("Output PDB file error!\n");

	lastid[0] = '\0';
	i = -1;
	
	for (;;)
	{
	    if (!fgets(buf, 1024, ifp))
		break;
	    
	    if (!strncmp(buf, "ATOM ", 5))
	    {
		if (strncmp(buf+20, lastid, 7))
		{
		    i++;
		    strncpy(lastid, buf+20, 7);
		}
		sprintf(bval, "%6.2f", sqrt(2.0 * rvar[i]/nmodels/(nmodels-1)));
		strncpy(buf+60, bval, 6);
	    }
	    
	    fprintf(ofp, "%s", buf);
	}
	
	fclose(ofp);
	fclose(ifp);
    }
}
