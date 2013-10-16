/* FILM3 Ensemble Model Recombination Program - by David Jones, November 2011 */
/* Modified by Tim Nugent, July 2013 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#ifdef WIN32
#include <windows.h>
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#endif

#include <iostream> 
#include "pdb.hpp"
#include "ga.hpp"
#include "direct.hpp"
#include <map>
#include <string>
#include <boost/thread.hpp>  

using namespace std;
using namespace boost;

typedef float  Vector[3];

typedef struct
{
    char            *atmnam;
    short           type;
    Vector          pos;
}
SIDECATM;

typedef struct
{
    Vector           n, h, ca, c, o, cb, sc_cg;
    float            phi, psi, omega;
    short            nscats, rotnum, builtflg;
    SIDECATM        *sc;
}
RESATM;

const char *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK", "UNK"
};

const int nscatoms[20] = {
    1, 7, 4, 4, 2, 5, 5, 0, 6, 4, 4, 5, 4, 7, 3, 2, 3, 10, 8, 3
};

/* Maximum chain length */
#define MAXSEQLEN 10000

/* Maximum number of structures per ensemble */
#define MAXNSTRUC 20000

/* Number of annealing steps - can be reduced for smaller proteins */
#define MAXSTEPS 5000000

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#define BIG ((float)1.0e10)
#define SMALL ((float)1.0e-3)
#define ZERO ((float)0.0)
#define HALF ((float)0.5)
#define ONE ((float)1.0)
#define TWO ((float)2.0)
#define THREE ((float)3.0)
#define FOUR ((float)4.0)
#define PI ((float)3.1415927)

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
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define SQR(x) ((x)*(x))
#define dist(a,b) sqrt(SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))

extern void    *
                calloc(), *malloc(), *realloc();

typedef float   Transform[4][4];
typedef float   Point[3];

int topology[100];
int res_type[MAXSEQLEN];
float zcoords[MAXSEQLEN];
float hhlipex[MAXSEQLEN];
int curr_helices = 0;
int helices = 0;
int n_term = 0;
int THREADS = 1;
int MAXGACALLS = 5000;
int GAPOOLSIZE = 500;
int curr_seqlen = 0;
int mem_debug = 0;
int seqlen;
char **seq;
double vbest, z_rot_best = 0, x_rot_best = 0, y_rot_best = 0, z_trans_best = 0, envwt = 0.1, topwt = 0.0, lipexwt = 0.1, x_shift = 0.0, y_shift = 0.0, z_shift = 0.0;

RESATM *chn;
map <const char*, int> aacodes_num; 
map <const char*, int>::iterator it;
int fill_res_type = 0;

typedef struct
{
    short           length;
    short           posn_a[MAXSEQLEN], posn_b[MAXSEQLEN];
}
ALNS;

struct pdbatm
{
    Point           pos;
    float           radius;	/* Radius of sphere */
    int             aac, seqnum, atmtyp, rostyp, dftyp;
    char            atmnam[5];
}              *atoms[MAXNSTRUC];

char            pdbfn[80], csdfn[80], logfn[80], keyword[40], buf[512];

/* Record names for decoding record types */
char           *tokstr[] =
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL",
 "HELIX", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM", "TER",
 "HETATM", "CONECT", "END", "JRNL", "TURN", "ENDMDL"
};

enum RESCODE
{
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys,
    Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val,
    Asx, Glx, Unk
};

struct distrec
{
    short           a, b;
}              *distmat, **sortlist;

int             nmodels, natoms[MAXNSTRUC];
float         **dmat;

char *atnamtab[20] = {
	"N  CA C  O  CB ",
	"N  CA C  O  CB CG CD NE CZ NH1NH2",
	"N  CA C  O  CB CG OD1ND2",
	"N  CA C  O  CB CG OD1OD2",
	"N  CA C  O  CB SG ",
	"N  CA C  O  CB CG CD OE1NE2",
	"N  CA C  O  CB CG CD OE1OE2",
	"N  CA C  O  ",
	"N  CA C  O  CB CG ND1CD2CE1NE2",
	"N  CA C  O  CB CG1CG2CD1",
	"N  CA C  O  CB CG CD1CD2",
	"N  CA C  O  CB CG CD CE NZ ",
	"N  CA C  O  CB CG SD CE ",
	"N  CA C  O  CB CG CD1CD2CE1CE2CZ ",
	"N  CA C  O  CB CG CD ",
	"N  CA C  O  CB OG ",
	"N  CA C  O  CB OG1CG2",
	"N  CA C  O  CB CG CD1CD2NE1CE2CE3CZ2CZ3CH2",
	"N  CA C  O  CB CG CD1CD2CE1CE2CZ OH ",
	"N  CA C  O  CB CG1CG2"
};

int rostypes[20][14] =
{
    { 17,18,19,20,5,0,0,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,4,4,11,1,11,11,0,0,0 },
    { 17,18,19,20,4,1,14,9,0,0,0,0,0,0 },
    { 17,18,19,20,4,2,15,15,0,0,0,0,0,0 },
    { 17,18,19,20,4,16,0,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,4,1,14,9,0,0,0,0,0 },
    { 17,18,19,20,4,4,2,15,15,0,0,0,0,0 },
    { 17,18,19,20,0,0,0,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,1,8,6,6,8,0,0,0,0 },
    { 17,18,19,20,3,4,5,5,0,0,0,0,0,0 },
    { 17,18,19,20,4,3,5,5,0,0,0,0,0,0 },
    { 17,18,19,20,4,4,4,4,10,0,0,0,0,0 },
    { 17,18,19,20,4,4,16,5,0,0,0,0,0,0 },
    { 17,18,19,20,4,1,6,6,6,6,6,0,0,0 },
    { 12,18,19,20,4,4,4,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,13,0,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,13,5,0,0,0,0,0,0,0 },
    { 17,18,19,20,4,1,6,1,7,1,6,6,6,6 },
    { 17,18,19,20,4,1,6,6,6,6,1,13,0,0 },
    { 17,18,19,20,3,5,5,0,0,0,0,0,0,0 }
};    


int isacceptor[20][14] =
{
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,1,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,1,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,1,1,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,1,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,1,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,1,0,0 },
    { 0,0,0,1,0,0,0,0,0,0,0,0,0,0 }
};

int isdonor[20][14] =
{
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,1,0,1,1,0,0,0 },
    { 1,0,0,0,0,0,0,1,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,1,0,0,1,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,1,0,0,0,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,1,0,0,0,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,1,0,0 },
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0 }
};


float rosetta_params[24][2] = {
    { 2.00, 0.1200 },
    { 2.00, 0.1200 },
    { 2.00, 0.0486 },
    { 2.00, 0.1142 },
    { 2.00, 0.1811 },
    { 2.00, 0.1200 },
    { 1.75, 0.2384 },
    { 1.75, 0.2384 },
    { 1.75, 0.2384 },
    { 1.75, 0.2384 },
    { 1.75, 0.2384 },
    { 1.75, 0.2384 },
    { 1.55, 0.1591 },
    { 1.55, 0.1591 },
    { 1.55, 0.2100 },
    { 1.90, 0.1600 },
    { 1.75, 0.2384 },
    { 2.00, 0.0486 },
    { 2.00, 0.1400 },
    { 1.55, 0.1591 },
    { 1.00, 0.0500 },
    { 1.20, 0.0500 },
    { 1.20, 0.0500 },
    { 1.00, 0.0500 }
};


int ambtypes[20][14] =
{
    {
	17, 1, 1,24, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1,19, 3,19,19, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 2,24,14, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 2,25,25, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1,26, 0, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 2,24,14, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 2,25,25, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 5,15, 7, 8,15, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1, 1,20, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1,26, 1, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 3, 3, 3, 3, 3, 3, 0, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1,22, 0, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1,22, 1, 0, 0, 0, 0, 0, 0, 0
    },
    {
	17, 1, 1,24, 1,10, 7, 9,15,11, 3, 3, 3, 3
    },
    {
	17, 1, 1,24, 1, 3, 3, 3, 3, 3, 2,22, 0, 0
    },
    {
	17, 1, 1,24, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0
    }
};


const float mindist[11][4][4] =
{
  {
    { 2.71, 2.28, 3.14, 1.91 },
    { 2.56, 2.73, 3.15, 2.05 },
    { 1.86, 1.73, 1.65, 1.86 },
    { 2.55, 1.39, 3.16, 2.29 },
  },
  {
    { 4.38, 3.62, 3.38, 3.32 },
    { 4.12, 3.68, 2.66, 3.26 },
    { 3.08, 2.61, 2.30, 2.37 },
    { 4.18, 3.55, 2.71, 3.47 },
  },
  {
    { 3.89, 2.90, 3.04, 3.84 },
    { 2.61, 2.32, 2.26, 3.03 },
    { 2.79, 2.06, 2.53, 2.63 },
    { 3.58, 2.69, 2.73, 3.51 },
  },
  {
    { 3.80, 2.77, 2.68, 3.61 },
    { 2.98, 2.49, 1.97, 3.53 },
    { 2.67, 1.73, 2.67, 2.49 },
    { 3.45, 2.33, 2.63, 3.46 },
  },
  {
    { 3.50, 2.51, 3.06, 3.22 },
    { 2.54, 2.50, 1.69, 2.79 },
    { 2.80, 2.03, 2.76, 2.53 },
    { 3.16, 2.94, 2.61, 3.37 },
  },
  {
    { 3.75, 2.78, 3.07, 3.66 },
    { 2.45, 2.22, 2.02, 2.61 },
    { 2.98, 1.83, 3.07, 2.66 },
    { 3.33, 2.90, 2.79, 3.72 },
  },
  {
    { 3.78, 2.45, 3.16, 3.51 },
    { 2.81, 2.83, 2.57, 2.89 },
    { 2.89, 1.83, 3.08, 2.72 },
    { 3.57, 3.25, 2.69, 3.80 },
  },
  {
    { 3.96, 3.12, 3.00, 3.71 },
    { 3.09, 2.51, 2.90, 2.87 },
    { 3.01, 1.86, 3.09, 2.66 },
    { 3.87, 3.41, 2.70, 3.75 },
  },
  {
    { 3.71, 3.30, 3.01, 3.80 },
    { 3.49, 2.64, 1.88, 3.36 },
    { 3.09, 2.62, 3.01, 2.69 },
    { 3.37, 2.52, 2.68, 3.56 },
  },
  {
    { 3.91, 2.63, 3.10, 3.81 },
    { 2.72, 2.62, 2.13, 2.84 },
    { 2.98, 1.93, 3.17, 2.65 },
    { 3.53, 2.34, 2.74, 3.73 },
  },
  {
    { 3.96, 2.86, 3.05, 3.69 },
    { 2.62, 2.80, 2.11, 2.90 },
    { 3.09, 1.85, 3.05, 2.72 },
    { 3.96, 3.24, 2.70, 3.60 },
  }
};

const float maxdist[11][4][4] =
{
  {
    { 4.11, 5.27, 6.18, 2.86 },
    { 5.11, 6.42, 7.40, 3.92 },
    { 3.68, 4.92, 5.81, 2.57 },
    { 5.25, 6.40, 7.24, 3.97 },
  },
  {
    { 7.71, 8.60, 9.55, 6.37 },
    { 8.63, 9.80, 10.78, 7.38 },
    { 7.06, 8.13, 9.03, 5.84 },
    { 8.87, 9.75, 10.64, 7.46 },
  },
  {
    { 11.00, 12.22, 12.86, 9.80 },
    { 12.15, 13.23, 13.86, 10.91 },
    { 10.46, 11.48, 11.98, 9.16 },
    { 12.10, 13.28, 14.00, 10.89 },
  },
  {
    { 14.52, 15.42, 16.19, 13.15 },
    { 15.52, 16.53, 17.22, 14.21 },
    { 13.65, 14.75, 15.29, 12.52 },
    { 15.60, 16.53, 17.41, 14.33 },
  },
  {
    { 17.75, 18.79, 19.68, 16.50 },
    { 18.96, 20.08, 20.66, 17.63 },
    { 16.98, 18.00, 18.86, 15.73 },
    { 19.00, 20.06, 20.78, 17.68 },
  },
  {
    { 21.15, 22.38, 23.02, 19.93 },
    { 22.09, 23.00, 23.59, 20.83 },
    { 20.18, 21.25, 21.92, 18.91 },
    { 22.32, 23.55, 24.07, 21.09 },
  },
  {
    { 24.59, 25.49, 26.34, 23.27 },
    { 25.54, 26.60, 27.03, 24.28 },
    { 23.59, 24.67, 25.39, 22.42 },
    { 25.69, 26.81, 27.46, 24.39 },
  },
  {
    { 27.98, 28.73, 29.73, 26.64 },
    { 28.76, 29.69, 30.34, 27.52 },
    { 26.89, 27.92, 28.55, 25.51 },
    { 28.97, 29.96, 30.88, 27.75 },
  },
  {
    { 31.29, 32.26, 33.10, 29.99 },
    { 32.35, 33.51, 33.44, 31.08 },
    { 30.22, 31.15, 32.11, 29.00 },
    { 32.38, 33.34, 34.30, 31.19 },
  },
  {
    { 34.39, 35.43, 35.84, 33.20 },
    { 34.98, 35.86, 36.51, 34.06 },
    { 33.60, 34.57, 34.84, 32.36 },
    { 35.76, 36.78, 36.96, 34.55 },
  },
  {
    { 37.75, 38.83, 39.41, 36.52 },
    { 38.74, 39.60, 40.25, 37.35 },
    { 36.54, 37.68, 38.52, 35.32 },
    { 38.50, 39.55, 39.50, 37.26 },
  }
};

/* Offsets for residue-specific atom types */
int resatofs[20] = {
    0, 5, 16, 24, 32, 38, 47, 56, 60, 70, 78, 86, 95, 103, 114, 121, 127, 134, 148, 160
};

float dfire[167][167][20];

int pat[MAXSEQLEN][MAXSEQLEN];

float    **conmat, **dcomat;

char *outpdbname;

/***************************************************************************/

void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = (void **) malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	for (i = 0; i < rows; i++)
	    if ((p[i] = calloc(columns, size)) == NULL)
		fail("allocmat: calloc [][] failed!");
    }
    else
	for (i = 0; i < rows; i++)
	    if ((p[i] = malloc(columns * size)) == NULL)
		fail("allocmat: malloc [][] failed!");

    return p;
}

void            freemat(void *p, int rows)
{
    int             i;

    for (i = rows - 1; i >= 0; i--)
	free(((void **) p)[i]);
    free(p);
}

/* PRNG based on Marsaglia's KISS design. This adds a 32-bit multiply-with-carry generator,
   a 3-shift register and a linear congruential generator to give a period of 1.7 x 10^38. */

/* Global KISS state variables: */
unsigned int zmwc=362436069, cmwc=521288629, jsr=123456789, jcong=380116160;

unsigned int kissrng(void)
{
    unsigned long long t;

    jsr ^= (jsr << 13);
    jsr ^= (jsr >> 17);
    jsr ^= (jsr << 5);

    jcong = 1099087573 * jcong + 1234567;
    
    t = 4294957665ULL * zmwc + cmwc;

    cmwc = (t >> 32);
    zmwc = t;
    
    return zmwc + jsr + jcong;
}


/* Calculate uniformly distributed hash of integer */
unsigned int bitmix(unsigned int a)
{
    int i;
    
    for (i=0; i<8; i++)
    {
	a = (a ^ 61) ^ (a >> 16);
	a = a + (a << 3);
	a = a ^ (a >> 4);
	a = a * 0x27d4eb2d;
	a = a ^ (a >> 15);
    }
    
    return a;
}


/*
 * Randomise RNG (initial state must be unique across a large cluster of machines)
 */
void
                randomise(void)
{
#ifdef WIN32
    DWORD i, hash, seed1, seed2, seed3, seed4;
    TCHAR  infoBuf[256];
    DWORD  bufCharCount = 256;
    SYSTEMTIME time_struct;
 
    if (!GetComputerName(infoBuf, &bufCharCount))
	fail("Cannot query computer name!");
    
    for (hash=i=0; i<bufCharCount; i++)
	hash = ((hash << 5) ^ infoBuf[i]) ^ (hash >> 27);
    
    GetLocalTime(&time_struct);
    
    seed1 = getpid();
    seed2 = time_struct.wMilliseconds;
    seed3 = hash;
    seed4 = time(NULL);
#else
    unsigned int i, seed1, seed2, seed3, seed4;
    struct timeval tv;
    FILE *ifp;

    if (gettimeofday(&tv, NULL))
	fail("randomise: cannot generate random number seeds!");
    
    seed1 = getpid();
    seed2 = tv.tv_usec;
    seed3 = gethostid();
    seed4 = tv.tv_sec;
#endif

    zmwc = bitmix(seed1 ^ seed2) ^ bitmix(seed3 ^ seed4);
    cmwc = bitmix(seed1 ^ seed2 ^ seed3 ^ seed4) % 4294957643U + 1;
    jcong = bitmix(seed1 + seed2) + bitmix(seed3 + seed4);
    jsr = bitmix(seed1 + seed2 + seed3 + seed4) % 4294957643U + 1;

    printf("Random seeds: %u %u %u %u\n", zmwc, cmwc, jsr, jcong);
}


/* Generate random number 0<=x<1 */
#define ran0()  (kissrng()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((low) + (int)(((high)-(low)+1) * ran0()))

/* Get atomic coordinates from ATOM records */
void
getcoord(float *x, float *y, float *z, char *chain, int *n, int *aacode, char *atmnam)
{
    int             i;
    char            id[4], temp[9];
    static char     oldchn = '\0';

    atmnam[4] = '\0';
    strncpy(atmnam, buf + 12, 4);
    sscanf(buf + 17, "%3c", id);
    id[3] = '\0';
    *chain = buf[21];

    buf[26] = ' ';
    sscanf(buf + 22, "%d", n);
    temp[8] = '\0';
    strncpy(temp, buf + 30, 8);
    sscanf(temp, "%f", x);
    strncpy(temp, buf + 38, 8);
    sscanf(temp, "%f", y);
    strncpy(temp, buf + 46, 8);
    sscanf(temp, "%f", z);

    for (i = 0; i < 20; i++)
	if (!strncmp(rnames[i], id, 3))
	    break;

    *aacode = i;
}

void read_atoms(FILE *ifp, int *natoms, struct pdbatm **atmptr)
{
    int             i, token=ENDENT, namino, aac, resnum = 0, blksize;
    float           x, y, z, u[3][3];
    char            atmnam[5], chain = '?';
    char *anp;
    struct pdbatm  *atom = *atmptr;

    blksize = 10 * MEMGRAIN;
    if (!(atom = (pdbatm*)malloc(blksize * sizeof(struct pdbatm))))
	fail("read_atoms : Out of memory!");

    while (!feof(ifp))
    {
	if (!fgets(buf, 160, ifp))
	    break;

	/* Read the record name */
	if (sscanf(buf, "%s", keyword) != 1)
	    break;

	if (!keyword[0])	/* Odd - there isn't a record name! Exit. */
	    break;

	token = 0;
	for (i = 1; i <= NTOKENS; i++)	/* Decode record type */
	    if (!strcmp(keyword, tokstr[i - 1]))
		token = i;

	if (token == ENDENT || token == ENDMDL)
	    break;

	switch (token)
	{
	case ATOM:
	    if (buf[16] != ' ' && buf[16] != 'A')
		continue;
	    buf[26] = ' ';
	    if (*natoms >= blksize)
	    {
		blksize += MEMGRAIN;
		atom = (pdbatm*)realloc(atom, blksize * sizeof(struct pdbatm));
		if (!atom)
		    fail("read_atoms : Out of Memory!");
	    }

	    if (!strncmp(buf+13, "CD  ILE", 7))
		strncpy(buf+13, "CD1 ILE", 7);
	    if (!strncmp(buf+13, "O1 ", 3))
		strncpy(buf+13, "O  ", 3);

	    getcoord(&x, &y, &z, &chain, &namino, &aac, atmnam);
	    
	    if (atmnam[0] != ' ' || atmnam[1] == 'H' || atmnam[1] == 'D')
		continue;
/*	    printf(">%s<", atmnam); */
	    if (!strncmp(atmnam, " N  ", 4))
		resnum++;
	    atom[*natoms].pos[0] = x;
	    atom[*natoms].pos[1] = y;
	    atom[*natoms].pos[2] = z;
	    atom[*natoms].seqnum = resnum;
	    atom[*natoms].aac = aac;
	    strcpy(atom[*natoms].atmnam, atmnam);

/*	    puts(atmnam); */

	    anp = strstr(atnamtab[aac], atmnam+1);
	    if (anp)
	    {
		atom[*natoms].atmtyp = (int)(anp - atnamtab[aac])/3;
		atom[*natoms].rostyp = rostypes[aac][(int)(anp - atnamtab[aac])/3];
		atom[*natoms].dftyp = resatofs[aac] + (int)(anp - atnamtab[aac])/3;
/*		printf("%3d %s", atom[*natoms].dftyp, buf); */
	    }
	    else
		continue;
	    ++*natoms;
	    break;
	default:		/* Ignore all other types in this version */
	    break;
	}
    }

    seqlen = MAX(seqlen, resnum);

    *atmptr = atom;
}

/* Write PDB file */
void writepdb(struct pdbatm *atom, int natoms, int steps)
{
    FILE           *ofp;
    int             i,j,atomn;
    float rmsdiff;
    static  int modeln;

    ofp = fopen(outpdbname, "w");
    if (ofp != NULL)
    {

    fprintf(ofp, "HEADER %d\n",steps);    
	for (i = 0; i < natoms; i++)
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
		    i+1, atom[i].atmnam, rnames[atom[i].aac], atom[i].seqnum, atom[i].pos[0], atom[i].pos[1], atom[i].pos[2]);
	fprintf(ofp, "TER\nEND\n");

	fclose(ofp);
    }
    else
	fail("Cannot open output file!");
}


/* Z > 3.5 */
float condmax[20][20] = {
    {  7.4, 8.8, 8.1, 8.8, 8.3, 8.3, 9.8, 7.4, 8.0, 8.5, 9.4, 9.1,10.0,10.0, 7.7, 7.2, 9.0, 8.8,10.0, 9.9 },
    {  8.8,10.0,10.0,10.0,10.0, 9.8,10.0, 8.8, 9.4, 8.0, 9.5,10.0,10.0, 9.9, 8.7, 9.5,10.0, 9.8,10.0, 9.1 },
    {  8.1,10.0, 8.0, 8.1,10.0, 8.5, 9.5, 8.1, 8.1, 9.2, 7.7, 8.7, 8.6,10.0, 7.4, 9.1, 7.5,10.0,10.0, 7.6 },
    {  8.8,10.0, 8.1, 8.8,10.0, 9.1, 9.4, 8.8, 9.1, 8.6, 9.8,10.0, 9.5, 8.8, 7.4, 7.6, 9.1,10.0, 9.8, 9.9 },
    {  8.3,10.0,10.0,10.0, 6.4, 7.9,10.0, 8.3,10.0, 7.9, 8.8,10.0,10.0, 8.7,10.0,10.0,10.0,10.0,10.0, 7.0 },
    {  8.3, 9.8, 8.5, 9.1, 7.9, 7.9, 9.7, 8.3, 9.9, 9.6,10.0, 9.7, 9.1,10.0, 9.7, 9.9, 9.2,10.0,10.0, 9.7 },
    {  9.8,10.0, 9.5, 9.4,10.0, 9.7, 9.5, 9.8, 9.3, 9.7, 9.3,10.0, 9.6, 9.7, 7.3, 8.9, 9.9, 9.9, 9.9, 8.8 },
    {  7.4, 8.8, 8.1, 8.8, 8.3, 8.3, 9.8, 7.4, 8.0, 8.5, 9.4, 9.1,10.0,10.0, 7.7, 7.2, 9.0, 8.8,10.0, 9.9 },
    {  8.0, 9.4, 8.1, 9.1,10.0, 9.9, 9.3, 8.0, 9.1, 8.9, 9.3, 8.8, 9.0,10.0,10.0, 9.8, 9.4,10.0, 9.4, 8.0 },
    {  8.5, 8.0, 9.2, 8.6, 7.9, 9.6, 9.7, 8.5, 8.9, 9.0, 8.7, 8.4, 9.3, 9.6, 8.0, 9.8, 9.3, 9.7,10.0, 8.7 },
    {  9.4, 9.5, 7.7, 9.8, 8.8,10.0, 9.3, 9.4, 9.3, 8.7, 9.2, 9.6, 9.5,10.0, 9.7, 8.6, 9.1,10.0, 9.6, 8.6 },
    {  9.1,10.0, 8.7,10.0,10.0, 9.7,10.0, 9.1, 8.8, 8.4, 9.6,10.0, 9.0, 9.1, 8.3,10.0, 8.9, 7.7,10.0, 9.5 },
    { 10.0,10.0, 8.6, 9.5,10.0, 9.1, 9.6,10.0, 9.0, 9.3, 9.5, 9.0,10.0,10.0,10.0, 8.2, 9.5, 8.6, 8.7, 9.0 },
    { 10.0, 9.9,10.0, 8.8, 8.7,10.0, 9.7,10.0,10.0, 9.6,10.0, 9.1,10.0,10.0, 9.2, 8.7, 9.9,10.0,10.0, 9.4 },
    {  7.7, 8.7, 7.4, 7.4,10.0, 9.7, 7.3, 7.7,10.0, 8.0, 9.7, 8.3,10.0, 9.2,10.0, 7.2, 8.4,10.0,10.0, 9.3 },
    {  7.2, 9.5, 9.1, 7.6,10.0, 9.9, 8.9, 7.2, 9.8, 9.8, 8.6,10.0, 8.2, 8.7, 7.2, 7.6, 8.3, 9.0, 9.1, 8.1 },
    {  9.0,10.0, 7.5, 9.1,10.0, 9.2, 9.9, 9.0, 9.4, 9.3, 9.1, 8.9, 9.5, 9.9, 8.4, 8.3, 8.1, 9.5,10.0, 9.6 },
    {  8.8, 9.8,10.0,10.0,10.0,10.0, 9.9, 8.8,10.0, 9.7,10.0, 7.7, 8.6,10.0,10.0, 9.0, 9.5,10.0,10.0, 9.8 },
    { 10.0,10.0,10.0, 9.8,10.0,10.0, 9.9,10.0, 9.4,10.0, 9.6,10.0, 8.7,10.0,10.0, 9.1,10.0,10.0, 8.5, 9.1 },
    {  9.9, 9.1, 7.6, 9.9, 7.0, 9.7, 8.8, 9.9, 8.0, 8.7, 8.6, 9.5, 9.0, 9.4, 9.3, 8.1, 9.6, 9.8, 9.1, 8.3 }
};

/* Convert AA letter to numeric code (0-22) */
int aanum(int ch)
{
    static int      aacvs[] =
    {
    999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
    22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return isalpha(ch) ? aacvs[ch & 31] : 22;
}


float e_calc(struct pdbatm *atom, int natoms)
{
    int i, j, k, t;
    float esum = 0.0, d, dsq, dviol;

    for (i=0; i<natoms; i++)
	if (atom[i].atmtyp == 1 || atom[i].atmtyp == 4)
	    for (j=i+1; j<natoms; j++)
		if (atom[j].atmtyp == atom[i].atmtyp)
		{
		    t = abs(atom[i].seqnum - atom[j].seqnum);
		    
//		    if (t < 5 && atom[i].atmtyp == 1)
		    if (atom[i].atmtyp == 1)
		    {
			dsq = distsq(atom[i].pos, atom[j].pos);
			if (t == 1 && dsq > 25.0F)
			    esum += 1.0;
			else if (t > 1 && dsq < 20.25F)
			    esum += 0.1;
		    }
		    else if (t >= 5)
		    {
			dsq = distsq(atom[i].pos, atom[j].pos);
			if (dsq < SQR(condmax[atom[i].aac][atom[j].aac]))
			    esum += conmat[atom[i].seqnum-1][atom[j].seqnum-1];
#if 0
			else
			{
			    dviol = condmax[atom[i].aac][atom[j].aac] - sqrtf(dsq);
			    esum += conmat[atom[i].seqnum-1][atom[j].seqnum-1] * expf(-dviol*dviol);
			}
#endif
		    }
		}


    float env_energy = 0;
    if (envwt > 0.0){
       
        PDB* protein = new PDB;
        for (i = 0; i < natoms; i++){  
            //printf("-%s- -%s- -%i- \n",atom[i].atmnam,rnames[atom[i].aac],atom[i].seqnum);
            if(strcmp(atom[i].atmnam," CB ") == 0){
                Point3d p(aacodes_num[rnames[atom[i].aac]],atom[i].seqnum,chn[atom[i].seqnum-1].cb[0],chn[atom[i].seqnum-1].cb[1],chn[atom[i].seqnum-1].cb[2]);
                protein->add_cb(p);            
                //cout << aacodes_num[rnames[atom[i].aac]] << "\t" << rnames[atom[i].aac] << endl;
            }else if(strcmp(atom[i].atmnam," CA ") == 0 && strcmp(rnames[atom[i].aac],"GLY") == 0){
                Point3d p(aacodes_num[rnames[atom[i].aac]],atom[i].seqnum,chn[atom[i].seqnum-1].ca[0],chn[atom[i].seqnum-1].ca[1],chn[atom[i].seqnum-1].ca[2]);
                protein->add_cb(p);  
                //cout << aacodes_num[rnames[atom[i].aac]] << "\t" << rnames[atom[i].aac] << endl;
            }  
        }

        vector<double> shifts = protein->origin_shift();
        x_shift = shifts[0];
        y_shift = shifts[1];
        z_shift = shifts[2];

        // Set search parameters
        vector<double> minparam(3), maxparam(3);

        minparam[0] = 0.0;
        maxparam[0] = 2*M_PI;
        minparam[1] = 0.0;
        maxparam[1] = 2*M_PI;
        minparam[2] = -15.0;
        maxparam[2] = protein->get_maxcdist()+15;

        GA* ga = new GA(minparam, maxparam);
        ga->set_verbose(false);
        ga->set_threads(THREADS);
        ga->set_poolsize(GAPOOLSIZE);
        ga->set_maxcalls(MAXGACALLS);
        ga->set_target(protein);
        vector<double> results = ga->run_ga();
        delete ga;        

        env_energy = protein->orientate(results[0],results[1],results[2]);
        env_energy *= envwt; 
        delete protein;

        //printf("env_energy = %f\n",env_energy);
    }     

    return esum+env_energy;
}

/* Ask Boltzmann if we should accept this energy change! */
int             boltzmann(float de, float t)
{
    float p;

    if (t > 0.0F && de >= 0.0F)
	p = exp(-de / t);
    else
	p = (de > 0.0F) ? 0.0F : 1.0F;

    return ran0() < p;
}

main(int argc, char **argv){

    int             i, j, k, from, start, finish, ii, jj, kk, l, n, nmax = 0, at1, at2, nequiv, **neqmat;
    int             blksize, hashval, nmatch, moda, modb, pass, nuphill, nupmax = 0;
    float           x, y, z, d, dme, r, rmsd, globrmsd, u[3][3], **rmsdmat, matchsum, naln, equivdsq = 64.0, energy[MAXNSTRUC], old_e, new_e, t, t0, ediff, e_min, dst = 0.0, dct = 0.0, maxd, prob;
    char predss[5000], buf[512];
    FILE           *ifp, *ofp;
    struct pdbatm *tmpatm, *curatm;
    Point           CG_a, CG_b;
    Transform       fr_xf;
    ALNS tplt;

    if (argc != 9)
	fail("usage : contactrecomb ensemble.pdb contactfile ssfile zcoords envwt gapoolsize threads out.pdb ");

    char *zfname = argv[4];

    envwt = atof(argv[5]);
    GAPOOLSIZE = atoi(argv[6]);
    THREADS = atoi(argv[7]);
    outpdbname = argv[8];
    
    printf("envwt = %f\n",envwt);
    printf("gapoolsize = %i\n",GAPOOLSIZE);
    printf("threads = %i\n",THREADS);

    randomise();

    ifp = fopen(argv[3], "r");
    if (!ifp)
	fail("main: unable to open SS file!");

    if (!fgets(buf, 160, ifp) || !fgets(buf, 160, ifp) || !fgets(predss, 5000, ifp))
	fail("Bad ESS file!");
	
    fclose(ifp);

    ifp = fopen(argv[1], "r");	/* Open PDB file in TEXT mode */
    if (!ifp)
    {
	printf("Target PDB file error!\n");
	exit(-1);
    }

    /* Read models */
    for (nmodels=0; nmodels<MAXNSTRUC; nmodels++)
    {
	read_atoms(ifp, natoms+nmodels, atoms+nmodels);

	if (!natoms[nmodels])
	    break;
    }

    fclose(ifp);

    fprintf(stderr, "%d models read from PDB files\n", nmodels);

    for (i=1; i<nmodels; i++)
	if (natoms[i] != natoms[0])
	{
	    writepdb(atoms[0], natoms[0],pass);
	    printf("%d %d %d\n", i, natoms[0], natoms[i]);
	    printf("Mismatching number of atoms!\n");
	}

//    nmodels--;

    float tempz;

    ifp = fopen(zfname, "r");
    if (!ifp)
        fail("main: unable to open z-file!");
    int tz = 0;
    for (i=0; i<seqlen; i++)
    {
        if (!fgets(buf, 256, ifp))
        fail("Cannot read from z-file!");
        if (sscanf(buf, "%*s%f", &tempz) == 1){
          zcoords[i] = tempz;
          if(tempz == 15 || tempz == -15){
            topology[tz] = i;
            if(tz == 0 && tempz == 15){
                // 0 == in, 1 == out
                n_term = 1;
            }
            tz++;
          }
        }
    }
    helices = tz/2;

    chn = (RESATM*)calloc(seqlen, sizeof(RESATM));
    curr_seqlen = seqlen;  
    curr_helices = helices;  

    seq = (char **) malloc(sizeof(char *));
    seq[0] = (char*)malloc(seqlen);

    const char* residues[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
    for(int d = 0; d < 20; d++){
        aacodes_num[residues[d]] = d;
    }

    conmat = (float**)allocmat(seqlen, seqlen, sizeof(float), TRUE);    
    dcomat = (float**)allocmat(seqlen, seqlen, sizeof(float), TRUE);    
    
    ifp = fopen(argv[2], "r");
    if (!ifp)
	fail("main: unable to open con file!");
    
    while (!feof(ifp))
    {
	if (!fgets(buf, sizeof(buf), ifp))
	    break;
	
	if (isalpha(buf[0]))
	    continue;
	
	if (sscanf(buf, "%d%d%*s%f%f", &i, &j, &maxd, &prob) != 4)
	    continue;
	
	if (i > seqlen || j > seqlen)
	    continue;
	
	conmat[i-1][j-1] += log(1.0-prob);
//		conmat[i-1][j-1] += prob;
	conmat[j-1][i-1] = conmat[i-1][j-1];
	dcomat[i-1][j-1] = MAX(dcomat[i-1][j-1], maxd);
	dcomat[j-1][i-1] = dcomat[i-1][j-1];
    }
	
    fclose(ifp);

    tmpatm = (pdbatm*)malloc(natoms[0] * sizeof(struct pdbatm));
    curatm = (pdbatm*)malloc(natoms[0] * sizeof(struct pdbatm));

#if 1
    e_min = 1e32;
    
    for (moda=0; moda<nmodels; moda++)
    {
	for (k=0; k<natoms[0]; k++)
	    curatm[k] = atoms[moda][k];
    
	new_e = e_calc(curatm, natoms[0]);

	if (new_e < e_min)
	{
	    writepdb(curatm, natoms[0],pass);
	    old_e = e_min = new_e;
	    printf("Emin = %f\n", e_min);
	    modb = moda;
	}
    }
    for (k=0; k<natoms[0]; k++)
	curatm[k] = atoms[modb][k];
#else
    modb = randint(0, nmodels-1);
    for (k=0; k<natoms[0]; k++)
	curatm[k] = atoms[modb][k];
    e_min = old_e = e_calc(curatm, natoms[0]);
#endif

    printf("Initial E = %f\n", old_e);

    writepdb(curatm, natoms[0],pass);

    t = t0 = 10.0;

    pass = 1;

    nuphill = 0;
    
    while (pass < MAXSTEPS)
    {
	for (k=0; k < natoms[0]; k++)
	    tmpatm[k] = curatm[k];

	n = randint(0, nmodels-1);

	do {
	    start = randint(1, seqlen);
	    finish = randint(start, seqlen);
	} while (predss[start-1] != 'C' || predss[finish-1] != 'C');
	
//	finish = start;
	
/*	finish = start+randint(0, 30);
  finish = MIN(finish, seqlen); */

/*	printf("%d %d %d\n", n, start, finish); */

	for (k=0; k<natoms[0]; k++)
	    if (atoms[n][k].seqnum >= start && atoms[n][k].seqnum <= finish)
		curatm[k] = atoms[n][k];
	
	new_e = e_calc(curatm, natoms[0]);

//	printf("%d %d %f\n", start, finish,  new_e);

	ediff = new_e - old_e;

	if (ediff > 0.0)
	{
	    nuphill++;
	    if (nuphill > nupmax)
		nupmax = nuphill;
	}
	else if (ediff < 0.0)
	    nuphill = 0;

//	t *= 0.9999;
	
	t = 0.001 + 0.1 * log(1.0 + nuphill);

//	t = 0;
	
	if (boltzmann(ediff, t))
	{
	    if (new_e < e_min)
	    {
		e_min = new_e;
		writepdb(curatm, natoms[0],pass);
		printf("Emin = %f\n", e_min);
	    }
	    old_e = new_e;
	}
	else
	    for (k=0; k < natoms[0]; k++)
		curatm[k] = tmpatm[k];

	pass++;
	
	if (pass % 1000 == 0 && pass)
	    printf("%d : T = %f E = %f\n", pass, t, old_e);

	fflush(stdout);
    }
}
