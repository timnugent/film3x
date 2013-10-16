/***************************************************
 *                      FILM3                      *
 *         By David T. Jones      Nov 2011         *
  *        Modified by T. Nugent  Jul 2013           *
 *             UCL Bioinformatics Group            *
 ***************************************************/

/* Parts of this software are Copyright (C) 1995 David T. Jones */
/* Parts of this software are Copyright (C) 2013 Tim Nugent     */

/* This software may not be distributed, copied or used without the permission of the Copyright holders */

#define FILMVersion	"3.03X"
#define Last_Edit_Date	"16th October 2013"

/*
 * This program attempts to fold a sequence into a plausible tertiary
 * structure by joining together protein fragments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#include <stdbool.h>
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

using namespace std;

#define NPAIRS 2

/* Define SAVETRAJ to save intermediate models */
#define noSAVETRAJ

/* Define RASMOL to monitor progress with RASMOL */
#define noRASMOL

/* SUPERSEC defined to use supersecondary fragments */
#define SUPERSEC

#define ELITIST

#define MAXSEQLEN 2000

#define SRDCUTOFF 15.0F
#define LRDCUTOFF 15.0F

#define SRTOPOMAX 8

/* Definitions for potentials */

#define TOPOMAX 11
#define INTERVALS 20
#define NACC 5
#define OOIMIN 6
#define OOIMAX 30
#define OOIDIV 1
#define NRR 11
#define RRWIDTH (1.0F)
#define RTCONST (0.582F)

#define OOICUTOFF (10.0F)


/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128F
#define CACBDIST 1.538F


/* Utility definitions */

#define FALSE 0
#define TRUE 1
#define BIG (1000000)
#define VBIG (1e32F)
#define PI (3.1415927F)


#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define CH alloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define veczero(v) memset(v, 0, sizeof(v))
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) ((a[0]=b[0]*(c)),(a[1]=b[1]*(c)),(a[2]=b[2]*(c)))
#define vecprod(a,b,c) ((a[0]=b[1]*c[2]-b[2]*c[1]),(a[1]=b[2]*c[0]-b[0]*c[2]),(a[2]=b[0]*c[1]-b[1]*c[0]))
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define veccopy(a,b) ((a[0]=b[0]),(a[1]=b[1]),(a[2]=b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define dist(a,b) (sqrtf(distsq(a,b)))

char alnfname[160], confname[160], zfname[160], lipexfname[160];

int rrcbflag;

float     srwt = 1.0, lrwt = 1.0, solvwt = 1.0, hbwt = 0.0, compactwt = 1.0, stericwt = 1.0, dswt = 0.0, rrwt = 0.0, targwt = 0.0;

/* Simulated annealing parameter (total steps) */
int MAXSTEPS = 1000000;

/* Simulated annealing parameter (initial "temperature") */
float INITEMP = 0.75;

/* Simulated annealing parameter (initial "temperature") */
float TRATIO = 0.75;

/* Parameter for thermodynamic annealing */
float KAVALUE = 100000.0;

/* Max. number of supersecondary fragments per residue */
int MAXFRAGS = 5;

/* Max. number of fixed length fragments per residue */
int MAXFRAGS2  = 25;

/* Min. length of fixed length fragments per residue */
int MINFRAGS2LEN  = 9;

/* Max. length of fixed length fragments per residue */
int MAXFRAGS2LEN  = 16;

/* Skip Z-coord test */
int SKIPZTEST  = 0;

/* Skip variable target function and fold one TM helix at a time */
int SEQHELIX = 0;

/* Orientation search type, 1 = GA, 0 = Direct search */
int ORIENTGA = 1;

int last_boundary = 0;

char     buf[MAXSEQLEN], *conseq;
char     tplt_ss[MAXSEQLEN];
short    tcbooi[MAXSEQLEN], gaps[MAXSEQLEN], svrflg[MAXSEQLEN], segidx[MAXSEQLEN], chnbrk[MAXSEQLEN];
int      firstid, nsvr, nseqs;
float    cn_energy, econtrib[MAXSEQLEN], global_e_min = VBIG;

float    param[10], exprad, expsurf, expmean, ecmin, ecmax;

float    ftotwt;

/* Verbosity flag */
int verboseflg = TRUE;

int poolsize = 5000;
int totalsteps = 0;

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

int topology[100];
int res_type[MAXSEQLEN];
float zcoords[MAXSEQLEN];
float hhlipex[MAXSEQLEN];
int curr_helices = 0;
int helices = 0;
int n_term = 0;
int THREADS = 1;
int MAXGACALLS = 1000000;
int GAPOOLSIZE = 500;
int curr_seqlen = 0;
int mem_debug = 0;
int seqlen;
char **seq;
double vbest, z_rot_best = 0, x_rot_best = 0, y_rot_best = 0, z_trans_best = 0, envwt = 0.1, topwt = 0.0, lipexwt = 0.1, x_shift = 0.0, y_shift = 0.0, z_shift = 0.0;

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

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

/* Known disulphide list */
int n_ds = 0;
int ds_from[10], ds_to[10];

char *outpdbn = "fold.pdb";
char *mqapcmd = NULL;

map <const char*, int> aacodes_num; 
map <const char*, int> aacodes_single_num; 
map <const char*, int>::iterator it;

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

enum sstypes
{
  COIL = 1, HELIX = 2, STRAND = 4
};

const char     *sscodes = ".CH.E";

enum PAIRS
{
    CA_CA, CB_CB, CB_N, N_CB, CB_O, O_CB, O_N, N_O
};

enum MODES
{
    FOLDM, REFINEM, MODELM, EVALM
};

const int       pairs[8][2] =
{
    {CAATOM, CAATOM},
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {NATOM, CBATOM},
    {CBATOM, OATOM},
    {OATOM, CBATOM},
    {OATOM, NATOM},
    {NATOM, OATOM}
};

const float mindist[TOPOMAX][4][4] =
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

const float maxdist[TOPOMAX][4][4] =
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

/* DSSP residue solvent accessible area in GGXGG extended pentapeptide */
const float     resacc[22] =
{
    113.0, 253.0, 167.0, 167.0, 140.0, 199.0, 198.0, 88.0, 194.0, 178.0,
    179.0, 215.0, 194.0, 226.0, 151.0, 134.0, 148.0, 268.0, 242.0, 157.0,
    88.0, 88.0
};

/* Amino acid molecular weights */
const float     molwt[22] =
{
    89.09, 174.20, 132.12, 133.10, 121.15, 146.15, 147.13, 75.07, 155.16, 131.17,
    131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.24, 181.19, 117.15,
    89.09, 89.09
};

/* Amino acid radii */
const float     aaradii[22] =
{
    0.77, 2.38, 1.45, 1.43, 1.22, 1.75, 1.77, 0.58, 1.78, 1.56,
    1.54, 2.08, 1.80, 1.90, 1.25, 1.08, 1.24, 2.21, 2.13, 1.29,
    0.77, 0.77
};

const float cgmin[20][20] =
{
    { 3.62, 3.93, 3.77, 3.72, 4.10, 4.03, 3.88, 3.70, 3.98, 4.02, 4.06, 4.18, 4.15, 3.77, 4.06, 3.68, 3.85, 3.84, 3.68, 3.85 },
    { 3.93, 4.40, 4.29, 3.88, 5.07, 4.38, 4.11, 4.15, 4.35, 4.48, 4.35, 4.81, 4.82, 4.31, 4.42, 3.93, 4.15, 4.30, 4.07, 4.32 },
    { 3.77, 4.29, 4.13, 4.03, 4.72, 4.34, 4.29, 3.91, 4.47, 4.49, 4.38, 4.32, 4.70, 4.45, 4.38, 3.92, 4.07, 4.28, 4.29, 4.30 },
    { 3.72, 3.88, 4.03, 4.06, 4.84, 4.36, 4.42, 3.79, 4.36, 4.48, 4.38, 3.95, 4.84, 4.68, 4.28, 3.58, 3.82, 5.32, 4.89, 4.42 },
    { 4.10, 5.07, 4.72, 4.84, 2.68, 4.97, 4.79, 4.13, 5.05, 4.55, 4.58, 5.42, 4.82, 4.40, 4.58, 4.50, 4.46, 5.51, 4.73, 4.34 },
    { 4.03, 4.38, 4.34, 4.36, 4.97, 4.38, 4.53, 4.15, 4.62, 4.65, 4.53, 4.45, 4.83, 4.35, 4.42, 3.98, 4.19, 4.58, 4.30, 4.47 },
    { 3.88, 4.11, 4.29, 4.42, 4.79, 4.53, 4.42, 4.05, 4.56, 4.49, 4.44, 4.05, 4.97, 4.49, 4.30, 3.78, 3.85, 5.00, 4.53, 4.26 },
    { 3.70, 4.15, 3.91, 3.79, 4.13, 4.15, 4.05, 3.82, 4.21, 4.21, 4.23, 4.24, 4.25, 3.73, 4.18, 3.59, 3.82, 3.95, 3.74, 4.06 },
    { 3.98, 4.35, 4.47, 4.36, 5.05, 4.62, 4.56, 4.21, 4.68, 4.55, 4.42, 4.71, 4.97, 4.81, 4.42, 4.17, 4.36, 4.86, 4.10, 4.44 },
    { 4.02, 4.48, 4.49, 4.48, 4.55, 4.65, 4.49, 4.21, 4.55, 4.43, 4.51, 4.81, 4.59, 4.17, 4.61, 4.30, 4.39, 4.47, 4.18, 4.34 },
    { 4.06, 4.35, 4.38, 4.38, 4.58, 4.53, 4.44, 4.23, 4.42, 4.51, 4.54, 4.65, 4.55, 4.28, 4.60, 4.10, 4.35, 4.23, 4.20, 4.36 },
    { 4.18, 4.81, 4.32, 3.95, 5.42, 4.45, 4.05, 4.24, 4.71, 4.81, 4.65, 4.88, 4.99, 4.32, 4.84, 4.18, 4.36, 4.35, 4.19, 4.67 },
    { 4.15, 4.82, 4.70, 4.84, 4.82, 4.83, 4.97, 4.25, 4.97, 4.59, 4.55, 4.99, 4.85, 4.41, 4.84, 4.40, 4.60, 4.98, 4.43, 4.50 },
    { 3.77, 4.31, 4.45, 4.68, 4.40, 4.35, 4.49, 3.73, 4.81, 4.17, 4.28, 4.32, 4.41, 4.46, 4.14, 4.51, 4.36, 5.03, 4.60, 4.07 },
    { 4.06, 4.42, 4.38, 4.28, 4.58, 4.42, 4.30, 4.18, 4.42, 4.61, 4.60, 4.84, 4.84, 4.14, 4.55, 4.11, 4.32, 4.26, 4.00, 4.45 },
    { 3.68, 3.93, 3.92, 3.58, 4.50, 3.98, 3.78, 3.59, 4.17, 4.30, 4.10, 4.18, 4.40, 4.51, 4.11, 3.42, 3.69, 4.61, 4.22, 4.01 },
    { 3.85, 4.15, 4.07, 3.82, 4.46, 4.19, 3.85, 3.82, 4.36, 4.39, 4.35, 4.36, 4.60, 4.36, 4.32, 3.69, 3.87, 5.01, 4.38, 4.26 },
    { 3.84, 4.30, 4.28, 5.32, 5.51, 4.58, 5.00, 3.95, 4.86, 4.47, 4.23, 4.35, 4.98, 5.03, 4.26, 4.61, 5.01, 5.20, 5.12, 4.24 },
    { 3.68, 4.07, 4.29, 4.89, 4.73, 4.30, 4.53, 3.74, 4.10, 4.18, 4.20, 4.19, 4.43, 4.60, 4.00, 4.22, 4.38, 5.12, 4.68, 4.10 },
    { 3.85, 4.32, 4.30, 4.42, 4.34, 4.47, 4.26, 4.06, 4.44, 4.34, 4.36, 4.67, 4.50, 4.07, 4.45, 4.01, 4.26, 4.24, 4.10, 4.21 }
};

/* Dirichlet Mixture Model Substitution Matrix (Crooks-Brenner, 2005) */
const short           dmmmat[23][23] =
{
    {  4, -2, -2, -2,  0, -1, -1, -1, -2, -2, -2, -2, -1, -2, -1,  0, -1, -3, -2, -1, -2, -1, -1 },
    { -2,  6, -1, -1, -4,  1,  0, -3,  0, -3, -3,  2, -2, -3, -2, -1, -1, -2, -2, -3, -1,  0, -1 },
    { -2, -1,  7,  1, -3,  0,  0, -1,  0, -5, -4,  0, -3, -4, -2,  0, -1, -3, -2, -4,  4,  0, -1 },
    { -2, -1,  1,  7, -4,  0,  1, -1, -1, -6, -5,  0, -4, -5, -1,  0, -1, -4, -3, -5,  4,  0, -1 },
    {  0, -4, -3, -4, 12, -3, -4, -3, -3, -1, -2, -4, -1, -2, -3, -2, -1, -2, -2,  0, -3, -3, -1 },
    { -1,  1,  0,  0, -3,  6,  1, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -2, -2, -3,  0,  3,  0 },
    { -1,  0,  0,  1, -4,  1,  5, -2, -1, -4, -4,  1, -3, -4, -1, -1, -1, -3, -3, -4,  0,  3, -1 },
    { -1, -3, -1, -1, -3, -2, -2,  7, -2, -6, -5, -2, -4, -5, -2, -1, -2, -4, -4, -5, -1, -2, -2 },
    { -2,  0,  0, -1, -3,  0, -1, -2,  9, -3, -3, -1, -2, -1, -2, -1, -1,  0,  0, -3,  0,  0,  0 },
    { -2, -3, -5, -6, -1, -3, -4, -6, -3,  5,  2, -4,  1,  0, -4, -4, -2, -1, -1,  3, -5, -3, -1 },
    { -2, -3, -4, -5, -2, -3, -4, -5, -3,  2,  5, -3,  2,  1, -3, -3, -2, -1, -1,  1, -4, -3, -1 },
    { -2,  2,  0,  0, -4,  1,  1, -2, -1, -4, -3,  5, -2, -4, -1, -1, -1, -3, -3, -3,  0,  1, -1 },
    { -1, -2, -3, -4, -1, -2, -3, -4, -2,  1,  2, -2,  7,  1, -3, -2, -1,  0,  0,  1, -3, -2,  0 },
    { -2, -3, -4, -5, -2, -3, -4, -5, -1,  0,  1, -4,  1,  7, -3, -3, -2,  3,  3,  0, -4, -3, -1 },
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -4, -3, -1, -3, -3,  8, -1, -2, -3, -3, -3, -1, -1, -1 },
    {  0, -1,  0,  0, -2,  0, -1, -1, -1, -4, -3, -1, -2, -3, -1,  4,  1, -3, -2, -3,  0,  0, -1 },
    { -1, -1, -1, -1, -1, -1, -1, -2, -1, -2, -2, -1, -1, -2, -2,  1,  5, -2, -2, -1, -1, -1,  0 },
    { -3, -2, -3, -4, -2, -2, -3, -4,  0, -1, -1, -3,  0,  3, -3, -3, -2, 12,  3, -2, -3, -2, -1 },
    { -2, -2, -2, -3, -2, -2, -3, -4,  0, -1, -1, -3,  0,  3, -3, -2, -2,  3,  8, -2, -2, -2, -1 },
    { -1, -3, -4, -5,  0, -3, -4, -5, -3,  3,  1, -3,  1,  0, -3, -3, -1, -2, -2,  5, -4, -3, -1 },
    { -2, -1,  4,  4, -3,  0,  0, -1,  0, -5, -4,  0, -3, -4, -1,  0, -1, -3, -2, -4,  5,  0,  0 },
    { -1,  0,  0,  0, -3,  3,  3, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -2, -2, -3,  0,  4,  0 },
    { -1, -1, -1, -1, -1,  0, -1, -2,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1 }
};

typedef float  Transform[4][4];

/* Side chain rotamer table */

int nrots[20];

SIDECATM rotamers[20][200][10];
float rotafreq[20][200];

struct chnentry
{
    short           length, *relacc;
    char           *seq, *sstruc;
    RESATM         *chain;
} chnlist[20000];

struct chnidx
{
    short           chain, pos, length;
    float           weight, rmsd;
};

struct flist
{
    struct chnidx  *frags;
    int             nfrags, ntot;
    short           nchn, from, bestlen;
    float           totwt, totbest, rmsd, totsum, totsumsq;
    char           *type;
}
fraglist[MAXSEQLEN + 1], fraglist2[MAXSEQLEN + 1];

struct SSTRUC
{
    short           start, length, type;
}
sstlist[100];

int      nchn, chnres;

int targrange = 10000;

float **cb_targ, **cb_last;

FILE    *rfp;

typedef struct
{
    RESATM         *chn;
    float           t, energy;
}
Replica;

float    diffdist[TOPOMAX][4][4];

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVX";

/* Offsets for residue-specific atom types */
int resatofs[20] = {
    0, 5, 16, 24, 32, 38, 47, 56, 60, 70, 78, 86, 95, 103, 114, 121, 127, 134, 148, 160
};

/* Table of mininum distances between 167 atom types */
float atomdistsq[167][167][3];

/* DFIRE potential table */
float dfire[167][167][20];

/* Local conformation tables */
float    sr_de[4][4][21][21][INTERVALS][TOPOMAX];

/* LR Residue-Residue interaction matrices */
float    lr_de[4][4][21][21][NRR];

/* Residue accessibility matrices */
float    acc_de[21][NACC];

/* Residue OOI number matrices */
float    ooi_de[21][(OOIMAX - OOIMIN) / OOIDIV + 1];

/* Structure generation arrays */

float    (**dsqmat)[NPAIRS];

float    ***cbcb_potmat, **solv_potmat, **conmat, **dcomat;

RESATM  *curchn, *bestchn, *oldchn, *targchn;

/* Flags for program modes */
int mod_mode, opt_mode, wt_mode = 1, fragsel;


/* Dump a rude message to standard error and exit */
void fail(char *fmt, ...)
{
    va_list ap;
    
    va_start(ap, fmt) ;
    fprintf(stderr, "*** ");
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    
    exit(-1);
}


/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    static int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return isalpha(ch) ? aacvs[ch & 31] : 22;
}


/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p;

    p = (void**)malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}


/* Free matrix */
void
                freemat(void *p, int rows)
{
    int             i;

    for (i = 0; i < rows; i++)
	free(((void **) p)[i]);
    free(p);
}


/* Byte-swap for little-endian systems */
void byteswap(void *f, register int n)
{
    register unsigned int *vp;

    vp = (unsigned int *)f;

    while (n--)
    {
        *vp = ((*vp >> 24) & 0xFF) | ((*vp >> 8) & 0x0000FF00) | ((*vp << 8) & 0x00FF0000) | ((*vp << 24) & 0xFF000000);
	vp++;
    }
}


/* Implementation of the WELL1024 PRNG by Panneton, L'Ecuyer & Matsumoto (Period 2^1024-1) */

unsigned int wstate[32], widx;

unsigned int WELL1024(void)
{
    unsigned int a, b, c;

    a = wstate[(widx+31) & 31];
    b = wstate[widx] ^ wstate[(widx+3) & 31] ^ (wstate[(widx+3) & 31]>>8);
    c = wstate[(widx+24) & 31] ^ (wstate[(widx+24) & 31]<<19) ^ wstate[(widx+10) & 31] ^ (wstate[(widx+10) & 31]<<14);

    wstate[widx] = b ^ c;
    widx = (widx + 31) & 31;

    return wstate[widx] = (a^(a<<11)) ^ (b^(b<<7)) ^ (c^(c<<13));
}


/* Generate random hash of integer */
unsigned int lcghash(unsigned int key)
{
    int i;
    
    for (i=0; i<4; i++)
    {
	key = 314527869 * key + 1;
	key ^= key >> 22;
    }

    return key;
}


/*
 * Randomise RNG (initial state must be unique across a large cluster of machines)
 */
void
                randomise(void)
{
#ifdef WIN32
    DWORD i, hash=2166136261, seed1, seed2, seed3, seed4;
    TCHAR  infoBuf[256];
    DWORD  bufCharCount = 256;
    SYSTEMTIME time_struct;
 
    if (!GetComputerName(infoBuf, &bufCharCount))
	fail("Cannot query computer name!");
    
    for (i=0; i<bufCharCount; i++)
	hash = (hash * 16777619) ^ infoBuf[i]; /* FNV hash */
    
    GetLocalTime(&time_struct);
    
    seed1 = getpid();
    seed2 = time_struct.wMilliseconds;
    seed3 = hash;
    seed4 = time(NULL);
#else
    unsigned int i, seed1, seed2, seed3, seed4;
    struct timeval tv;

    if (gettimeofday(&tv, NULL))
	fail("randomise: cannot generate random number seeds!");
    
    seed1 = getpid();
    seed2 = tv.tv_usec;
    seed3 = gethostid();
    seed4 = tv.tv_sec;
#endif

    for (i=0; i<32; i++)
	wstate[i] = lcghash(seed1+i) + lcghash(seed2+i) + lcghash(seed3+i) + lcghash(seed4+i);

    widx = 0;

    printf("Random seeds: %u %u %u %u\n", seed1, seed2, seed3, seed4);
}


/* Generate random number 0<=x<1 */
#define ran0()  (WELL1024()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((low) + (int)(((high)-(low)+1) * ran0()))

/* 
   Apply a transformation matrix to a point:
   Transform       transform; transformation to apply to the point
   Vector           p;           the point to transformation
   Vector           tp;  the returned point after transformation
 */
void
                transform_point(Transform transform, Vector p, Vector tp)
{
    Vector           temp;

    temp[0] = p[0] + transform[0][3];
    temp[1] = p[1] + transform[1][3];
    temp[2] = p[2] + transform[2][3];
    tp[0] = dotprod(transform[0], temp);
    tp[1] = dotprod(transform[1], temp);
    tp[2] = dotprod(transform[2], temp);
}


/* Calculate transformation matrix for coordinate frame defined by 3 points */
void
                calcxf(Vector p1, Vector p2, Vector p3, Transform xf)
{
    int             i;
    Vector          rx, ry, rz, temp;
    float           m;

    vecsub(rz, p2, p1);
    m = 1.0F / sqrtf(dotprod(rz, rz));
    vecscale(rz, rz, m);
    vecsub(temp, p3, p1);
    vecprod(rx, temp, rz);
    m = 1.0F / sqrtf(dotprod(rx, rx));
    vecscale(rx, rx, m);
    vecprod(ry, rz, rx);
    for (i = 0; i < 3; i++)
    {
	xf[0][i] = rx[i];
	xf[1][i] = ry[i];
	xf[2][i] = rz[i];
	xf[3][i] = 0.0F;
	xf[i][3] = -p1[i];
    }
    xf[3][3] = 1.0F;
}


/* Position atom d relative to abc with bond length, bond angle phi, and torsion angle theta */
void setdihedral(Vector a, Vector b, Vector c, Vector d, float blencd, float theta, float phi)
{
    float m;
    Vector AB, bc, n, nbc, d2;

    /* Unit vector bc */
    vecsub(bc, c, b);
    m = 1.0F / sqrtf(dotprod(bc, bc));
    vecscale(bc, bc, m);

    d2[0] = blencd * cosf(theta);
    d2[1] = blencd * cosf(phi) * sinf(theta);
    d2[2] = blencd * sinf(phi) * sinf(theta);
    
    /* Vector AB */
    vecsub(AB, b, a);

    /* Unit normal */
    vecprod(n, AB, bc);
    m = 1.0F / sqrtf(dotprod(n, n));
    vecscale(n, n, m);

    vecprod(nbc, n, bc);
 
    /* NeRF method by Parsons et al. */
    d[0] = bc[0] * d2[0] + nbc[0] * d2[1] + n[0] * d2[2];
    d[1] = bc[1] * d2[0] + nbc[1] * d2[1] + n[1] * d2[2];
    d[2] = bc[2] * d2[0] + nbc[2] * d2[1] + n[2] * d2[2];
    
    vecadd(d, d, c);
}


/* Splice a fragment from src to dest+didx */
void splice(RESATM *dest, RESATM *src, int didx, int length)
{
    int             i;
    Transform       xf1, xf2;

/*    printf("splice: %d %d %d\n", didx, length, svrflg[didx]); */
    
    /* Calculate transformation for source fragment based on the first residue (C-alpha at the origin) */
    calcxf(src->n, src->ca, src->c, xf2);

    /* Unless splice is at the beginning we need to transform the residues before insertion point */
    if (didx)
    {
	calcxf(dest[didx].n, dest[didx].ca, dest[didx].c, xf1);
	
	for (i = 0; i < didx; i++)
	{
	    transform_point(xf1, dest[i].n, dest[i].n);
	    transform_point(xf1, dest[i].ca, dest[i].ca);
	    transform_point(xf1, dest[i].c, dest[i].c);
	    transform_point(xf1, dest[i].o, dest[i].o);
	    transform_point(xf1, dest[i].cb, dest[i].cb);
	}
    }
    
    if (didx + length < seqlen)
	   calcxf(dest[didx + length - 1].n, dest[didx + length - 1].ca, dest[didx + length - 1].c, xf1);

    for (i = 0; i < length; i++)
    {
	if (didx + i >= seqlen)
	{
	    printf("*** Cannot carry out splice (%d %i %d %d) !!\n", didx, i, length, seqlen);
        exit(1);
	}

	transform_point(xf2, src[i].n, dest[didx + i].n);
	transform_point(xf2, src[i].ca, dest[didx + i].ca);
	transform_point(xf2, src[i].c, dest[didx + i].c);
	transform_point(xf2, src[i].o, dest[didx + i].o);
	transform_point(xf2, src[i].cb, dest[didx + i].cb);
    }

    /* If the fragment isn't placed at the end then we must transform the following chain segment */ 
    if (didx + length < seqlen){

	   calcxf(dest[didx + length - 1].n, dest[didx + length - 1].ca, dest[didx + length - 1].c, xf2);

	   for (i = didx + length; i < seqlen; i++){
	       transform_point(xf1, dest[i].n, dest[i].n);
	       transform_point(xf1, dest[i].ca, dest[i].ca);
    	   transform_point(xf1, dest[i].c, dest[i].c);
    	   transform_point(xf1, dest[i].o, dest[i].o);
    	   transform_point(xf1, dest[i].cb, dest[i].cb);
	   }

	   for (i = 0; i < didx + length; i++){
    	   transform_point(xf2, dest[i].n, dest[i].n);
    	   transform_point(xf2, dest[i].ca, dest[i].ca);
    	   transform_point(xf2, dest[i].c, dest[i].c);
    	   transform_point(xf2, dest[i].o, dest[i].o);
    	   transform_point(xf2, dest[i].cb, dest[i].cb);
	   }
    }
}


/* Check that copied fragment agrees with sec. struc. */
int chksst(int chain, int from, int to, int len)
{
    int i;

    if (mod_mode == MODELM)
    {
	for (i=0; i<len; i++)
	    if (svrflg[to+i] && !(tplt_ss[to + i] & chnlist[chain].sstruc[from + i]))
		break;
    }
    else
    {
	for (i=0; i<len; i++)
	    if (!(tplt_ss[to + i] & chnlist[chain].sstruc[from + i]))
		break;
    }
    
    return i == len;
}

/* Read in rotamers file */
void readrots(void)
{
    int aa, i, atomn = 0;
    float x, y, z, perc;
    char atmnam[10], buf[160];
    FILE *rfp;

    if (!(rfp = fopen("rotalib.dat", "r")))
	fail("Cannot read rotalib.dat!");

    while (!feof(rfp))
    {
	if (!fgets(buf, 160, rfp))
	    break;

	if (buf[0] == '*')
	{
	    atomn = 0;
	    nrots[aa]++;
	    continue;
	}

	if (buf[1] == ' ')
	{
	    aa = aanum(buf[0]);
	    if (aa >= 20)
		fail("readrots: Unknown aa type!");
	    if (sscanf(buf+2, "%f", &perc) != 1)
		fail("readrots: Bad file format!");
	    rotafreq[aa][nrots[aa]] = perc;
	    continue;
	}
	
	if (sscanf(buf+3, "%f%f%f", &x, &y, &z) != 3)
	    break;

	if (aa != GLY && atomn < nscatoms[aa])
	{
	    rotamers[aa][nrots[aa]][atomn].pos[0] = x;
	    rotamers[aa][nrots[aa]][atomn].pos[1] = y;
	    rotamers[aa][nrots[aa]][atomn].pos[2] = z;
	    sscanf(buf, "%s", atmnam);
	    if ((rotamers[aa][nrots[aa]][atomn].atmnam = strdup(atmnam)) == NULL)
		fail("readrots: Out of memory!");
	    switch(atmnam[0])
	    {
	    case 'C':
		rotamers[aa][nrots[aa]][atomn].type = 0;
		break;
	    case 'N':
		rotamers[aa][nrots[aa]][atomn].type = 1;
		break;
	    case 'O':
		rotamers[aa][nrots[aa]][atomn].type = 2;
		break;
	    case 'S':
		rotamers[aa][nrots[aa]][atomn].type = 3;
		break;
	    default:
		puts(buf);
		fail("readrots: Unknown atom type in rotamer file!");
	    }
	}
	atomn++;
    }

    if (verboseflg)
	for (i=0; i<20; i++)
	    printf("%c %d rotamers\n", rescodes[i], nrots[i]);
}


/* Construct side chains */
void
buildschn(RESATM *chn, int first, int last)
{
    int             aa, i, j;
    Vector           rx, ry, rz, temp, vsum;
    float          m;

    for (i=first; i<=last; i++)
    {
	aa = seq[0][i];

	if (aa >= 20)
	    fail("buildschn: cannot build side chain for unknown residue type!");
	
	vecsub(rz, chn[i].cb, chn[i].ca);
	m = 1.0F / sqrtf(dotprod(rz, rz));
	vecscale(rz, rz, m);
	vecsub(temp, chn[i].n, chn[i].ca);
	vecprod(rx, temp, rz);
	m = 1.0F / sqrtf(dotprod(rx, rx));
	vecscale(rx, rx, m);
	vecprod(ry, rz, rx);

	if (chn[i].sc == NULL && nscatoms[aa])
	{
	    if ((chn[i].sc = (SIDECATM*)malloc(sizeof(SIDECATM) * nscatoms[aa])) == NULL)
		fail("buildschn: cannot allocate memory for side chains!");
	    chn[i].nscats = nscatoms[aa];
	}

	/* Calculate side chain centroid position */
	if (nscatoms[aa])
	{
	    veczero(chn[i].sc_cg);
	    for (j=0; j<nscatoms[aa]; j++)
	    {
		chn[i].sc[j] = rotamers[aa][chn[i].rotnum][j];
		vecscale(vsum, rx, chn[i].sc[j].pos[0]);
		vecscale(temp, ry, chn[i].sc[j].pos[1]);
		vecadd(vsum, vsum, temp);
		vecscale(temp, rz, chn[i].sc[j].pos[2]);
		vecadd(vsum, vsum, temp);
		vecadd(chn[i].sc[j].pos, vsum, chn[i].ca);
		vecadd(chn[i].sc_cg, chn[i].sc_cg, chn[i].sc[j].pos);
	    }
	
	    vecscale(chn[i].sc_cg, chn[i].sc_cg, 1.0F/nscatoms[aa]);
	}
	else
	    veccopy(chn[i].sc_cg, chn[i].ca);
    }
}

/* Main chain hydrogen bond evaluation */
float calchb(float (**dsqmat)[NPAIRS], RESATM * chn){

    int i, j, ndonor[MAXSEQLEN], nacceptor[MAXSEQLEN], nhbs=0, bestj;
    float hb_energy = 0.0F, costheta, e, ebest;
    float dho, dcn, dnca, dno, hvlen;
    Vector c_nvec, ca_nvec, n_hvec;

    /* Build backbone amide hydrogens (Pauling & Corey rule) */
    for (i = 0; i < seqlen; i++){

	   ndonor[i] = nacceptor[i] = 0;

	    /* Generate main chain amide H position bisecting C-N-CA (mean of C->N and CA->N directions) */
    	if (i){
    	    dcn = dist(chn[i-1].c, chn[i].n);
    	    dnca = dist(chn[i].n, chn[i].ca);

    	    vecsub(c_nvec, chn[i].n, chn[i-1].c);
    	    vecscale(c_nvec, c_nvec, 1.0F / dcn);
    	    vecsub(ca_nvec, chn[i].n, chn[i].ca);
    	    vecscale(ca_nvec, ca_nvec, 1.0F / dnca);
    	    
    	    vecadd(n_hvec, c_nvec, ca_nvec);

    	    hvlen = sqrtf(dotprod(n_hvec, n_hvec));
    	    vecscale(n_hvec, n_hvec, 1.0F / hvlen);

    	    vecadd(chn[i].h, chn[i].n, n_hvec);
    	}
    }
    
    for (i = 0; i < seqlen; i++)
	if (seq[0][i] != PRO)
	{
	    /* Find best acceptor for this donor */
	    bestj = -1;
	    ebest = VBIG;
	    for (j = 0; j < seqlen; j++)
		if (abs(j-i) > 2 && dsqmat[i][j][CA_CA] < 81.0F && nacceptor[j] < 2)
		{
		    dno = dist(chn[i].n, chn[j].o);
		    if (dno > 5.0F)
			continue;

		    if (i)
		    {
			/* Dot product N->H . O->H (N.B. N-H vector is already unit length) */
			dho = dist(chn[i].h, chn[j].o);
			
			costheta = ((chn[i].h[0] - chn[i].n[0]) * (chn[i].h[0] - chn[j].o[0]) + (chn[i].h[1] - chn[i].n[1]) * (chn[i].h[1] - chn[j].o[1]) + (chn[i].h[2] - chn[i].n[2]) * (chn[i].h[2] - chn[j].o[2])) / dho;
			
			/* DREIDING hydrogen bond potential with cos^2 weighting */
			e = 5.5F * (5.0F * powf(2.9F / dno, 12.0F) - 6.0F * powf(2.9F / dno, 10.0F)) * SQR(costheta);
		    }
		    else
		    {
			/* Non-directional potential for first N-H */
			e = 5.5F * (5.0F * powf(2.9F / dno, 12.0F) - 6.0F * powf(2.9F / dno, 10.0F));
		    }

		    if (e < ebest)
		    {
			ebest = e;
			bestj = j;
		    }
		}
		    
	    if (ebest < -0.5F)
	    {
		ndonor[i]++;
		nacceptor[bestj]++;
		nhbs++;
		
		if (abs(bestj-i) <= targrange && abs(bestj-i) > 5 && (tplt_ss[i] & STRAND) && (tplt_ss[bestj] & STRAND))
		    hb_energy -= 1.0;
	    }
	}

    //if (verboseflg) printf("NHBs = %d\n", nhbs);
    
    return hb_energy;
}


/* Return total CPU time from first call in seconds */
#ifdef WIN32
float           cputime()
{
    return clock() / (float) CLK_TCK;
}
#else
float           cputime()
{
    float           tot;
    static float    ftime = -1.0, tickspersec = -1.0;
    struct tms      t;

    if (tickspersec < 0.0F)
	tickspersec = sysconf(_SC_CLK_TCK);

    (void) times(&t);

    tot = (float) (t.tms_utime + t.tms_stime) / tickspersec;

    if (ftime < 0.0F)
	ftime = tot;

    return tot - ftime;
}
#endif


double sr_sum, sr_sumsq, lr_sum, lr_sumsq, comp_sum, comp_sumsq;
double steric_sum, steric_sumsq, solv_sum, solv_sumsq, ds_sum, ds_sumsq, rr_sum, rr_sumsq;
double hb_sum, hb_sumsq, targ_sum, targ_sumsq;

float           last_sr, last_lr, last_compact, last_steric, last_hbond, last_solv, last_ds, last_rr, last_rmsd;

/* Compute steric energy contribution */
float calc_steric(float (**dsqmat)[NPAIRS], RESATM * chn, int residx)
{
    int             i, j, k, l, t, t2, atmtyp1, atmtyp2;
    float           steric;
    float           d, dsq, radsum;
    float          *at1, *at2;

    steric = 0.0F;

    for (i = 0; i < seqlen; i++)
    {
	for (j = i+1; j < seqlen; j++)
	{
	    if (residx >= 0 && i != residx && j != residx)
		continue;

	    t = j - i;
	    
	    t2 = MIN(2, t);
	    
	    /* SC-SC */
	    if (dsqmat[i][j][CB_CB] < 14.0F)
		for (k = 0; k < nscatoms[seq[0][i]]; k++)
		{
		    atmtyp1 = resatofs[seq[0][i]]+k+4;
		    
		    for (l = 0; l < nscatoms[seq[0][j]]; l++)
		    {
			atmtyp2 = resatofs[seq[0][j]]+l+4;
			dsq = distsq(chn[i].sc[k].pos, chn[j].sc[l].pos);
			if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			    steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    }
		}
	    
	    if (dsqmat[i][j][CA_CA] < 13.0F)
	    {
		/* SC-MC */
		for (k = 0; k < nscatoms[seq[0][i]]; k++)
		{
		    atmtyp1 = resatofs[seq[0][i]]+k+4;
		    
		    atmtyp2 = resatofs[seq[0][j]];
		    
		    dsq = distsq(chn[i].sc[k].pos, chn[j].n);
		    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp2++;
		    
		    dsq = distsq(chn[i].sc[k].pos, chn[j].ca);
		    if (dsq <= atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp2++;
		    
		    dsq = distsq(chn[i].sc[k].pos, chn[j].c);
		    if (dsq <= atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp2++;
		    
		    dsq = distsq(chn[i].sc[k].pos, chn[j].o);
		    if (dsq <= atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		}
		
		/* MC-SC */
		for (l = 0; l < nscatoms[seq[0][j]]; l++)
		{
		    atmtyp2 = resatofs[seq[0][j]]+l+4;
		    
		    atmtyp1 = resatofs[seq[0][i]];
		    
		    dsq = distsq(chn[i].n, chn[j].sc[l].pos);
		    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp1++;
		    
		    dsq = distsq(chn[i].ca, chn[j].sc[l].pos);
		    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp1++;
		    
		    dsq = distsq(chn[i].c, chn[j].sc[l].pos);
		    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		    
		    atmtyp1++;
		    
		    dsq = distsq(chn[i].o, chn[j].sc[l].pos);
		    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
			steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
		}
		
		/* MC-MC */
		if (dsqmat[i][j][CA_CA] < 9.0F)
		    for (k = 0; k < 4; k++)
		    {
			atmtyp1 = resatofs[seq[0][i]] + k;
			
			switch (k)
			{
			case 0:
			    at1 = chn[i].n;
			    break;
			    
			case 1:
			    at1 = chn[i].ca;
			    break;
			    
			case 2:
			    at1 = chn[i].c;
			    break;
			    
			case 3:
			    at1 = chn[i].o;
			    break;
			}
			
			for (l = 0; l < 4; l++)
			{
			    atmtyp2 = resatofs[seq[0][j]] + l;
			    
			    switch (l)
			    {
			    case 0:
				at2 = chn[j].n;
				break;
				
			    case 1:
				at2 = chn[j].ca;
				break;
				
			    case 2:
				at2 = chn[j].c;
				break;
				
			    case 3:
				at2 = chn[j].o;
				break;
			    }
			    
			    dsq = distsq(at1, at2);
			    if (dsq < atomdistsq[atmtyp1][atmtyp2][t2])
				steric += atomdistsq[atmtyp1][atmtyp2][t2] - dsq;
			}
		    }
	    }
	}
    }

    return steric;
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

float condmin[20][20] = {
    {  3.6, 3.7, 4.8, 3.7, 3.9, 3.6, 3.8, 3.6, 3.6, 4.0, 3.8, 4.0, 3.8, 3.9, 3.8, 3.6, 4.0, 4.0, 3.6, 3.7 },
    {  3.7, 5.6, 4.2, 4.6, 4.2, 4.1, 3.8, 3.7, 4.5, 5.1, 4.5, 4.9, 4.2, 5.2, 4.8, 3.6, 4.1, 4.6, 4.2, 4.4 },
    {  4.8, 4.2, 4.6, 4.7, 4.2, 4.3, 4.0, 4.8, 4.5, 5.1, 4.8, 4.8, 4.7, 4.7, 4.2, 4.4, 4.5, 4.2, 4.5, 5.0 },
    {  3.7, 4.6, 4.7, 5.5, 4.2, 5.2, 4.7, 3.7, 4.2, 4.4, 4.5, 3.9, 4.9, 5.8, 4.1, 3.8, 4.0, 4.2, 3.9, 4.3 },
    {  3.9, 4.2, 4.2, 4.2, 3.7, 4.8, 4.2, 3.9, 4.2, 4.7, 4.3, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.4 },
    {  3.6, 4.1, 4.3, 5.2, 4.8, 4.3, 4.4, 3.6, 3.9, 4.6, 4.4, 4.0, 4.8, 4.6, 4.5, 4.2, 3.9, 4.2, 4.9, 4.2 },
    {  3.8, 3.8, 4.0, 4.7, 4.2, 4.4, 5.9, 3.8, 4.3, 4.5, 4.2, 4.1, 5.0, 4.3, 4.5, 3.8, 4.4, 5.3, 4.3, 3.9 },
    {  3.6, 3.7, 4.8, 3.7, 3.9, 3.6, 3.8, 3.6, 3.6, 4.0, 3.8, 4.0, 3.8, 3.9, 3.8, 3.6, 4.0, 4.0, 3.6, 3.7 },
    {  3.6, 4.5, 4.5, 4.2, 4.2, 3.9, 4.3, 3.6, 4.5, 4.6, 4.1, 4.3, 5.0, 4.1, 4.3, 4.4, 4.3, 4.2, 5.2, 4.7 },
    {  4.0, 5.1, 5.1, 4.4, 4.7, 4.6, 4.5, 4.0, 4.6, 4.3, 4.1, 4.8, 4.7, 3.9, 4.4, 4.3, 4.1, 4.2, 4.4, 4.1 },
    {  3.8, 4.5, 4.8, 4.5, 4.3, 4.4, 4.2, 3.8, 4.1, 4.1, 4.1, 4.2, 4.2, 4.1, 4.1, 3.8, 4.3, 4.5, 4.1, 3.9 },
    {  4.0, 4.9, 4.8, 3.9, 4.2, 4.0, 4.1, 4.0, 4.3, 4.8, 4.2, 4.6, 4.9, 4.2, 4.6, 3.8, 4.4, 5.0, 4.1, 4.3 },
    {  3.8, 4.2, 4.7, 4.9, 4.2, 4.8, 5.0, 3.8, 5.0, 4.7, 4.2, 4.9, 4.2, 5.4, 4.2, 4.3, 4.2, 5.0, 4.4, 4.1 },
    {  3.9, 5.2, 4.7, 5.8, 4.2, 4.6, 4.3, 3.9, 4.1, 3.9, 4.1, 4.2, 5.4, 4.1, 4.3, 3.9, 4.3, 5.2, 4.7, 4.1 },
    {  3.8, 4.8, 4.2, 4.1, 4.2, 4.5, 4.5, 3.8, 4.3, 4.4, 4.1, 4.6, 4.2, 4.3, 4.2, 4.5, 4.8, 4.2, 4.6, 4.1 },
    {  3.6, 3.6, 4.4, 3.8, 4.2, 4.2, 3.8, 3.6, 4.4, 4.3, 3.8, 3.8, 4.3, 3.9, 4.5, 4.1, 4.0, 3.9, 4.2, 4.0 },
    {  4.0, 4.1, 4.5, 4.0, 4.2, 3.9, 4.4, 4.0, 4.3, 4.1, 4.3, 4.4, 4.2, 4.3, 4.8, 4.0, 4.0, 4.5, 4.9, 4.2 },
    {  4.0, 4.6, 4.2, 4.2, 4.2, 4.2, 5.3, 4.0, 4.2, 4.2, 4.5, 5.0, 5.0, 5.2, 4.2, 3.9, 4.5, 4.2, 4.8, 4.8 },
    {  3.6, 4.2, 4.5, 3.9, 4.2, 4.9, 4.3, 3.6, 5.2, 4.4, 4.1, 4.1, 4.4, 4.7, 4.6, 4.2, 4.9, 4.8, 4.7, 4.4 },
    {  3.7, 4.4, 5.0, 4.3, 4.4, 4.2, 3.9, 3.7, 4.7, 4.1, 3.9, 4.3, 4.1, 4.1, 4.1, 4.0, 4.2, 4.8, 4.4, 4.3 }
};

float e_model(float (**dsqmat)[NPAIRS], RESATM * chn, int debug)
{
    int             i, j, k, t, n_sr, n_lr, n_rr;
    float           p, sr, lr, compact, steric, hbond = 0.0, env_energy = 0.0, top_energy = 0.0, solv, ds, rrcon, lclosure;
    float           dmax, dmean, rmsd, dviol;
    short           dspart[MAXSEQLEN];

    rrcon = 0;

#if 0
    for (n_sr = i = 0; i < seqlen-4; i++)
	if (tplt_ss[i] == HELIX && tplt_ss[i+1] == HELIX && tplt_ss[i+2] == HELIX && tplt_ss[i+3] == HELIX && tplt_ss[i+4] == HELIX)
	{
	    dviol = fabsf(sqrtf(dsqmat[i][i+4][CA_CA]) - 6.2F);
	    if (dviol > 0.5F)
		rrcon += (dviol - 0.5F);
	    n_sr++;
	}

//    rrcon /= n_sr;
#endif

    for (i = 0; i < MIN(curr_seqlen,seqlen); i++)
	for (j = MIN(i+targrange, seqlen-1); j >= i+5; j--)
	{	    
	    if (conmat[i][j] != 0.0F)
	    {
		if (j-i > 5 && (tplt_ss[i] & STRAND) && (tplt_ss[j] & STRAND))
		{
		    dviol = sqrtf(dsqmat[i][j][CA_CA]) - 6.5F;
		}
		else
		    dviol = dist(chn[i].cb, chn[j].cb) - condmax[conseq[i]][conseq[j]];
		
		if (dviol <= 0.0F)
		    rrcon += conmat[i][j];
		else
		    rrcon += conmat[i][j] * expf(-dviol*dviol);
	    }
	}

    int z_test_score = 0;

    if (hbwt > 0.0)
        hbond = calchb(dsqmat, chn);
 
    if (envwt > 0.0){

        PDB* protein = new PDB;
        for (i = 0; i < MIN(curr_seqlen,seqlen); i++){
            if(res_type[i] != 7){
                Point3d p(res_type[i],i,chn[i].cb[0],chn[i].cb[1],chn[i].cb[2]);
                protein->add_cb(p);
            }else{
                Point3d p(res_type[i],i,chn[i].ca[0],chn[i].ca[1],chn[i].ca[2]);
                protein->add_cb(p);                
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

        //cout << "Orientating.... energy = ";

        if(ORIENTGA){

            GA* ga = new GA(minparam, maxparam);
            ga->set_verbose(false);
            ga->set_threads(THREADS);
            ga->set_poolsize(GAPOOLSIZE);
            ga->set_maxcalls(MAXGACALLS);
            ga->set_target(protein);
            vector<double> results = ga->run_ga();
            delete ga;        

            x_rot_best = results[0];
            y_rot_best = results[1];
            z_trans_best = results[2];

        }else{

            double startpt[3], results[3];
            double prevbest = 1000000;
            // Quick grid search to get starting values
            for(double x = 0; x < maxparam[0]; x+= 0.104719755){
                for(double y = 0; y < maxparam[0]; y+= 0.1047197555){
                    for(double z = minparam[2]; z < maxparam[2]; z+= 5){
                        double lowest_e = protein->orientate(x,y,z);
                        if(lowest_e < prevbest){
                            prevbest = lowest_e;
                            startpt[0] = x;
                            startpt[1] = y;
                            startpt[2] = z;
                        }   
                    }
                }
            }
            Pattern* pa = new Pattern();
            pa->set_target(protein);
            pa->set_verbose(false);
            int pa_calls = pa->hooke(3, startpt, results);
            delete pa;

            x_rot_best = results[0];
            y_rot_best = results[1];
            z_trans_best = results[2];

        }

        env_energy = protein->orientate(x_rot_best,y_rot_best,z_trans_best);
        env_energy *= envwt; 
        delete protein;
        //printf("%3.3f\n",env_energy);
        
    }    

    if (debug > 0 && verboseflg){
    	printf("Target range = %d\n", targrange);
        printf("Current seqlen = %d\n", curr_seqlen);
    	printf("Econtact = %f\n", rrcon);
        printf("Eenv = %f\n", env_energy);
    	printf("Ehb = %f\n", hbond);
    	printf("TOT = %f\n\n", rrcon + hbond + env_energy + top_energy);
    	fflush(stdout);
    }

    return rrcon + hbond + env_energy + top_energy;
}


/* Generate distance matrix and Ooi numbers - return TRUE immediately if chain overlaps itself */
int calcdist(const RESATM * chain, const int check_clash)
{
    int             i, j, t, tt;
    static const int tpc[6] = {4, 5, 6, 3, 7, 2}; /* List of most likely clashing seq. separations in order */

    /* If checking for clashes scan residues close in sequence first */
    if (check_clash)
    {
	for (tt = 2; tt < seqlen; tt++)
	{
	    if (tt >= 8)
		t = tt;
	    else
		t = tpc[tt-2];
		
	    for (i = 0; i < seqlen-t; i++)
	    {
		j = i+t;
		dsqmat[i][j][CA_CA] = distsq(chain[i].ca, chain[j].ca);
		if (dsqmat[i][j][CA_CA] < 20.25F)
		    return TRUE;
	    }
	}
	
	for (i = 0; i < seqlen-1; i++)
	    dsqmat[i][i+1][CA_CA] = distsq(chain[i].ca, chain[i+1].ca);
    }
    else
	for (i = 0; i < seqlen; i++)
	    for (j = i + 1; j < seqlen; j++)
		dsqmat[i][j][CA_CA] = distsq(chain[i].ca, chain[j].ca);

    for (i = 0; i < seqlen; i++)
	for (j = i + 1; j < seqlen; j++)
	    dsqmat[j][i][CA_CA] = dsqmat[i][j][CA_CA];

    return FALSE;
}


/* Count severe CA-CA clashes */
int clashcount(const RESATM * chain, int length)
{
    int             i, j, nclash=0;
    float dsq;

    for (i = 0; i < length; i++)
	for (j = i + 2; j < length; j++)
	{
	    dsq = distsq(chain[i].ca, chain[j].ca);
 	    if (dsq < 12.25F)
		nclash++;
	}

    return nclash;
}

int ztest(float (**dsqmat)[NPAIRS]){

    if(SKIPZTEST) return FALSE;

    int i, j;    
    for (i=0; i<seqlen; i++)
	if (zcoords[i] < 1000.0F)
	    for (j=i+5; j<seqlen; j++)
		if (zcoords[j] < 1000.F && dsqmat[i][j][CA_CA] < SQR(zcoords[i] - zcoords[j]) - 36.0F)
		    return TRUE;
    
    return FALSE;
}

/* Copy chain segment with checking for side chain atom memory allocation */
void chaincpy(RESATM *dest, RESATM *src, int n)
{
    SIDECATM        *oldsc;

    while (n--)
    {
	if (dest->sc != NULL && src->sc == NULL)
	    continue;
	oldsc = dest->sc;
	memcpy(dest, src, sizeof(RESATM));
	dest->sc = oldsc;

	if (src->sc != NULL)
	{
	    if (dest->sc == NULL)
		if ((dest->sc = (SIDECATM*)malloc(sizeof(SIDECATM) * dest->nscats)) == NULL)
		    fail("chaincpy: cannot allocate memory for side chains!");
	
	    memcpy(dest->sc, src->sc, dest->nscats * sizeof(SIDECATM));
	}
	
	dest++;
	src++;
    }
}


void
                readchn(void)
{
    FILE           *dfp, *tfp;
    char            brkid[10], brkidlst[100000], tdbname[512], buf[512], seq[MAXSEQLEN], sstruc[MAXSEQLEN];
    int             i, n, helix = 0, strand = 0, tot = 0, relacc[MAXSEQLEN];
    RESATM          chain[MAXSEQLEN];
    //extern char    *getenv();

    if (!(tfp = fopen("tor.lst", "r")))
	fail("readchn: cannot open tor.lst!");

    brkidlst[0] = '\0';
    while (!feof(tfp))
    {
	if (fscanf(tfp, "%s", brkid) != 1)
	    break;
	if (brkid[0] == '#')
	    continue;
	strcat(brkidlst, brkid);
    }

    fclose(tfp);

    /* Read coords from PDB file */
    if (getenv("TDBDAT"))
	strcpy(tdbname, getenv("TDBDAT"));
    else
	strcpy(tdbname, "tdb.dat");
    dfp = fopen(tdbname, "r");
    if (!dfp)
	fail("readchn: cannot open tdb.dat file!");

    if (!fgets(buf, 160, dfp)){
	   printf("tdb error1:\n%s\n",buf);
       //fail("readchn: bad tdb.dat file!");
    }

    while (!feof(dfp))
    {
	if (buf[0] != '#'){
        printf("tdb error2:\n%s\n",buf);
	    //fail("readchn: bad tdb.dat file!");
    }

	if (sscanf(buf, "%*s%s", brkid) != 1)
	    break;

	n = 0;
	
	while (!feof(dfp))
	{
	    if (!fgets(buf, 512, dfp))
		break;

	    if (buf[0] == '#')
		break;

	    if (!strstr(brkidlst, brkid))
		continue;
	    
	    if (n >= MAXSEQLEN)
		fail("readchn: MAXSEQLEN exceeded!");
	    
	    if (buf[7] == 'H' || buf[7] == 'G')
	    {
		helix++;
		sstruc[n] = HELIX;
	    }
	    else if (toupper(buf[7]) == 'E' || toupper(buf[7]) == 'A' || toupper(buf[7]) == 'P')
	    {
		strand++;
		sstruc[n] = STRAND;
	    }
	    else
		sstruc[n] = COIL;
	    tot++;
	    sscanf(buf + 9, "%d", &relacc[n]);
	    seq[n] = aanum(buf[5]);
	    if (sscanf(buf + 13, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &chain[n].phi, &chain[n].psi, &chain[n].omega, chain[n].n, chain[n].n + 1, chain[n].n + 2, chain[n].ca, chain[n].ca + 1, chain[n].ca + 2, chain[n].c, chain[n].c + 1, chain[n].c + 2, chain[n].o, chain[n].o + 1, chain[n].o + 2, chain[n].cb, chain[n].cb + 1, chain[n].cb + 2) != 18)
		  printf("readchn: bad tdb file!!\n");
        //fail("readchn: bad tdb file!!");
	    chain[n].sc = NULL;
	    chain[n].rotnum = 0;
	    
	    n++;
	}

#if 0
/* Check for mistakes in chain geometry - use only when creating new tdb.dat file */
	if (n < 20)
	    continue;

	for (i = 0; i < n; i++)
	    if (seq[i] >= 20)
		break;
	
	if (i != n)
	{
	    printf("Unknown residue in chain!\n");
	    continue;
	}
	
	for (i = 0; i < n - 1; i++)
	{
	    bondl = dist(chain[i].c, chain[i + 1].n);
	    if (bondl > 1.45F)
		break;

#if 0
	    printf("! ");
	    for (j=1; j<8; j++)
	    {
		if (i+j < n)
		{
		    printf(" %7.3f", dist(chain[i].ca, chain[i + j].ca));
		}
		else
		    printf(" %7.3f", 0.0);
	    }
	    putchar('\n');
#endif
	}

	if (i != n - 1)
	{
	    printf("Bondlen = %f!\n", bondl);
	    continue;
	}

	/* Check backbone geometry */
	for (i = 0; i < n; i++)
	{
	    bondl = dist(chain[i].n, chain[i].c);
	    if (fabsf(bondl - 2.46F) / 0.0475F > 3.0F)
	    {
		printf("Bad N..C dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].n, chain[i].ca);
	    if (fabsf(bondl - 1.46F) / 0.0162F > 3.0F)
	    {
		printf("Bad N-CA dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].ca, chain[i].c);
	    if (fabsf(bondl - 1.52F) / 0.0159F > 3.0F)
	    {
		printf("Bad CA-C dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	    bondl = dist(chain[i].c, chain[i].o);
	    if (fabsf(bondl - 1.23F) / 0.0103F > 4.0F)
	    {
		printf("Bad C=O dist at %d = %f!\n", i+1, bondl);
		break;
	    }
	}

	if (i != n)
	    continue;
	
	for (i = 0; i < n - 1; i++)
	{
	    bondl = dist(chain[i].ca, chain[i + 1].ca);

	    /* Check virtual bond distance for cis and trans peptides */
	    if (bondl < 2.8F || bondl > 3.9F || (bondl >= 3.2F && bondl < 3.7F))
	    {
		printf("Bad virtual bond (%f A)!\n", bondl);
		break;
	    }
	}

	if (i != n-1)
	    continue;

	if (clashcount(chain, n) > 0)
	{
	    printf("Bad CA-CA bumps in %s chain!\n", brkid);
	    continue;
	}

	printf("OKCHAIN: %s\n", brkid);
#endif

	chnlist[nchn].length = n;
	if (!(chnlist[nchn].chain = (RESATM*)malloc(n * sizeof(RESATM))) ||
	    !(chnlist[nchn].relacc = (short int*)malloc(n * sizeof(short))) ||
	    !(chnlist[nchn].seq = (char*)malloc(n)) ||
	    !(chnlist[nchn].sstruc = (char*)malloc(n)))
	        fail("readchn: malloc failed!");

	for (i = 0; i < n; i++)
	{
	    chaincpy(chnlist[nchn].chain + i, chain + i, 1);
	    chnlist[nchn].relacc[i] = relacc[i];
	    chnlist[nchn].seq[i] = seq[i];
	    chnlist[nchn].sstruc[i] = sstruc[i];
	}

	nchn++;
	chnres += n;
    }
    
    fclose(dfp);

    if (verboseflg)
    {
	printf("%d chains read (%d residues)\n", nchn, chnres);
	printf("%5.1f %% helix\n", 100.0 * helix / tot);
	printf("%5.1f %% strand\n", 100.0 * strand / tot);
    }
}


/* Generate fully extended chain */
void
                fullyextended(RESATM *chn)
{
    int             i;
    float           sx, sy;
    Vector          cca, nca, xx, yy;

    /* Generate fully extended chain to start without clashes */
    chn[0].n[0] = chn[0].n[1] = chn[0].n[2] = 0.0;
    chn[0].ca[0] = 1.45;
    chn[0].ca[1] = chn[0].ca[2] = 0.0;
    chn[0].c[0] = chn[0].ca[0] + cos(1.2) * 1.52;
    chn[0].c[1] = sin(1.2) * 1.52;
    chn[0].c[2] = 0.0;
    
    setdihedral(chn[0].n, chn[0].ca, chn[0].c, chn[0].o, 1.23, 58.9 * PI / 180.0, 0.0);
    
    for (i=1; i<seqlen; i++)
    {
	setdihedral(chn[i-1].n, chn[i-1].ca, chn[i-1].c, chn[i].n, 1.33, 64.4 * PI / 180.0, PI);
	setdihedral(chn[i-1].ca, chn[i-1].c, chn[i].n, chn[i].ca, 1.45, 58.1 * PI / 180.0, PI);
	setdihedral(chn[i-1].ca, chn[i-1].c, chn[i].n, chn[i].h, 1.0, 58.1 * PI / 180.0, 0.0);
	setdihedral(chn[i-1].c, chn[i].n, chn[i].ca, chn[i].c, 1.52, 68.8 * PI / 180.0, PI);
	setdihedral(chn[i].n, chn[i].ca, chn[i].c, chn[i].o, 1.23, 58.9 * PI / 180.0, 0.0);
    }
    
    /* Generate CB positions */
    for (i=0; i<seqlen; i++)
    {
	vecsub(nca, chn[i].ca, chn[i].n);
	vecsub(cca, chn[i].ca, chn[i].c);
	vecadd(xx, nca, cca);
	vecprod(yy, nca, cca);
	sx = CACBDIST * cosf(TETH_ANG) / sqrtf(dotprod(xx, xx));
	sy = CACBDIST * sinf(TETH_ANG) / sqrtf(dotprod(yy, yy));
	vecscale(xx, xx, sx);
	vecscale(yy, yy, sy);
	vecadd(chn[i].cb, chn[i].ca, xx);
	vecadd(chn[i].cb, chn[i].cb, yy);
    }
}

/* Estimate coordinate error from energy contribution (NOT IMPLEMENTED) */
#define errorest(x) (9.999)

void transform_atom(float *x, float *y, float *z){

    float x1 = *x - x_shift;
    float y1 = *y - y_shift;
    float z1 = *z - z_shift;
    float xp, yp, zp;


    // X-axis rotation
    yp = (y1*cos(x_rot_best)) - (z1*sin(x_rot_best));
    zp = (y1*sin(x_rot_best)) + (z1*cos(x_rot_best));
    y1 = yp;
    z1 = zp;

    // Y-axis rotation
    zp = (z1*cos(y_rot_best)) - (x1*sin(y_rot_best));
    xp = (z1*sin(y_rot_best)) + (x1*cos(y_rot_best));
    z1 = zp;
    x1 = xp;

    // Z-axis translation
    z1 += z_trans_best;

    *x = x1;
    *y = y1;
    *z = z1;
    return;
}

void fit_chain(RESATM * chain){

    int i,j,atomn;
    for (atomn = i = 0; i < seqlen; i++){
        //printf("XYZ before transformation: %f %f %f\n",chain[i].ca[0],chain[i].ca[1],chain[i].ca[2]);
        transform_atom(&chain[i].n[0], &chain[i].n[1], &chain[i].n[2]);
        transform_atom(&chain[i].ca[0], &chain[i].ca[1], &chain[i].ca[2]);
        //printf("XYZ after transformation: %f %f %f\n",chain[i].ca[0],chain[i].ca[1],chain[i].ca[2]);
        transform_atom(&chain[i].c[0], &chain[i].c[1], &chain[i].c[2]);
        transform_atom(&chain[i].o[0], &chain[i].o[1], &chain[i].o[2]);
        transform_atom(&chain[i].cb[0], &chain[i].cb[1], &chain[i].cb[2]);

        //cout << i << "\t" << chain[i].cb[2] << endl;

        if (chain[i].sc && seq[0][i] != GLY && seq[0][i] != ALA){
            for (j = 1; j < nscatoms[seq[0][i]]; j++){
                transform_atom(&chain[i].sc[j].pos[0], &chain[i].sc[j].pos[1], &chain[i].sc[j].pos[2]);  
            }      
        }
    }
    return;    
}    

/* Write PDB file */
void writepdb(RESATM * chain, float e){

    FILE           *ofp;
    int             i,j,atomn;
    float rmsdiff;
    char outpdbn2[160];
    static  int modeln;
    float max_x = -BIG, max_y = -BIG, min_x = BIG, min_y = BIG;

    if (envwt > 0.0){
        fit_chain(chain);
    }

    ofp = fopen(outpdbn, "w");
    if (ofp != NULL)
    {
	if (mqapcmd == NULL)
	    fprintf(ofp, "HEADER Potential terms: %f %f %i\n", e, vbest, totalsteps);
	else
	    fprintf(ofp, "HEADER %f\n", e);

	for (atomn = i = 0; i < MIN(curr_seqlen,seqlen); i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " N  ", rnames[seq[0][i]], i + 1, chain[i].n[0], chain[i].n[1], chain[i].n[2], errorest(econtrib[i]));
#if 0
	    if (hbwt > 0.0 && i > 0 && seq[0][i] != PRO)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " H  ", rnames[seq[0][i]], i + 1, chain[i].h[0], chain[i].h[1], chain[i].h[2], errorest(econtrib[i]));
#endif
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " CA ", rnames[seq[0][i]], i + 1, chain[i].ca[0], chain[i].ca[1], chain[i].ca[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " C  ", rnames[seq[0][i]], i + 1, chain[i].c[0], chain[i].c[1], chain[i].c[2], errorest(econtrib[i]));
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
		    ++atomn, " O  ", rnames[seq[0][i]], i + 1, chain[i].o[0], chain[i].o[1], chain[i].o[2], errorest(econtrib[i]));
	    if (seq[0][i] != GLY)
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			++atomn, " CB ", rnames[seq[0][i]], i + 1, chain[i].cb[0], chain[i].cb[1], chain[i].cb[2], errorest(econtrib[i]));

        if(chain[i].cb[0] > max_x) max_x = chain[i].cb[0];
        if(chain[i].cb[1] > max_y) max_y = chain[i].cb[1];
        if(chain[i].cb[0] < min_x) min_x = chain[i].cb[0];
        if(chain[i].cb[1] < min_y) min_y = chain[i].cb[1];

	    if (chain[i].sc && seq[0][i] != GLY && seq[0][i] != ALA)
	    {
		for (j=1; j<nscatoms[seq[0][i]]; j++)
		    fprintf(ofp, "ATOM   %4d  %-3.3s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
			    ++atomn, chain[i].sc[j].atmnam, rnames[seq[0][i]], i + 1, chain[i].sc[j].pos[0], chain[i].sc[j].pos[1], chain[i].sc[j].pos[2], errorest(econtrib[i]));
	    }
	}

    // Add membrane
    if (envwt > 0.0){
        max_x += 8;
        max_y += 8;
        min_x -= 8;
        min_y -= 8;
        int count = 10000;
        for(i = 0; i <= (int)(max_x-min_x)/2;i++){
            for(j = 0; j <= (int)(max_y-min_y)/2;j++){
                double z = 15.000;
                double x = min_x + (i * 2);
                double y = min_y + (j * 2);
                fprintf(ofp, "HETATM%5d  O   DUM  0000    %8.3f%8.3f%8.3f\n", count,x,y,z);
                count++;
                if(count > 99999) count = 10000;
                z = -15.000;
                fprintf(ofp, "HETATM%5d  N   DUM  0000    %8.3f%8.3f%8.3f\n", count,x,y,z);
                count++;
                if(count > 99999) count = 10000;
            }
        }
    }

	fprintf(ofp, "END\n");

	fclose(ofp);
    }
    //exit(1);

    return;
}


/* Evaluate given set of coordinates */
float eval(RESATM * t, int reinit){

    int i;
    float e, mq, rmsd;
    static unsigned int ncalls;

    if (!reinit)
	   ncalls++;

    buildschn(t, 0, seqlen-1);

    e = e_model(dsqmat, t, FALSE);
    if (e < global_e_min - 0.01F && !reinit){

    	chaincpy(bestchn, t, seqlen);
    	global_e_min = e;

        /*
    	if (mqapcmd == NULL)
    	    e = e_model(dsqmat, t, TRUE);
    	else
    	    printf("RMSD = %f MQAP = %f  E = %f\n", rmsd, mq, e);
        */

    	if (verboseflg)
    	    printf("%f conformations per CPU second.\n", (float) ncalls / cputime());

        /*
    	if (fabsf(e - global_e_min) > 0.0001F * fabsf(e))
    	{
    	    printf("!!! %f %f\n", e, global_e_min);
    	    printf("%f\n", e_model(dsqmat, t, FALSE));
    	    printf("%f\n", e_model(dsqmat, t, FALSE));
    	    printf("%f\n", e_model(dsqmat, t, FALSE));
    	    printf("%f\n", e_model(dsqmat, t, FALSE));
    	    printf("%f\n", e_model(dsqmat, t, TRUE));
    	    fail("eval: inconsistent results!");
    	}
        */

    	if (verboseflg)
    	    writepdb(t, e);

    	if (rfp){

    	    fprintf(rfp, "zap\n");
    	    fflush(rfp);
    	    fprintf(rfp, "load %s\n", outpdbn);
    	    fflush(rfp);
    	    fprintf(rfp, "restrict *.CA\n");
    	    fflush(rfp);
    	    fprintf(rfp, "backbone 20\n");
    	    fflush(rfp);
    	}
    }
    return e;
}


/* Iteratively optimize side-chain packing */
void sc_opt(float (**dsqmat)[NPAIRS], RESATM * chn)
{
    int             i, nr, oldnr, optflg, pass=1;
    float e, e_min;

    buildschn(chn, 0, seqlen-1);

    printf("sc_opt : before = %f\n", e_min = calc_steric(dsqmat, chn, -1));

    do
    {
	for (optflg=i=0; i<seqlen; i++)
	    if (nrots[seq[0][i]] > 1)
	    {
		oldnr = chn[i].rotnum;
		for (nr=0; nr<nrots[seq[0][i]]; nr++)
		{
		    chn[i].rotnum = nr;
		    buildschn(chn, i, i);
		    e = calc_steric(dsqmat, chn, -1);
		    if (e < e_min - 0.001F)
		    {
			e_min = e;
			optflg = TRUE;
			oldnr = nr;
		    }
		}

		chn[i].rotnum = oldnr;
		buildschn(chn, i, i);
	    }
    } while (optflg && pass <= 10);

    printf("sc_opt : after = %f\n", e_min);
}


/* Ask Boltzmann if we should accept this energy change! */
int             boltzmann(float de, float t)
{
    if (t > 0.0F && de >= 0.0F)
	return ran0() < exp(-de / t);
    else
	return de < 0.0F;
}


/* Locate secondary structure elements */
int
                loc_sst(char *tsstruc, int len)
{
    int             i, istart, l, nsst = 0;

    for (i = 0; i < len; i++)
	if (tsstruc[i] != COIL)
	{
	    istart = i;
	    while (i < len && tsstruc[i] == tsstruc[istart])
		i++;
	    l = i - istart;
	    if (l < 3)
		continue;
	    sstlist[nsst].start = istart;
	    sstlist[nsst].length = l;
	    switch (tsstruc[istart])
	    {
	    case HELIX:
		sstlist[nsst].type = HELIX;
		break;
	    case STRAND:
		sstlist[nsst].type = STRAND;
		break;
	    default:
		printf("Code = %d\n", tsstruc[istart]);
		fail("loc_sst: unknown secondary structure code!");
	    }
	    nsst++;
	}

    return nsst;
}



/* Build supersecondary fragment lookup table */
void
                mkfragtb(void)
{
    int             i, j, k, l, m, n, ns, npair, nrms, len, sim, nclose, nfrag,
	nsst, from, to, ncon, ooi[MAXSEQLEN];
    float         **cb_d, pair, seqsim, tot, av_n = 0.0F, minwt, rmsd, rmsdtot,
	minrms, maxrms, av, sd;
    char           *type;
    
    if (verboseflg)
	puts("Building supersecondary motif lists...");
    
    for (i = 0; i < seqlen; i++)
    {
	fraglist[i].totbest = -1000.0F;
	fraglist[i].type = "";
	fraglist[i].frags = (chnidx*)malloc(MAXFRAGS * sizeof(struct chnidx));

	if (fraglist[i].frags == NULL)
	    fail("mkfragtb: cannot allocate fraglist!");
    }

    for (i = 0; i < nchn; i++)
    {
	cb_d = (float**)allocmat(chnlist[i].length, chnlist[i].length, sizeof(float));
	
	for (j = 0; j < chnlist[i].length; j++)
	    ooi[j] = 0;

	for (j = 0; j < chnlist[i].length; j++)
	    for (k = j + 1; k < chnlist[i].length; k++)
	    {
		cb_d[j][k] = cb_d[k][j] = dist(chnlist[i].chain[j].cb, chnlist[i].chain[k].cb);
		if (j-i > 1 && cb_d[j][k] <= OOICUTOFF)
		{
		    ooi[j]++;
		    ooi[k]++;
		}
	    }
	
	nsst = loc_sst(chnlist[i].sstruc, chnlist[i].length);
/*	printf("nsst = %d\n", nsst); */
	for (n = 0; n < nsst - 1; n++)
	{
	    from = to = -1;

	    if (sstlist[n].type == HELIX && sstlist[n + 1].type == HELIX)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 8.0F)
			    ncon++;
		if (ncon > 4)
		{
		    if (verboseflg)
			printf("Alpha hairpin at %d\n", sstlist[n].start);
		    type = "ALPHA HAIRPIN";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
		else
		{
		    if (verboseflg)
			printf("Alpha corner at %d\n", sstlist[n].start);
		    type = "ALPHA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
	    }
	    else if (sstlist[n].type == STRAND && sstlist[n + 1].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 3)
		{
		    if (verboseflg)
			printf("Beta hairpin at %d\n", sstlist[n].start);
		    type = "BETA HAIRPIN";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 0
		else
		{
		    if (verboseflg)
			printf("Beta corner at %d\n", sstlist[n].start);
		    type = "BETA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }
	    else if (n < nsst - 2 && sstlist[n].type == STRAND && sstlist[n + 1].type == HELIX && sstlist[n + 2].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 2].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 2].start + k] < 9.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("BAB unit at %d\n", sstlist[n].start);
		    type = "BAB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 2].start + sstlist[n + 2].length - 1;
		}
		else
		{
		    if (verboseflg)
			printf("Split BAB unit at %d\n", sstlist[n].start);
		    type = "Split BAB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 2].start + sstlist[n + 2].length - 1;
		}
	    }
	    else if (sstlist[n].type == STRAND && sstlist[n + 1].type == HELIX)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("BA unit at %d\n", sstlist[n].start);
		    type = "BA UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 0
		else
		{
		    if (verboseflg)
			printf("BA corner at %d\n", sstlist[n].start);
		    type = "BA CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }
	    else if (sstlist[n].type == HELIX && sstlist[n + 1].type == STRAND)
	    {
		ncon = 0;
		for (j = 0; j < sstlist[n].length; j++)
		    for (k = 0; k < sstlist[n + 1].length; k++)
			if (cb_d[sstlist[n].start + j][sstlist[n + 1].start + k] < 7.0F)
			    ncon++;
		if (ncon > 2)
		{
		    if (verboseflg)
			printf("AB unit at %d\n", sstlist[n].start);
		    type = "AB UNIT";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#if 0
		else
		{
		    if (verboseflg)
			printf("AB corner at %d\n", sstlist[n].start);
		    type = "AB CORNER";
		    from = sstlist[n].start;
		    to = sstlist[n + 1].start + sstlist[n + 1].length - 1;
		}
#endif
	    }

	    if (from >= 0)
	    {
		len = to - from + 1;
		if (len > 1)
		    for (k = 0; k <= seqlen - len; k++)
		    {
			/* Only pick fragments within SVRs */
			if (mod_mode == MODELM)
			{
			    for (l = 1; l < len; l++)
				if (!svrflg[k+l])
				    break;

			    if (l != len)
				continue;
			}
			
			for (sim = l = 0; l < len; l++)
			    if (tplt_ss[k + l] & chnlist[i].sstruc[from + l])
				sim++;

			if (sim < len-1)
			    continue;
			
			for (rmsd = pair = seqsim = 0.0F, npair = nrms = l = 0; l < len; l++)
			{
			    seqsim += dmmmat[chnlist[i].seq[from + l]][seq[0][k + l]];
			    
			    for (m = l + 1; m < len; m++)
			    {
				for (ns = 0; ns < n_ds; ns++)
				    if (ds_from[ns] == k+l && ds_to[ns] == k+m && cb_d[from+l][from+m] > 6.0F)
					pair -= 1000.0F;

				if (cb_d[from + l][from + m] < 10.0F)
				{
				    pair -= conmat[k + l][k + m];
				    npair++;
				}
			    }
			}
				
/*		    tot = ran0(); */
/*		    tot = pair; */
			
		    tot = 0.1F * seqsim / len + pair / MAX(1, npair);

		    /* Make use of target structure for selecting fragments if modelling */
		    if (fragsel && nrms > 1)
			tot = -rmsd / nrms;

#if 0
		    for (l = 0; l < len; l++)
			printf("%c %f - %c %f\n", rescodes[seq[0][k+l]], bestchn[k+l].cb[0], rescodes[chnlist[i].seq[from+l]], chnlist[i].chain[from+l].cb[0]);

		    printf("%d %d %d %d %d %f\n", k, len, seqlen, from, to, sqrt(-tot));
#endif

#if 0
		    printf("%g %g %g\n", seqsim, pair, tot);
#endif
		    
		    if (tot > fraglist[k].totbest)
		    {
			fraglist[k].totbest = tot;
			fraglist[k].bestlen = len;
			fraglist[k].nchn = i;
			fraglist[k].from = from;
			fraglist[k].type = type;
			fraglist[k].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
		    }

		    if (fraglist[k].nfrags < MAXFRAGS)
		    {
			fraglist[k].frags[fraglist[k].nfrags].chain = i;
			fraglist[k].frags[fraglist[k].nfrags].pos = from;
			fraglist[k].frags[fraglist[k].nfrags].weight = tot;
			fraglist[k].frags[fraglist[k].nfrags].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
			fraglist[k].frags[fraglist[k].nfrags].length = len;
			fraglist[k].totwt += tot;
			fraglist[k].totsum += tot;
			fraglist[k].totsumsq += SQR(tot);
			fraglist[k].nfrags++;
			fraglist[k].ntot++;
		    }
		    else
		    {
			minwt = fraglist[k].frags[0].weight;
			j = 0;
			for (l = 1; l < MAXFRAGS; l++)
			    if (fraglist[k].frags[l].weight < minwt)
			    {
				minwt = fraglist[k].frags[l].weight;
				j = l;
			    }
			if (minwt < tot)
			{
			    fraglist[k].frags[j].chain = i;
			    fraglist[k].frags[j].pos = from;
			    fraglist[k].frags[j].length = len;
			    fraglist[k].totwt -= minwt;
			    fraglist[k].frags[j].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 0.0F;
			    fraglist[k].frags[j].weight = tot;
			    fraglist[k].totwt += tot;
			    fraglist[k].totsum += tot;
			    fraglist[k].totsumsq += SQR(tot);
			    fraglist[k].ntot++;
			}
		    }
		}
	    }
	}
	freemat(cb_d, chnlist[i].length);
    }

    for (rmsdtot = nrms = av_n = nfrag = nclose = i = 0; i < seqlen; i++)
	if (fraglist[i].nfrags)
	{
	    for (len = j = 0; j < fraglist[i].nfrags; j++)
		len = MAX(fraglist[i].frags[j].length, len);
	    av = fraglist[i].totsum / fraglist[i].ntot;
	    sd = sqrt(fraglist[i].totsumsq / fraglist[i].ntot - SQR(av));
	    if (verboseflg)
		printf("%4d %2d %4d %4d %4d %8.3f %8.3f %5.2f %-14s ", i, fraglist[i].bestlen, fraglist[i].nchn, fraglist[i].from, fraglist[i].nfrags, fraglist[i].totbest, (fraglist[i].totbest - av) / sd, fraglist[i].rmsd, fraglist[i].type);
	    minrms = 1000.0F;
	    maxrms = 0.0F;
	    for (j = 0; j < fraglist[i].nfrags; j++)
	    {
		minrms = MIN(minrms, fraglist[i].frags[j].rmsd);
		maxrms = MAX(maxrms, fraglist[i].frags[j].rmsd);
		if (fraglist[i].frags[j].rmsd < 2.0F)
		    nclose++;
		nfrag++;
	    }
	    if (verboseflg)
	    {
		printf(" %4.1f %4.1f", minrms, maxrms);
		putchar('\n');
	    }
	    av_n += fraglist[i].nfrags;
	    rmsdtot += minrms;
	    nrms++;
	}

    if (verboseflg)
    {
	printf("Average Min RMSD = %f\n", rmsdtot / nrms);
	printf("p(RMSD < 2.0) = %f\n", (float) nclose / nfrag);
	printf("Average frags per position = %f\n", av_n / seqlen);
    }
}

/* Build fixed length fragment lookup table */
void
                mkfragtb2(void)
{
    int             i, j, k, l, m, ns, npair, nrms, len, fraglen, from, to, ooi[MAXSEQLEN], nfrag, nclose;
    float         **cb_d, pair, seqsim, tot, av_n = 0.0F, minwt, rmsd, rmsdtot, minrms, maxrms, av, sd;
    char           *type;
    Transform       xf1, xf2;
    Vector point1, point2;

    if (verboseflg)
	puts("Building fixed length fragment lists...");

    for (i = 0; i < seqlen; i++)
    {
	fraglist2[i].totbest = -1000.0F;
	fraglist2[i].type = "";
	fraglist2[i].frags = (chnidx*)malloc(MAXFRAGS2 * sizeof(struct chnidx));

	if (fraglist2[i].frags == NULL)
	    fail("mkfragtb: cannot allocate fraglist2!");
    }

    for (i = 0; i < nchn; i++)
    {
	cb_d = (float**)allocmat(chnlist[i].length, chnlist[i].length, sizeof(float));
	
	for (j = 0; j < chnlist[i].length; j++)
	    ooi[j] = 0;

	for (j = 0; j < chnlist[i].length; j++)
	    for (k = j + 1; k < chnlist[i].length; k++)
	    {
		cb_d[j][k] = cb_d[k][j] = dist(chnlist[i].chain[j].cb, chnlist[i].chain[k].cb);
		if (j-i > 1 && cb_d[j][k] <= OOICUTOFF)
		{
		    ooi[j]++;
		    ooi[k]++;
		}
	    }

	for (fraglen=MINFRAGS2LEN; fraglen<=MAXFRAGS2LEN; fraglen++)
	    for (from = 0; from <= chnlist[i].length - fraglen; from++)
	    {
		to = from + fraglen - 1;
		len = to - from + 1;
		
		type = (char*)malloc(len + 1);
		if (type == NULL)
		    fail("mkfragtb2: malloc failed!");
		
		for (k = 0; k < len; k++)
		    type[k] = sscodes[chnlist[i].sstruc[from + k]];
		type[len] = '\0';
		
		if (from >= 0)
		{
		    if (len > 1)
			for (k = 0; k <= seqlen - len; k++)
			{
			    /* Only pick fragments within SVRs */
			    if (mod_mode == MODELM)
			    {
				for (l = 1; l < len; l++)
				    if (!svrflg[k+l])
					break;
				
				if (l != len)
				    continue;
			    }
#if 1
			    for (l = 0; l < len; l++)
				if (!(tplt_ss[k + l] & chnlist[i].sstruc[from + l]))
				    break;
			    
			    if (l != len)
				continue;
#endif			
			    for (rmsd = pair = seqsim = 0.0F, npair = nrms = l = 0; l < len; l++)
			    {
				seqsim += dmmmat[chnlist[i].seq[from + l]][seq[0][k + l]];
				
				for (m = l + 1; m < len; m++)
				{
				    for (ns = 0; ns < n_ds; ns++)
					if (ds_from[ns] == k+l && ds_to[ns] == k+m && cb_d[from+l][from+m] > 6.0F)
					    pair -= 1000.0F;
				    
				    if (cb_d[from + l][from + m] < 10.0F)
				    {
					pair -= conmat[k + l][k + m];
					npair++;
				    }
				}
			    }
			    
			    
/*		    tot = ran0(); */
/*		    tot = pair; */
			    
			    tot = 0.05F * seqsim / len + pair / MAX(1, npair);

			    /* Try to use target structure for selecting fragments if modelling or refining */
			    if (mod_mode == REFINEM && nrms > 1 && !svrflg[k] && !svrflg[k+len-1] && segidx[k] == segidx[k+len-1])
			    {
				    /* Find fragments which will cause small global RMSD changes */
				    calcxf(targchn[k].n, targchn[k].ca, targchn[k].c, xf1);
				    calcxf(chnlist[i].chain[from].n, chnlist[i].chain[from].ca, chnlist[i].chain[from].c, xf2);
				    
				    /* See where fragment will "land" when spliced */
				    transform_point(xf1, targchn[k+len-1].n, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].n, point2);
				    tot = -distsq(point1, point2);
				    transform_point(xf1, targchn[k+len-1].ca, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].ca, point2);
				    tot -= distsq(point1, point2);
				    transform_point(xf1, targchn[k+len-1].c, point1);
				    transform_point(xf2, chnlist[i].chain[from+len-1].c, point2);
				    tot -= distsq(point1, point2);
			    }
			    else if (fragsel && nrms > 1)
				tot = -rmsd / nrms;

#if 0
			    printf("%g %g\n", pair, tot);
#endif
			    
			    if (tot > fraglist2[k].totbest)
			    {
				fraglist2[k].totbest = tot;
				fraglist2[k].bestlen = len;
				fraglist2[k].nchn = i;
				fraglist2[k].from = from;
				fraglist2[k].type = type;
				fraglist2[k].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
			    }
			    
			    if (fraglist2[k].nfrags < MAXFRAGS2)
			    {
				fraglist2[k].frags[fraglist2[k].nfrags].chain = i;
				fraglist2[k].frags[fraglist2[k].nfrags].pos = from;
				fraglist2[k].frags[fraglist2[k].nfrags].weight = tot;
				fraglist2[k].frags[fraglist2[k].nfrags].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
				fraglist2[k].frags[fraglist2[k].nfrags].length = len;
				fraglist2[k].totwt += tot;
				fraglist2[k].totsum += tot;
				fraglist2[k].totsumsq += SQR(tot);
				fraglist2[k].nfrags++;
				fraglist2[k].ntot++;
			    }
			    else
			    {
				minwt = fraglist2[k].frags[0].weight;
				j = 0;
				for (l = 1; l < MAXFRAGS2; l++)
				    if (fraglist2[k].frags[l].weight < minwt)
				    {
					minwt = fraglist2[k].frags[l].weight;
					j = l;
				    }
				if (minwt < tot)
				{
				    fraglist2[k].frags[j].chain = i;
				    fraglist2[k].frags[j].pos = from;
				    fraglist2[k].frags[j].length = len;
				    fraglist2[k].totwt -= minwt;
				    fraglist2[k].frags[j].rmsd = rmsd > 0.0F ? sqrt(rmsd / nrms) : 999.0F;
				    fraglist2[k].frags[j].weight = tot;
				    fraglist2[k].totwt += tot;
				    fraglist2[k].totsum += tot;
				    fraglist2[k].totsumsq += SQR(tot);
				    fraglist2[k].ntot++;
				}
			    }
			}
		}
	    }
	freemat(cb_d, chnlist[i].length);
    }
    
    for (rmsdtot = nrms = av_n = nfrag = nclose = i = 0; i < seqlen; i++)
	if (fraglist2[i].nfrags)
	{
	    for (len = j = 0; j < fraglist2[i].nfrags; j++)
		len = MAX(fraglist2[i].frags[j].length, len);
	    av = fraglist2[i].totsum / fraglist2[i].ntot;
	    sd = sqrt(fraglist2[i].totsumsq / fraglist2[i].ntot - SQR(av));
	    if (verboseflg)
		printf("%4d %2d %4d %4d %4d %8.3f %8.3f %5.2f %-17s ", i, fraglist2[i].bestlen, fraglist2[i].nchn, fraglist2[i].from, fraglist2[i].nfrags, fraglist2[i].totbest, (fraglist2[i].totbest - av) / sd, fraglist2[i].rmsd, fraglist2[i].type);
	    minrms = 1000.0F;
	    maxrms = 0.0F;
	    for (j = 0; j < fraglist2[i].nfrags; j++)
	    {
		minrms = MIN(minrms, fraglist2[i].frags[j].rmsd);
		maxrms = MAX(maxrms, fraglist2[i].frags[j].rmsd);
	    }

	    if (minrms < 1.0F)
		nclose++;
	    nfrag++;

	    if (verboseflg)
	    {
		printf(" %4.1f %4.1f", minrms, maxrms);
		putchar('\n');
	    }
	    av_n += fraglist2[i].nfrags;
	    rmsdtot += minrms;
	    nrms++;
	}
    
    if (verboseflg)
    {
	printf("Average Min RMSD = %f\n", rmsdtot / nrms);
	printf("p(RMSD < 1.0) = %f\n", (float) nclose / nfrag);
	printf("Average frags per position = %f\n", av_n / seqlen);
    }
}

unsigned int movefreq[3] = { 5000, 5000, 5000 };
unsigned int lastmove;
float moveprob[3] = { 0.333333, 0.333333, 0.333333 };

/* Choose large fragment from supersec. or fixed-length list */
void fragsamp(int *to, int *src, int *from, int *len)
{
    int             fn = -1;

    do
    {
	*to = randint(0 + last_boundary, MIN(curr_seqlen,seqlen)-1);

	if (MAXFRAGS > 0 && ran0() < moveprob[0] / (moveprob[0] + moveprob[1]))
	{
	    if (fraglist[*to].nfrags)
	    {
		fn = randint(0, fraglist[*to].nfrags - 1);
		*len = fraglist[*to].frags[fn].length;
		*src = fraglist[*to].frags[fn].chain;
		*from = fraglist[*to].frags[fn].pos;
		lastmove = 0;
	    }
	}
	else
	{
	    if (fraglist2[*to].nfrags)
	    {
		fn = randint(0, fraglist2[*to].nfrags - 1);
		*len = fraglist2[*to].frags[fn].length;
		*src = fraglist2[*to].frags[fn].chain;
		*from = fraglist2[*to].frags[fn].pos;
		lastmove = 1;
	    }
	}
    }
    while (fn < 0);
}


void sa_randmove(){
    
    int             i, len, src, from, to;

    do {
	if (ran0() < moveprob[0] + moveprob[1])
	{
	    fragsamp(&to, &src, &from, &len);
        //printf("1 len %d to %d\n",len,to);
	}
	else
	{
	    do
	    {
		len = randint(2, 3);
		to = randint(0 + last_boundary, MIN(curr_seqlen,seqlen)-len);
		src = randint(0, nchn - 1);
		from = randint(0, chnlist[src].length - len);



	    }
	    while (!chksst(src, from, to, len));
	    lastmove = 2;

        //printf("2 len %d to %d\n",len,to);
	}

	/* Only pick fragments within SVRs */
	if (mod_mode == MODELM)
	{
	    for (i = 0; i < len; i++)
		if (!svrflg[to+i] || (i > 0 && svrflg[to+i] != svrflg[to+i-1]))
		    break;
	    
	    if (i == len)
		break;
	}
    } while (mod_mode == MODELM);

    //printf("3 len %d to %d\n",len,to);
    //printf("In sa_randmove, calling splice with %d %d\n",to, len);
    splice(curchn, chnlist[src].chain + from, to, len);
}


/* Replica Exchange Monte Carlo */
void repmc(void){

    float           ediff, av_ediff = 0.0, max_ediff = 0.0, old_e, new_e, tmax, chi0,
	e_min, delta, rmsd, lastrmsd, *rep_e, etemp, ttemp;
    int             i, j, cycle=0, *nsucc, *nuphill, nswaps, src, minsteps;
    Replica         *ensemble;
    static unsigned optcount[100], swapcount[100];
    FILE *ifp;

    curr_helices = helices;
    curr_seqlen = seqlen;
    if(SEQHELIX){
        curr_helices = MIN(2,helices);
        targrange = seqlen;
    }else{
        targrange = 6;
    }    

    fullyextended(curchn);

    /* Generate starting conformation */
    std::cout << "Generate starting conformation" << std::endl;
    for (i = 0; i < 500; i++){

	   chaincpy(oldchn, curchn, seqlen);
	   for (;;){
	       sa_randmove();
            //pre_position(curchn);
            //if(!calcdist(curchn, TRUE) && !mean_tm_ztest(curchn)){
            //    cout << i << " calcdist and mean_tm_ztest OK" << endl;
            //    break;
            //}    
	       if (!calcdist(curchn, TRUE) && (!zfname[0] || !ztest(dsqmat)))
		      break;
	       chaincpy(curchn, oldchn, seqlen);
	   }  
    }
    std::cout << "done" << std::endl; 

    /* Evaluate average energy change for random moves */
    for (i = 0; i < 500; i++)
    {
	chaincpy(oldchn, curchn, seqlen);
	for (;;)
	{
	    sa_randmove();
        //pre_position(curchn);
        //if(!calcdist(curchn, TRUE) && !mean_tm_ztest(curchn)){
        //    break;
        //}    
	    if (!calcdist(curchn, TRUE) && (!zfname[0] || !ztest(dsqmat)))
		  break;
	    chaincpy(curchn, oldchn, seqlen);
	}

	new_e = eval(curchn, FALSE);

	if (i)
	{
	    ediff = fabs(new_e - old_e);

//	    printf("ediff = %f\n", ediff);
	    
	    av_ediff += fabs(new_e);
	    max_ediff = MAX(max_ediff, ediff);
	}

	old_e = new_e;
    }

    av_ediff /= (i+1.0);

    ensemble = (Replica *) calloc(poolsize, sizeof(Replica));
    nsucc = (int *) calloc(poolsize, sizeof(int));
    nuphill = (int *) calloc(poolsize, sizeof(int));

    if (mod_mode == 1)
	chaincpy(curchn, targchn, seqlen);

    if (verboseflg)
	puts("Replica Exchange Annealing...");

    minsteps = 0;
    
    printf("Calc global min....\n");
    mem_debug = 1;
    global_e_min = e_min = eval(curchn, FALSE);
    printf("done....\n");
    //exit(1);

    for (i=0; i<poolsize; i++)
    {
	ensemble[i].chn = (RESATM *) calloc(seqlen, sizeof(RESATM));

	chaincpy(ensemble[i].chn, curchn, seqlen);

	ensemble[i].energy = e_min;
//	ensemble[i].t = av_ediff * (poolsize - i) / poolsize;
	ensemble[i].t = pow(TRATIO, (float)i) * av_ediff * INITEMP;
//	ensemble[i].t = pow(TRATIO, (float)i) * max_ediff * INITEMP;
    }

    if(SEQHELIX && helices > 1){
        curr_seqlen = 1+topology[3]; 
        printf("\n\nFolding helix-by-helix. Starting length is %d\n",curr_seqlen);   
    }

    printf("Begin....\n\n");

    while (totalsteps < MAXSTEPS)
    {

    if(!SEQHELIX){    
    	if (minsteps > MAXSTEPS / seqlen / 2){
    	    printf("New Target Range = %d\n", targrange = 6 + seqlen * 1.1 * totalsteps / MAXSTEPS);
    //	    printf("New Target Range = %d\n", ++targrange);
    	    minsteps = 0;
    	}
    }else{
        if(curr_helices < helices && totalsteps >= (curr_helices*((0.5*MAXSTEPS)/(helices-2)))){
            curr_helices++;
            curr_seqlen = 1+topology[(curr_helices*2)-1];
            last_boundary = 1+topology[(curr_helices*2)-3];
            printf("\n***** Folding helices 1 - %d - Sequence length = %d\n", curr_helices,curr_seqlen); 
            e_min = eval(curchn, FALSE);
            z_rot_best = 0, x_rot_best = 0, y_rot_best = 0, z_trans_best = 0;
            writepdb(curchn, e_min);  
            minsteps = 0;     
        }else if(curr_helices == helices){
            curr_seqlen = seqlen; 
            last_boundary = 0;    
            if (minsteps > MAXSTEPS / seqlen / 4){
                printf("New Target Range = %d\n", targrange = 6 + seqlen * 1.1 * totalsteps / MAXSTEPS);
                minsteps = 0;
            }
        }
    }    


	for (i=0; i<poolsize; i++)
	{
	    chaincpy(curchn, ensemble[i].chn, seqlen);

	    for (;;){
    		sa_randmove();
            //pre_position(curchn);
            //if(!calcdist(curchn, TRUE) && !mean_tm_ztest(curchn)){
            //    break;
            //}    
    		if (!calcdist(curchn, TRUE) && (!zfname[0] || !ztest(dsqmat)))
    		    break;
    		chaincpy(curchn, ensemble[i].chn, seqlen);
	    }

//	    sc_randmove();
	    	
	    totalsteps++;

	    new_e = eval(curchn, FALSE);

	    if (new_e > ensemble[i].energy)
		nuphill[i]++;

	    minsteps++;

	    if (boltzmann(new_e - ensemble[i].energy, ensemble[i].t))
	    {
		if (new_e < e_min)
		{
		    e_min = new_e;
		    optcount[i]++;
		    minsteps = 0;
		}

		if (new_e > ensemble[i].energy)
		    nsucc[i]++;

		/* printf("new_e : %d %f\n", nsucc+1, new_e); */
		ensemble[i].energy = new_e;
		chaincpy(ensemble[i].chn, curchn, seqlen);
		movefreq[lastmove]++;
	    }
	}

	/* Replica exchange... */
	nswaps = 0;

//	i = randint(0, poolsize-2);
//	j = randint(i+1, poolsize-1);

	for (i=0; i<poolsize; i++){
	    for (j=i+1; j<poolsize; j++){
		  delta = (1.0/ensemble[i].t - 1.0/ensemble[j].t) * (ensemble[i].energy - ensemble[j].energy);
		  if (delta >= 0.0 || ran0() <= exp(delta)){
		    ttemp = ensemble[i].t;
		    ensemble[i].t = ensemble[j].t;
		    ensemble[j].t = ttemp;
		    nswaps++;
		    swapcount[i]++;
		    swapcount[j]++;
		  }
	    }
	}
	
	if (verboseflg && !(++cycle % 100))
	{
	    printf("\nEner:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", ensemble[i].energy);
	    printf("\nTemp:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", ensemble[i].t);
	    printf("\nchi0:");
	    for (i=0; i<poolsize; i++)
		printf(" %f", (float) nsucc[i] / (float) nuphill[i]);
	    printf("\nOptC:");
	    for (i=0; i<poolsize; i++)
		printf(" %d", optcount[i]);
	    printf("\nSwap:");
	    for (i=0; i<poolsize; i++)
		printf(" %d", swapcount[i]);
	    printf("\nCycle %3d : Swaps = %d Steps = %d Helices = %d Min = %d Length = %d Emin = %g\n", cycle, nswaps, totalsteps, curr_helices, last_boundary, MIN(curr_seqlen,seqlen), e_min);
	}

	
	/* Calculate current move accept probabilities */

#if 0
	moveprob[0] = (float) movefreq[0] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[1] = (float) movefreq[1] / (movefreq[0] + movefreq[1] + movefreq[2]);
	moveprob[2] = (float) movefreq[2] / (movefreq[0] + movefreq[1] + movefreq[2]);

/*	printf("New move probabilities: %f %f %f\n", moveprob[0], moveprob[1], moveprob[2]); */
#endif

	fflush(stdout);
    }
}

void
                readxyz(char *buf, float *x, float *y, float *z)
{
    char            temp[9];

    temp[8] = '\0';
    strncpy(temp, buf, 8);
    *x = atof(temp);
    strncpy(temp, buf + 8, 8);
    *y = atof(temp);
    strncpy(temp, buf + 16, 8);
    *z = atof(temp);
}

/* Read target structure coords */
void
                maketemplate(char *pdbname, char chainid, int model, RESATM *chn)
{
    FILE           *pfp;
    char            buf[160], whichatm[MAXSEQLEN], inscode;
    int             atc, i, j, k, nres, resid, lastid, nseg;
    float           dv, maxd, sx, sy;
    Vector          cca, nca, xx, yy;
    short           at1, at2;
    float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
    static const float caca_maxdist[5] = { 3.90, 7.34, 10.81, 14.19, 17.58 };

    if (verboseflg)
	printf("Reading %s chain %c ...\n", pdbname, chainid);
    pfp = fopen(pdbname, "r");
    if (!pfp)
	fail("maketemplate: Cannot open PDB file!");

    for (i = 0; i < MAXSEQLEN; i++)
	whichatm[i] = 0;

    if (model >= 1)
    {
	while (!feof(pfp))
	{
	    if (fgets(buf, 160, pfp) == NULL)
		fail("maketemplate: cannot find MODEL!");
	    if (strncmp(buf, "MODEL", 5))
		continue;
	    sscanf(buf + 5, "%d", &i);
	    if (i == model)
		break;
	}
    }

    i = -1;
    nsvr = nseg = 0;
    while (!feof(pfp))
    {
	if (fgets(buf, 160, pfp) == NULL)
	    break;
	if (!strncmp(buf, "END", 3))
	    break;
	if (!strncmp(buf, "TER", 3))
	{
	    nseg++;
	    continue;
	}
	/* printf("%d %s\n",i,buf); */
	if (strncmp(buf, "ATOM", 4) || (buf[21] != chainid && !(buf[21] == 'A' && chainid == ' ')) || (buf[16] != ' ' && buf[16] != 'A'))
	    continue;
	for (atc = 0; atc <= CATOM; atc++)
	    if (!strncmp(buf + 13, atmnames[atc], 2))
	    {
		inscode = buf[26];
		buf[26] = ' ';
		sscanf(buf + 22, "%d", &resid);
		if (atc == NATOM)
		{
		    ++i;
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		    if (strlen(buf) > 56)
			svrflg[i] = buf[56] != '1';
		    else
			svrflg[i] = 0;
		    if (svrflg[i])
			nsvr++;
		    segidx[i] = nseg;
		    if (!i)
			firstid = resid;
		    else if (inscode == ' ' && resid - lastid > 1)
		    {
			if (SQR(x[i][NATOM] - x[i - 1][NATOM]) + SQR(y[i][NATOM] - y[i - 1][NATOM]) + SQR(z[i][NATOM] - z[i - 1][NATOM]) > 15.0)
			{
			    printf("WARNING: Skipping %d missing residues!\n", resid - lastid - 1);
			    for (k = 0; k < resid - lastid - 1; k++)
				i++;
			}
		    }
		    lastid = resid;
		}
		else if (i >= 0)
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		whichatm[i] |= 1 << atc;
	    }
    }
    fclose(pfp);

    nres = i + 1;

    /* Check atoms */
    for (i = 0; i < nres; i++)
	if (!svrflg[i])
	{
	    if (!(whichatm[i] & (1 << NATOM)))
	    {
		printf("FATAL: Missing N atom in %d!\n", i + 1);
		exit(1);
	    }
	    if (!(whichatm[i] & (1 << CAATOM)))
	    {
		printf("WARNING: Missing CA atom in %d!\n", i + 1);
	    }
	    if (!(whichatm[i] & (1 << CBATOM)))
	    {
		if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)))
		{
		    /* Not much left of this residue! */
		    printf("WARNING: Missing main-chain atom in %d!\n", i + 1);
		    continue;
		}
		
		/* Reconstruct CB atom */
		nca[0] = x[i][CAATOM] - x[i][NATOM];
		nca[1] = y[i][CAATOM] - y[i][NATOM];
		nca[2] = z[i][CAATOM] - z[i][NATOM];
		cca[0] = x[i][CAATOM] - x[i][CATOM];
		cca[1] = y[i][CAATOM] - y[i][CATOM];
		cca[2] = z[i][CAATOM] - z[i][CATOM];
		vecadd(xx, nca, cca);
		vecprod(yy, nca, cca);
		sx = CACBDIST * cosf(TETH_ANG) / sqrtf(dotprod(xx, xx));
		sy = CACBDIST * sinf(TETH_ANG) / sqrtf(dotprod(yy, yy));
		x[i][CBATOM] = x[i][CAATOM] + xx[0] * sx + yy[0] * sy;
		y[i][CBATOM] = y[i][CAATOM] + xx[1] * sx + yy[1] * sy;
		z[i][CBATOM] = z[i][CAATOM] + xx[2] * sx + yy[2] * sy;
		whichatm[i] |= 1 << CBATOM;
	    }
	}

    if (nres != seqlen)
	fail("Sequence length mismatch in PDB file!");
    
    for (i = 0; i < nres; i++)
    {
	chn[i].n[0] = x[i][NATOM];
	chn[i].n[1] = y[i][NATOM];
	chn[i].n[2] = z[i][NATOM];
	chn[i].ca[0] = x[i][CAATOM];
	chn[i].ca[1] = y[i][CAATOM];
	chn[i].ca[2] = z[i][CAATOM];
	chn[i].c[0] = x[i][CATOM];
	chn[i].c[1] = y[i][CATOM];
	chn[i].c[2] = z[i][CATOM];
	chn[i].o[0] = x[i][OATOM];
	chn[i].o[1] = y[i][OATOM];
	chn[i].o[2] = z[i][OATOM];
	chn[i].cb[0] = x[i][CBATOM];
	chn[i].cb[1] = y[i][CBATOM];
	chn[i].cb[2] = z[i][CBATOM];
    }

    for (i = 0; i < nres; i++)
	chn[i].rotnum = 0;

    buildschn(chn, 0, seqlen-1);

    if (verboseflg)
	printf("Target NRES = %d\n", nres);
    
    /* Calculate interatomic distance templates */
    at1 = at2 = CBATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = sqrtf(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
	    dsqmat[i][j][CB_CB] = dsqmat[j][i][CB_CB] = dv;
	}

    at1 = at2 = CAATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = sqrtf(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
	    dsqmat[i][j][CA_CA] = dsqmat[j][i][CA_CA] = dv;
	}

    for (i = 0; i < nres; i++)
	tcbooi[i] = 0;

    /* Compute CB-Ooi numbers */
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	    if (dsqmat[i][j][CB_CB] > 0.0F && dsqmat[i][j][CB_CB] <= OOICUTOFF)
	    {
		tcbooi[i]++;
		tcbooi[j]++;
	    }
}


char *getnextp(FILE *pfp)
{
    static char pv[256];

    while (!feof(pfp))
    {
	if (fscanf(pfp, "%s", pv) != 1)
	    return NULL;
	if (pv[0] == '#')
	{
	    fgets(pv, 256, pfp);
	    continue;
	}

	if (verboseflg)
	    printf("Read parameter: %s\n", pv);
	return pv;
    }

    return NULL;
}

/* Read parameter file if present */
void readparams(char *fname)
{
    char *vname, *value;
    FILE *pfp;

    pfp = fopen(fname, "r");

    if (!pfp)
	fail("readparams: cannot open parameter file!");

    if (verboseflg)
	puts("\nReading parameter file...");

    while (vname = getnextp(pfp))
    {
	value = "N/A";

	if (!strcasecmp(vname, "ALNFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(alnfname, value);
	    if (verboseflg)
		puts(value);
	}
	else if (!strcasecmp(vname, "CONFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(confname, value);
	    if (verboseflg)
		puts(value);
	}
	else if (!strcasecmp(vname, "ZFILE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		strcpy(zfname, value);
	    if (verboseflg)
		puts(value);
	}
    else if (!strcasecmp(vname, "LIPEXFILE"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        strcpy(lipexfname, value);
        if (verboseflg)
        puts(value);
    }
	else if (!strcasecmp(vname, "HBWT"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		hbwt = atof(value);
	}
	else if (!strcasecmp(vname, "MAXFRAGS"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXFRAGS2"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS2 = atoi(value);
	}
	else if (!strcasecmp(vname, "MINFRAGS2LEN"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MINFRAGS2LEN = atoi(value);
	}
	else if (!strcasecmp(vname, "MAXFRAGS2LEN"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXFRAGS2LEN = atoi(value);
	}
	else if (!strcasecmp(vname, "POOLSIZE"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		poolsize = atoi(value);
	}
	else if (!strcasecmp(vname, "TRATIO"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		TRATIO = atof(value);
	}
	else if (!strcasecmp(vname, "MAXSTEPS"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		MAXSTEPS = atoi(value);
	}
	else if (!strcasecmp(vname, "INITEMP"))
	{
	    value = getnextp(pfp);
	    if (value != NULL)
		INITEMP = atof(value);
	}
    else if (!strcasecmp(vname, "THREADS"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        THREADS = atoi(value);
    }    
    else if (!strcasecmp(vname, "MAXGACALLS"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        MAXGACALLS = atoi(value);
    }    
    else if (!strcasecmp(vname, "GAPOOLSIZE"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        GAPOOLSIZE = atoi(value);
    }    
    else if (!strcasecmp(vname, "ENVWT"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        envwt = atof(value);
    }  
    else if (!strcasecmp(vname, "LIPEXWT"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        lipexwt = atof(value);
    }  
    else if (!strcasecmp(vname, "TOPWT"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        topwt = atof(value);
    }  
    else if (!strcasecmp(vname, "FOLDBYHELIX"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        SEQHELIX = atoi(value);
    }
    else if (!strcasecmp(vname, "ORIENTGA"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        ORIENTGA = atoi(value);
    }
    else if (!strcasecmp(vname, "SKIPZTEST"))
    {
        value = getnextp(pfp);
        if (value != NULL)
        SKIPZTEST = atoi(value);
    }
	else
	    fail("Unknown keyword (%s) in parameter file!", vname);

	if (value == NULL)
	{
	    fprintf(stderr, "%s - missing parameter value in parameter file!\n", vname);
	    exit(-1);
	}
    }

    fclose(pfp);
}


int main(int argc, char **argv)
{
    int             a, aa, b, i, j, k, optstart = 0;
    float eqd, dsd, maxd, prob, e;
    const int       bigendv=1;
    FILE           *ifp;
    
    printf("FILM3 Version %s (Last edit %s)\n", FILMVersion, Last_Edit_Date);
    printf("Motif-constrained Transmembrane Protein Folding Program\n");
    printf("Build date : %s\n",__DATE__);
    printf("Copyright (C) 1995 David T. Jones\nCopyright (C) 2011 University College London\nModified by Tim Nugent 2013\n\n");

    if (argc < 2)
	fail("usage : %s paramfile [outpdbfn] [pdbfname] [chainid] [options]\n", argv[0]);

    randomise();

    /* Parse options at end of command line */

    for (i=2; i<argc; i++)
    {
	if (argv[i][0] == '-')
	{
	    if (!optstart)
		optstart = i;
	    switch(argv[i][1])
	    {
	    case 'q':
		verboseflg = FALSE;
		break;
	    }
	}
    }

    if (optstart)
	argc = optstart;

    readparams(argv[1]);

    ifp = fopen("atomdist.dat", "rb");
    if (!ifp)
	fail("main: cannot open atomdist.dat!");
    fread(atomdistsq, 1, sizeof(atomdistsq), ifp);
    fclose(ifp);

    /* Swap bytes if running on bigendian machine */
    if (*(char *)&bigendv != 1)
	byteswap(atomdistsq, sizeof(atomdistsq)/4);

    for (i = 0; i < TOPOMAX; i++)
	for (a = 0; a < 4; a++)
	    for (b = 0; b < 4; b++)
		diffdist[i][a][b] = (maxdist[i][a][b] - mindist[i][a][b]);

    for (i=0; i<167; i++)
	for (j=0; j<167; j++)
	    for (k=0; k<3; k++)
		atomdistsq[i][j][k] = SQR(atomdistsq[i][j][k]);

    cout << "Reading aln file...\n";
    ifp = fopen(alnfname, "r");
    if (!ifp)
	fail("main: unable to open aln file!");
    fgets(buf, sizeof(buf), ifp);
    fscanf(ifp, "%d\n", &nseqs);
    fgets(tplt_ss, sizeof(tplt_ss), ifp);
    seqlen = strlen(tplt_ss) - 1;
    cout << "done.\n";

    if (seqlen < 10)
	fail("Target sequence length < 10!!");

    for (i=0; i<seqlen; i++)
    {
	if (tplt_ss[i] == 'C')
	    tplt_ss[i] = COIL;
	else if (tplt_ss[i] == 'H')
	    tplt_ss[i] = HELIX;
	else if (tplt_ss[i] == 'E')
	    tplt_ss[i] = STRAND;
	else if (tplt_ss[i] == 'h')
	    tplt_ss[i] = HELIX|COIL;
	else if (tplt_ss[i] == 'e')
	    tplt_ss[i] = STRAND|COIL;
	else
	    tplt_ss[i] = HELIX|STRAND|COIL;
    }
	
    fgets(buf, sizeof(buf), ifp);

    if (strlen(buf) != seqlen+1)
	fail("main: mismatching line length in aln file!");

    for (i=0; i<9; i++)
	ds_from[i] = -1;
    for (i=0; i<seqlen; i++)
	if (isdigit(buf[i]))
	{
	    j = buf[i] - '1';
	    buf[i] = 'C';
	    if (ds_from[j] < 0)
		ds_from[j] = i;
	    else
		ds_to[j] = i;
	    if (j+1 > n_ds)
		n_ds = j+1;
	}

    seq = (char **) malloc(nseqs * sizeof(char *));

    if (!seq)
	fail("main: out of memory!");
    seq[0] = (char*)malloc(seqlen);
    if (!seq[0])
	fail("main: out of memory!");

    memcpy(seq[0], buf, seqlen);

    const char* residues[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
    const char* residues_single[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};

    for(int d = 0; d < 20; d++){
        aacodes_num[residues[d]] = d;
        aacodes_single_num[residues_single[d]] = d;
        //cout << residues_single[d] << " ----> " << d << endl;
     }

    for (i = 0; i < seqlen; i++){
        res_type[i] = aanum(seq[0][i]);
    }    

    for (i = 1; i < nseqs; i++)
    {
	seq[i] = (char*)malloc(seqlen);
	if (!seq[i])
	    fail("main: out of memory!");
	if (!fgets(buf, MAXSEQLEN, ifp))
	    fail("main: error while reading aln file!");
	if (strlen(buf) != seqlen+1)
	    fail("main: mismatching line length in aln file!");
	memcpy(seq[i], buf, seqlen);
    }

    conseq = (char*)malloc(seqlen);
    if (!conseq)
	fail("main: out of memory!");

    for (i = 0; i < nseqs; i++)
	for (j = 0; j < seqlen; j++)
	{
	    if (!isalpha(seq[i][j]))
	    {
		seq[i][j] = GAP;
		gaps[j]++;
	    }
	    else
	    {
		aa = aanum(seq[i][j]);
		switch (aa)
		{
		case 20:
		    /* ASP is more common than ASN! */
		    seq[i][j] = ASP;
		    break;
		case 21:
		    /* GLU is more common than GLN! */
		    seq[i][j] = GLU;
		    break;
		case 22:
		    seq[i][j] = UNK;
		    break;
		default:
		    seq[i][j] = aa;
		    break;
		}
	    }
	}

    fclose(ifp);

    /* Calculate consensus sequence */
    for (i = 0; i < seqlen; i++)
    {
	int aacount[20], aacmax = 0, caa;
	
	memset(aacount, 0, sizeof(aacount));
	
	for (j = 0; j < nseqs; j++)
	{
	    if (seq[j][i] < 20)
	    {
		aacount[seq[j][i]]++;
		if (aacount[seq[j][i]] > aacmax)
		{
		    caa = seq[j][i];
		    aacmax = aacount[caa];
		}
	    }
	}

	conseq[i] = caa;
    }

    for (i=0; i<seqlen; i++)
	zcoords[i] = 9999.9;
    
    cout << "Reading z-file...\n";
    if (zfname[0])
    {
	float tempz;
	
	i = 0;
    
	
	ifp = fopen(zfname, "r");
	if (!ifp)
	    fail("main: unable to open z-file!");
    int t = 0;
	for (i=0; i<seqlen; i++)
	{
	    if (!fgets(buf, 256, ifp))
		fail("Cannot read from z-file!");
	    if (sscanf(buf, "%*s%f", &tempz) == 1){
		  zcoords[i] = tempz;
          if(tempz == 15 || tempz == -15){
            topology[t] = i;
            if(t == 0 && tempz == 15){
                // 0 == in, 1 == out
                n_term = 1;
            }
            t++;
          }
        }
	}
    helices = t/2;
    /*
    printf("TM helices:\t%i\n",helices);
    for(t=0;t<helices*2;t+=2){
        printf("%i - %i\n",topology[t],topology[t+1]);
    }
    exit(1);
    */
    }    
    cout << "done.\n";

    cout << "Reading lipex file...\n";
    if(lipexfname[0]){
        float tempz;

        ifp = fopen(lipexfname, "r");
        if (!ifp)
            fail("main: unable to open lipid exposure file!");

        for (i=0; i<seqlen; i++){
            if (!fgets(buf, 256, ifp))
                fail("Cannot read from lipid exposure file!");
            if (sscanf(buf, "%f", &tempz) == 1){
                hhlipex[i] = tempz;
            }
        }
    }
    cout << "done.\n";

    if (confname[0])
    {
	conmat = (float**)allocmat(seqlen, seqlen, sizeof(float));    
	dcomat = (float**)allocmat(seqlen, seqlen, sizeof(float));    

	rrcbflag = 0;
	ifp = fopen(confname, "r");
	if (!ifp)
	    fail("main: unable to open con file!");
	
	while (!feof(ifp))
	{
	    if (!fgets(buf, sizeof(buf), ifp))
		break;

	    if (!strncmp(buf, "PFRMAT RR", 9))
	    {
		rrcbflag = 1;
		continue;
	    }

	    if (!strncmp(buf, "RR CONSTR", 9))
	    {
		rrcbflag = 2;
		continue;
	    }

	    if (isalpha(buf[0]))
		continue;

	    if (rrcbflag == 2)
	    {
		if (sscanf(buf, "%d%d%f%f", &i, &j, &eqd, &dsd) != 4)
		    continue;
		conmat[i-1][j-1] = conmat[j-1][i-1] = eqd;
		dcomat[i-1][j-1] = dcomat[j-1][i-1] = dsd;
	    }
	    else
	    {
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
	}
	
	fclose(ifp);
    }

    if (zfname[0])
    {
	for (k=i=0; i<seqlen; i++)
	    if (zcoords[i] < 1000.0F)
		for (j=i+5; j<seqlen; j++)
		    if (zcoords[j] < 1000.F && dcomat[i][j] > 0.0F && fabsf(zcoords[i]-zcoords[j]) > 14.0F)
		    {
			dcomat[i][j] = conmat[i][j] = dcomat[j][i] = conmat[j][i] = 0.0F;
			k++;
		    }
	
	if (verboseflg)
	    printf("%d contacts excluded by Z-coordinate filter\n", k);
    }
    
    dsqmat = (float (**)[NPAIRS])allocmat(seqlen, seqlen, sizeof(**dsqmat));

    bestchn = (RESATM*)calloc(seqlen, sizeof(RESATM));

    if (bestchn == NULL)
	fail("main: calloc bestchn failed!");

    curchn = (RESATM*)calloc(seqlen, sizeof(RESATM));

    for (i=0; i<seqlen; i++)
	curchn[i].sc = NULL;

    if (curchn == NULL)
	fail("main: calloc curchn failed!");

    oldchn = (RESATM*)calloc(seqlen, sizeof(RESATM));

    if (oldchn == NULL)
	fail("main: calloc oldchn failed!");

    targchn = (RESATM*)calloc(seqlen, sizeof(RESATM));
    
    if (targchn == NULL)
	fail("main: calloc targchn failed!");

    if (verboseflg)
	puts("\nReading fold library...");
    readchn();

    if (argc > 2)
    {
	outpdbn = argv[2];
	if (fopen(argv[2], "r") != NULL)
        printf("main: output PDB file exists!\n");
	    //fail("main: output PDB file exists!");
    }

    readrots();

    if (MAXFRAGS > 0)
	mkfragtb();

    if (MAXFRAGS2 > 0)
	mkfragtb2();

    (void) cputime();

    repmc();

    calcdist(bestchn, FALSE);
    sc_opt(dsqmat, bestchn); /* Optimize side chains to reduce steric clashes */

    /* Print final energy report */
    verboseflg = TRUE;

    printf("Final energy report:\n");
    //MAXGACALLS = BIG;
    //GAPOOLSIZE = 10000;
    e = e_model(dsqmat, bestchn, TRUE);

    puts("FINISHED!");

    /* Orientate and write out final best structure */

    PDB* protein = new PDB;
    for (i = 0; i < MIN(curr_seqlen,seqlen); i++){
        if(res_type[i] != 7){
            Point3d p(res_type[i],i,bestchn[i].cb[0],bestchn[i].cb[1],bestchn[i].cb[2]);
            protein->add_cb(p);
        }else{
            Point3d p(res_type[i],i,bestchn[i].ca[0],bestchn[i].ca[1],bestchn[i].ca[2]);
            protein->add_cb(p);                
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
    ga->set_poolsize(10000);
    ga->set_maxcalls(1000000);
    ga->set_target(protein);
    vector<double> results = ga->run_ga();
    delete ga;              

    x_rot_best = results[0];
    y_rot_best = results[1];
    z_trans_best = results[2];
    e = protein->orientate(x_rot_best,y_rot_best,z_trans_best);
    delete protein;      

    writepdb(bestchn, e);

    return 0;
}
