/* Copyright (c) Signal Analysis and Imaging Group (SAIG), University of Alberta, 2013.*/
/* All rights reserved.                       */
/* sudecimate  :  $Date: May 2013- Last version May 2013  */

#include "su.h"
#include "cwp.h"
#include "segy.h"
#include "header.h"
#include <time.h>

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   							                                  ",
  " SUDECIMATE  decimate a 2d dataset randomly or regularly           ",
  "                                                                   ",
  " User provides:                                                    ",
  "           <in.su, >out.su                                         ",
  "                                                                   ",
  " Other parameters: (parameter and its default setting)             ",
  "           verbose=0; (=1 to show messages)                        ",
  "           method=0; (0=random decimation, 1=regular decimation)   ",
  "           dec=0.25; (random decimation of 25% of traces)          ",
  "           j=4;      (regular decimation of every 4th trace)       ",
  "                                                                   ",
  " Example:                                                          ",
  "  # remove 75% of traces randomly                                  ",
  "  sudecimate < in.su > out.su method=0 dec=0.75                    ",
  "  # remove every 4th trace                                         ",
  "  sudecimate < in.su > out.su method=1 j=4                         ",
  "                                                                   ",
  " Future development:                                               ",
  "  headers other than offset are not maintained                     ",
  "  choice of seed for random number generator should be added       ",
  "  add conditions of decimation to be met (max gap, etc)            ",
  "  types of decimation could be added (jitter sampling, max entropy,...)",
  "  higher dimensions (3d,5d)                                        ",
  "                                                                   ",
 NULL};
/* Credits:
 * Aaron Stanton
 * Trace header fields accessed: 
 * Last changes: May : 2013 
 */
/**************** end self doc ***********************************/

segy tr;
int main(int argc, char **argv)
{
  int verbose;
  time_t start,finish;
  double elapsed_time;
  int it,ix,ix_dec,nt,nx,nx_dec,method,j,counter;
  float dt;
  float *h,*h_dec;
  float **din, **dout;
  float rndnum, dec;

  /********/    
  fprintf(stderr,"*******SUDECIMATE*********\n");
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  start=time(0);    
  /* Get parameters */
  if (!getparint("method", &method)) method = 0; /* 0=random decimation, 1=regular decimation */
  if (!getparint("verbose", &verbose)) verbose = 0;
  if (!getparint("nx", &nx)) nx = 10000;
  if (!getparfloat("dec", &dec)) dec = 0.55; /* remove 55% of input traces */
  if (!getparint("j", &j)) j = 2; /* remove every 2nd trace */
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;

  din   = ealloc2float(nt,nx);
  dout  = ealloc2float(nt,nx);
  h     = ealloc1float(nx);
  h_dec = ealloc1float(nx);
  /* ***********************************************************************
  input data
  *********************************************************************** */
  ix=0;
  do {
    h[ix]=(float)  tr.offset;
    memcpy((void *) din[ix],(const void *) tr.data,nt*sizeof(float));
    ix++;
    if (ix > nx) err("Number of traces > %d\n",nx); 
  } while (gettr(&tr));
  erewind(stdin);
  nx=ix;
  if (verbose) fprintf(stderr,"processing %d traces \n", nx);
  ix_dec = 0;
  
  if (method==0){
    for (ix=0;ix<nx;ix++){    
      rndnum = franuni();/* generate random number distributed between 0 and 1 */
      if (rndnum>dec){
      h_dec[ix_dec] = h[ix];
      for (it=0;it<nt;it++) dout[ix_dec][it] = din[ix][it];
      ix_dec++;
      }
    }
    nx_dec = ix_dec;
  }
  else {
    counter = 1;
    for (ix=0;ix<nx;ix++){
      if (counter!=j){
      h_dec[ix_dec] = h[ix];
      for (it=0;it<nt;it++) dout[ix_dec][it] = din[ix][it];
      ix_dec++;
      }
      else counter = 0;
      counter++;
    }
    nx_dec = ix_dec;
  }
  if (verbose) fprintf(stderr,"removed %d traces \n", nx - nx_dec);

  /* ***********************************************************************
  output data
  *********************************************************************** */
  rewind(stdin);
  for (ix=0;ix<nx_dec;ix++){ 
    memcpy((void *) tr.data,(const void *) dout[ix],nt*sizeof(float));
    tr.offset=(int) h_dec[ix];
    tr.ntr=nx_dec;
    tr.ns=nt;
    tr.dt = NINT(dt*1000000.);
    fputtr(stdout,&tr);
  }
  /******** End of output **********/
 return EXIT_SUCCESS;
}

