#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "longnam.h"
#include "arrays.h"

long Globalidum;

#define BLOCKSIZE 1e7

typedef struct {
  int X1,Y1;
  int X2,Y2;
  int Width;
} SatTrackRec;

typedef struct {
  int Num;
  SatTrackRec *PSatTrackProp;
} SatTrackSetRec;

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float quick_select(float arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]);
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]);
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]);
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]);

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]);

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]);
        do hh--; while (arr[hh]  > arr[low]);

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]);
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]);

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

extern void printerror( int status)
{

  if (status) {
    fits_report_error(stderr, status); 
    
    exit( status ); 
  }
  return;
}

extern void SReadFits(char *FN, float ***ptfs, int *dimx, int *dimy,
		      int WhImage)
{
  int HDU;
  /*  char *comv[MAX_COMMENTS];*/
  int status = 0;
  long firstpix = 1, npixels;
  long anaxes[] = {0,0};
  float *apix;
  fitsfile *afptr;
  int hdutype,TotPass,YBlock,YLine,Off;
  int I,J,K;

  fits_open_file(&afptr, FN, READONLY, &status);

  HDU = WhImage+1;

  if (fits_movabs_hdu(afptr,HDU,&hdutype,&status))
    printerror( status );

  fits_get_img_size(afptr, 2, anaxes, &status);
  *dimx = anaxes[0];
  *dimy = anaxes[1];
  allocFloatArray(ptfs,*dimx,*dimy);
  if (anaxes[0] > 0) {
    YBlock = (int)(BLOCKSIZE / anaxes[0]);
    
    TotPass = (anaxes[1] / YBlock) + 1;
    for (K=0;K<TotPass;K++) {
      if (K < TotPass-1)
	YLine = YBlock;
      else
	YLine = *dimy - YBlock*(TotPass-1);
      
      npixels = YLine * anaxes[0];
      firstpix = anaxes[0] * YBlock * K + 1;
      apix = (float *)calloc(npixels,sizeof(float));
      fits_read_img(afptr, TFLOAT, firstpix, npixels,
		    NULL, apix, NULL, &status);
      Off = YBlock * K;
      for (I=0;I<*dimx;I++)
	for (J=0;J<YLine;J++)
	  (*ptfs)[J+Off][I] = apix[I+J*(*dimx)];
      free(apix);
    }
  }
  fits_close_file(afptr, &status);
}

extern void WriteFloatFits(char *FN, float **tfs, int dimx, int dimy)
{
  int status,naxis,J;
  float *ffs;
  fitsfile *fptr;
  long fpixel,naxes[2];
  char OFN[200];

  status = 0;
  sprintf(OFN,"!%s",FN);
  if (fits_create_file(&fptr,OFN,&status))
    printerror( status );

  naxis = 2;
  fpixel = 1;
  naxes[0] = dimx;
  naxes[1] = dimy;

  if (fits_create_img(fptr,FLOAT_IMG,naxis,naxes,&status))
    printerror( status );

  ffs = (float *)calloc(dimx*dimy,sizeof(float));
  for (J=0;J<dimx*dimy;J++) 
    ffs[J] = tfs[J/dimx][J%dimx];
  if (fits_write_img(fptr,TFLOAT,fpixel,dimx*dimy,ffs,&status))
    printerror( status );
  free(ffs);
  if (fits_close_file(fptr,&status))
    printerror( status );
}

extern int MakeInteger(float F)
{
  int V;
  if (F >= 0) V = (int)F;
  else V = (int)F - 1;
  return V;
}

extern int LineSearch(int X1, int Y1, int X2, int Y2, float **nfs, 
		      float lowlim)
{
  int Found;
  
  if (fabs(Y1-Y2)>fabs(X1-X2)) {
    int Y,NGood,X;

    NGood = 0;
    if (Y2 > Y1) {
      for (Y=Y1;Y<Y2;Y++) {
	X = MakeInteger(X1 + (float)(X2-X1)*(Y-Y1)/(Y2-Y1)+0.5);
	if (nfs[Y][X] > lowlim) 
	  NGood++;
      }
    }
    else {
      for (Y=Y2;Y<Y1;Y++) {
	X = MakeInteger(X1 + (float)(X2-X1)*(Y-Y1)/(Y2-Y1)+0.5);
	if (nfs[Y][X] > lowlim) 
	  NGood++;
      }
    }
    if (NGood > 0.75*fabs(Y2-Y1))
      Found = 1;
    else
      Found = 0;
  }
  else {
    int X,NGood,Y;

    NGood = 0;
    if (X2 > X1) {
      for (X=X1;X<X2;X++) {
	Y = MakeInteger(Y1 + (float)(Y2-Y1)*(X-X1)/(X2-X1)+0.5);
	if (nfs[Y][X] > lowlim) 
	  NGood++;
      }
    }
    else {
      for (X=X2;X<X1;X++) {
	Y = MakeInteger(Y1 + (float)(Y2-Y1)*(X-X1)/(X2-X1)+0.5);
	if (nfs[Y][X] > lowlim) 
	  NGood++;
      }      
    }
    if (NGood > 0.75*fabs(X2-X1))
      Found = 1;
    else
      Found = 0;
  }
  return Found;
}

extern int WithinSatTrackProp(int X, int Y, SatTrackRec *PSatTrackProp)
{
  int OK;
  float PredX,PredY;

  if (fabs(PSatTrackProp->X2 - PSatTrackProp->X1) <
      fabs(PSatTrackProp->Y2 - PSatTrackProp->Y1)) {
    PredX = PSatTrackProp->X1 +
      (PSatTrackProp->X2-PSatTrackProp->X1)*
      (Y-PSatTrackProp->Y1)/
      (PSatTrackProp->Y2-PSatTrackProp->Y1);
    if (fabs(X - PredX) < PSatTrackProp->Width)
      OK = 1;
    else
      OK = 0;
  }
  else {
    PredY = PSatTrackProp->Y1 +
      (PSatTrackProp->Y2-PSatTrackProp->Y1)*
      (X-PSatTrackProp->X1)/
      (PSatTrackProp->X2-PSatTrackProp->X1);
    if (fabs(Y - PredY) < PSatTrackProp->Width)
      OK = 1;
    else
      OK = 0;    
  }
  return OK;
}

extern int LineAlreadyExist(int X1, int Y1, int X2, int Y2, 
			    SatTrackSetRec *PSatTrackSetProp)
{
  int Found,I;

  Found = 0;
  for (I=0;I<PSatTrackSetProp->Num;I++) {
    SatTrackRec *PSatTrackProp;

    PSatTrackProp = &PSatTrackSetProp->PSatTrackProp[I];
    if (WithinSatTrackProp(X1,Y1,PSatTrackProp)&&
	WithinSatTrackProp(X2,Y2,PSatTrackProp))
      Found = 1;
    if (!Found) {
      if (WithinSatTrackProp(X1,Y1,PSatTrackProp)+
	  WithinSatTrackProp(X2,Y2,PSatTrackProp)+
	  WithinSatTrackProp((X1+X2)/2,(Y1+Y2)/2,PSatTrackProp)>1)
	Found = 1;
    }
  }
  return Found;
}

extern float CalcScoreDX(int X1, int Y1, int X2, int Y2,
			 float **nfs, float lowlim, int ndimx, int ndimy,
			 float *PIntScore)
{
  int X,Y,NGood,Tot;
  float *POrd,Score;

  POrd = (float *)calloc(fabs(X2-X1),sizeof(float));
  NGood = 0;
  Tot = 0;
  if (X2 > X1) {
    for (X=X1;X<X2;X++) {   
      Y = MakeInteger(Y1 + (float)(X-X1)*(Y2-Y1)/(X2-X1)+0.5);
      
      if ((Y>=0)&&(Y<ndimy)) {
	if (nfs[Y][X] > lowlim)
	  NGood++;
	POrd[Tot] = nfs[Y][X];
	Tot++;
      }
    }
  }
  else {
    for (X=X2;X<X1;X++) {   
      Y = MakeInteger(Y1 + (float)(X-X1)*(Y2-Y1)/(X2-X1)+0.5);
      
      if ((Y>=0)&&(Y<ndimy)) {
	if (nfs[Y][X] > lowlim)
	  NGood++;
	POrd[Tot] = nfs[Y][X];
	Tot++;
      }
    }
  }
  *PIntScore = quick_select(POrd,Tot);
  if (NGood > Tot*0.6)
    Score = *PIntScore;
  else
    Score = 0;    
  free(POrd);

  return Score;
}

extern float CalcScoreDY(int X1, int Y1, int X2, int Y2,
			 float **nfs, float lowlim, int ndimx, int ndimy,
			 float *PIntScore)
{
  int X,Y,NGood,Tot;
  float *POrd,Score;  

  POrd = (float *)calloc(fabs(Y2-Y1),sizeof(float));
  NGood = 0;
  Tot = 0;
  if (Y2 > Y1) {
    for (Y=Y1;Y<Y2;Y++) {   
      X = MakeInteger(X1 + (float)(Y-Y1)*(X2-X1)/(Y2-Y1)+0.5);
      
      if ((X>=0)&&(X<ndimx)) {
	if (nfs[Y][X] > lowlim)
	  NGood++;
	POrd[Tot] = nfs[Y][X];
	Tot++;
      }
    }
  }
  else {
    for (Y=Y2;Y<Y1;Y++) {   
      X = MakeInteger(X1 + (float)(Y-Y1)*(X2-X1)/(Y2-Y1)+0.5);
      
      if ((X>=0)&&(X<ndimx)) {
	if (nfs[Y][X] > lowlim)
	  NGood++;
	POrd[Tot] = nfs[Y][X];
	Tot++;
      }
    }    
  }
  *PIntScore = quick_select(POrd,Tot);
  if (NGood > Tot*0.6)
    Score = *PIntScore;
  else
    Score = 0;

  free(POrd);

  return Score;
}

extern float CalcScore(int X1, int Y1, int X2, int Y2,
		       float **nfs, float lowlim, int ndimx, int ndimy,
		       float *PIntScore)
{
  if (fabs(X2-X1)<fabs(Y2-Y1))
    return CalcScoreDY(X1,Y1,X2,Y2,nfs,lowlim,ndimx,ndimy,PIntScore);
  else
    return CalcScoreDX(X1,Y1,X2,Y2,nfs,lowlim,ndimx,ndimy,PIntScore);
}

extern void TryExtensionDX(int X1, int Y1, int X2, int Y2,
			   float **nfs, int ndimx, int ndimy, float lowlim, 
			   int Tweak, int Direction, int Size,
			   float *PScore, int *PX, int *PY)
{
  int NX,NY;
  float TScore;

  if (Direction==0) {
    NX = X1 - 25;
    if (NX < 0)
      NX = 0;
  }
  else {
    NX = X2 + 25;
    if (NX >= ndimx)
      NX = ndimx-1;
  }

  NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)/(X2-X1)+0.5)+Tweak;
  if (NY < 0) {
    NY = 0;
    if (fabs(Y2-Y1)>0)
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)/(Y2-Y1)+0.5)+Tweak;
    else if (Direction == 0)
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)+0.5)+Tweak;
    else
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)/(-1.0)+0.5)+Tweak;
  }
  else if (NY >= ndimy) {
    NY = ndimy-1;
    if (fabs(Y2-Y1)>0)
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)/(Y2-Y1)+0.5)+Tweak;
    else if (Direction == 0)
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)/(-1.0)+0.5)+Tweak;
    else
      NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)+0.5)+Tweak;
  }
  if (Direction == 0)
    *PScore = CalcScore(NX,NY,X1,Y1,nfs,lowlim,ndimx,ndimy,&TScore);
  else 
    *PScore = CalcScore(X2,Y2,NX,NY,nfs,lowlim,ndimx,ndimy,&TScore);  
  *PX = NX;
  *PY = NY;
}

extern void TryExtensionDY(int X1, int Y1, int X2, int Y2,
			   float **nfs, int ndimx, int ndimy, float lowlim, 
			   int Tweak, int Direction, int Size,
			   float *PScore, int *PX, int *PY)
{
  int NX,NY;
  float TScore;

  if (Direction==0) {
    NY = Y1 - 25;
    if (NY < 0)
      NY = 0;
  }
  else {
    NY = Y2 + 25;
    if (NY >= ndimy)
      NY = ndimy-1;
  }

  NX = MakeInteger(X1 + (float)(NY-Y1)*(X2-X1)/(Y2-Y1)+0.5)+Tweak;
  if (NX < 0) {
    NX = 0;
    if (fabs(X1-X2)>0)
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)/(X2-X1)+0.5)+Tweak;
    else if (Direction == 0)
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)+0.5)+Tweak;
    else
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)/(-1.0)+0.5)+Tweak;
  }
  else if (NX >= ndimx) {
    NX = ndimx-1;
    if (fabs(Y1-Y2)>0)
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)/(X2-X1)+0.5)+Tweak;
    else if (Direction == 0)
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)/(-1.0)+0.5)+Tweak;
    else
      NY = MakeInteger(Y1 + (float)(NX-X1)*(Y2-Y1)+0.5)+Tweak;
  }
  if (Direction == 0)
    *PScore = CalcScore(NX,NY,X1,Y1,nfs,lowlim,ndimx,ndimy,&TScore);
  else 
    *PScore = CalcScore(X2,Y2,NX,NY,nfs,lowlim,ndimx,ndimy,&TScore);
  *PX = NX;
  *PY = NY;
}

extern void TryExtension(int X1, int Y1, int X2, int Y2, int BX,
			 float **nfs, int ndimx, int ndimy, float lowlim, 
			 int Tweak, int Direction, int Size,
			 float *PScore, int *PX, int *PY)
{
  if (!BX)
    TryExtensionDY(X1,Y1,X2,Y2,nfs,ndimx,ndimy,lowlim,Tweak,Direction,Size,
		   PScore,PX,PY);
  else
    TryExtensionDX(X1,Y1,X2,Y2,nfs,ndimx,ndimy,lowlim,Tweak,Direction,Size,
		   PScore,PX,PY);
}

extern void WidthLXY(int X2, int Y2, int X1, int Y1, int ndimx, int ndimy,
		     int *PNX2, int *PNY2, int Width)
{
  if (X2-Width < 0) {
    *PNX2 = 0;
    *PNY2 = MakeInteger(Y2 + (float)(Y2-Y1)*(Width-X2)/(X2-X1)+0.5);
  }
  else {
    if (X2 < ndimx) {
      *PNX2 = X2-Width;
      *PNY2 = Y2;
    }
    else {
      *PNX2 = ndimx;
      *PNY2 = MakeInteger(Y2 + (float)(Y2-Y1)*Width/(X2-X1)+0.5);
    }
  }  
}

extern void WidthHXY(int X2, int Y2, int X1, int Y1, int ndimx, int ndimy,
		     int *PNX2, int *PNY2, int Width)
{
  if (X2+Width >= ndimx) {
    *PNX2 = ndimx-1;
    *PNY2 = MakeInteger(Y2 - (float)(Y2-Y1)*(X2+Width-ndimx+1)/(X2-X1)+0.5);
  }
  else {
    if (X2 > 0) {
      *PNX2 = X2+Width;
      *PNY2 = Y2;
    }
    else {
      *PNX2 = 0;
      *PNY2 = MakeInteger(Y2 - (float)(Y2-Y1)*Width/(X2-X1)+0.5);      
    }
  }  
}

extern int CalculateLineWidthVarX(float **nfs, int ndimx, int ndimy, 
				  float lowlim, 
				  int X1, int Y1, int X2, int Y2)
{
  float Score;
  int WidthL,WidthH,Done;

  WidthL = 0;
  Done = 0;
  while (!Done) {
    int NX1,NY1,NX2,NY2;

    WidthLXY(X1,Y1,X2,Y2,ndimx,ndimy,&NX1,&NY1,WidthL);
    WidthLXY(X2,Y2,X1,Y1,ndimx,ndimy,&NX2,&NY2,WidthL);

    CalcScore(NX1,NY1,NX2,NY2,nfs,lowlim,ndimx,ndimy,&Score);
    if (Score < lowlim)
      Done = 1;
    else
      WidthL++;
  }

  WidthH = 0;
  Done = 0;
  while (!Done) {
    int NX1,NY1,NX2,NY2;

    WidthHXY(X1,Y1,X2,Y2,ndimx,ndimy,&NX1,&NY1,WidthH);
    WidthHXY(X2,Y2,X1,Y1,ndimx,ndimy,&NX2,&NY2,WidthH);

    CalcScore(NX1,NY1,NX2,NY2,nfs,lowlim,ndimx,ndimy,&Score);
    if (Score < lowlim)
      Done = 1;
    else
      WidthH++;
  }
  return MakeInteger((WidthL+WidthH)/2.0+0.5);
}

extern int CalculateLineWidthVarY(float **nfs, int ndimx, int ndimy, 
				  float lowlim, 
				  int X1, int Y1, int X2, int Y2)
{
  float Score;
  int WidthL,WidthH,Done;

  WidthL = 0;
  Done = 0;
  while (!Done) {
    int NX1,NY1,NX2,NY2;

    WidthLXY(Y1,X1,Y2,X2,ndimy,ndimx,&NY1,&NX1,WidthL);
    WidthLXY(Y2,X2,Y1,X1,ndimy,ndimx,&NY2,&NX2,WidthL);

    CalcScore(NX1,NY1,NX2,NY2,nfs,lowlim,ndimx,ndimy,&Score);
    if (Score < lowlim)
      Done = 1;
    else
      WidthL++;
  }

  WidthH = 0;
  Done = 0;
  while (!Done) {
    int NX1,NY1,NX2,NY2;

    WidthHXY(Y1,X1,Y2,X2,ndimy,ndimx,&NY1,&NX1,WidthH);
    WidthHXY(Y2,X2,Y1,X1,ndimy,ndimx,&NY2,&NX2,WidthH);

    CalcScore(NX1,NY1,NX2,NY2,nfs,lowlim,ndimx,ndimy,&Score);
    if (Score < lowlim)
      Done = 1;
    else
      WidthH++;
  }
  return (WidthL+WidthH)/2.0;
}

extern float CalculateLineWidth(float **nfs, int ndimx, int ndimy, 
				float lowlim, 
				int X1, int Y1, int X2, int Y2)
{
  int Width;

  if (fabs(X2-X1)<fabs(Y2-Y1))
    Width = CalculateLineWidthVarX(nfs,ndimx,ndimy,lowlim,X1,Y1,X2,Y2);
  else
    Width = CalculateLineWidthVarY(nfs,ndimx,ndimy,lowlim,X1,Y1,X2,Y2);

  Width = (int)(Width + 1);
  return Width;
}

extern void AddSatTrackProp(int X1, int Y1, int X2, int Y2,
			    float Width, SatTrackSetRec *PSatTrackSetProp)
{
  SatTrackRec *PSatTrackProp;

  if (PSatTrackSetProp->Num) {
    PSatTrackSetProp->PSatTrackProp = (SatTrackRec *)
      realloc(PSatTrackSetProp->PSatTrackProp,(PSatTrackSetProp->Num+1)*
	      sizeof(SatTrackRec));
  }
  else {
    PSatTrackSetProp->PSatTrackProp = (SatTrackRec *)
      calloc(1,sizeof(SatTrackRec));    
  }
  PSatTrackProp = &PSatTrackSetProp->
    PSatTrackProp[PSatTrackSetProp->Num];
  PSatTrackProp->X1 = X1;
  PSatTrackProp->Y1 = Y1;
  PSatTrackProp->X2 = X2;
  PSatTrackProp->Y2 = Y2;
  PSatTrackProp->Width = Width;
  (PSatTrackSetProp->Num)++;
}

extern void OptimizeLineVarY1(float **nfs, int ndimx, int ndimy, 
			      float lowlim,
			      int *PX1, int *PX2, int *PY1, int *PY2)
{
  int Direction;
  float OldScore,Score,TScore;

  OldScore = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
		       &TScore);

  for (Direction=0;Direction<2;Direction++) {
    int Done,Life,Cur;

    Done = 0;
    Life = 7;
    Cur = *PY1;

    while (!Done) {
      if (((Direction==0)&&(*PY1==0))||
	  ((Direction==1)&&(*PY1==ndimy-1)))
	Done = 1;
      else {
	if (Direction==0) (*PY1)--;
	else if (Direction==1) (*PY1)++;
	Score = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
			  &TScore);
	if (Score < OldScore*1.0001+0.1) {
	  Life--;
	  if (Life==0) {
	    Done = 1;
	  }
	}
	else {
	  Life = 7;
	  OldScore = Score;
	  Cur = *PY1;
	}
      }
    }
    (*PY1) = Cur;
  }
}

extern void OptimizeLineVarY2(float **nfs, int ndimx, int ndimy, 
			      float lowlim,
			      int *PX1, int *PX2, int *PY1, int *PY2)
{
  int Direction;
  float OldScore,Score,TScore;

  OldScore = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
		       &TScore);

  for (Direction=0;Direction<2;Direction++) {
    int Done,Life,Cur;

    Done = 0;
    Life = 7;
    Cur = *PY2;

    while (!Done) {
      if (((Direction==0)&&(*PY2==0))||
	  ((Direction==1)&&(*PY2==ndimy-1)))
	Done = 1;
      else {
	if (Direction==0) (*PY2)--;
	else if (Direction==1) (*PY2)++;
	Score = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
			  &TScore);
	if (Score < OldScore*1.0001+0.1) {
	  Life--;
	  if (Life==0) {
	    Done = 1;
	  }
	}
	else {
	  Life = 7;
	  OldScore = Score;
	  Cur = *PY2;
	}
      }
    }
    (*PY2) = Cur;
  }
}

extern void OptimizeLineVarX1(float **nfs, int ndimx, int ndimy, 
			      float lowlim,
			      int *PX1, int *PX2, int *PY1, int *PY2)
{
  int Direction;
  float OldScore,Score,TScore;

  OldScore = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
		       &TScore);

  for (Direction=0;Direction<2;Direction++) {
    int Done,Life,Cur;

    Done = 0;
    Life = 7;
    Cur = *PX1;

    while (!Done) {
      if (((Direction==0)&&(*PX1==0))||
	  ((Direction==1)&&(*PX1==ndimx-1))) 
	Done = 1;
      else {
	if (Direction==0) (*PX1)--;
	else if (Direction==1) (*PX1)++;
	Score = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
			  &TScore);
	if (Score < OldScore*1.0001+0.1) {
	  Life--;
	  if (Life==0) {
	    Done = 1;
	  }
	}
	else {
	  Life = 7;
	  OldScore = Score;
	  Cur = *PX1;
	}
      }
    }
    (*PX1) = Cur;
  }
}

extern void OptimizeLineVarX2(float **nfs, int ndimx, int ndimy, 
			      float lowlim,
			      int *PX1, int *PX2, int *PY1, int *PY2)
{
  int Direction;
  float OldScore,Score,TScore;

  OldScore = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
		       &TScore);

  for (Direction=0;Direction<2;Direction++) {
    int Done,Life,Cur;

    Done = 0;
    Life = 7;
    Cur = *PX2;

    while (!Done) {
      if (((Direction==0)&&(*PX2==0))||
	  ((Direction==1)&&(*PX2==ndimx-1))) 
	Done = 1;
      else {
	if (Direction==0) (*PX2)--;
	else if (Direction==1) (*PX2)++;
	Score = CalcScore(*PX1,*PY1,*PX2,*PY2,nfs,lowlim,ndimx,ndimy,
			  &TScore);
	if (Score < OldScore*1.0001+0.1) {
	  Life--;
	  if (Life==0) {
	    Done = 1;
	  }
	}
	else {
	  Life = 7;
	  OldScore = Score;
	  Cur = *PX2;
	}
      }
    }
    (*PX2) = Cur;
  }
}

extern void OptimizeLine(float **nfs, int ndimx, int ndimy, 
			 float lowlim,
			 int *PX1, int *PX2, int *PY1, int *PY2)
{
  if (fabs(*PX2-*PX1)<fabs(*PY2-*PY1)) {
    if ((*PX1 == 0)||(*PX1 == ndimx-1))
      OptimizeLineVarY1(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
    else
      OptimizeLineVarX1(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);

    if ((*PX2 == 0)||(*PX2 == ndimx-1))
      OptimizeLineVarY2(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
    else
      OptimizeLineVarX2(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
  }
  else {
    if ((*PY1 == 0)||(*PY1 == ndimy-1))
      OptimizeLineVarX1(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
    else
      OptimizeLineVarY1(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);

    if ((*PY2 == 0)||(*PY2 == ndimy-1))
      OptimizeLineVarX2(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
    else
      OptimizeLineVarY2(nfs,ndimx,ndimy,lowlim,PX1,PX2,PY1,PY2);
  }
}

extern int TouchingEdges(int ndimx, int ndimy, int X1, int X2, int Y1,
			 int Y2)
{
  int Touch;

  if ((X1 < 1) || (X1 > ndimx-2) || (Y1 < 1) || (Y1 > ndimy - 2))
    Touch = 1;
  else
    Touch = 0;
  if ((X2 > 0) && (X2 < ndimx-1) && (Y2 > 0) && (Y2 < ndimy - 1))
    Touch = 0;
  if ((X1 < 1) && (X2 < 1))
    Touch = 0;
  if ((X1 > ndimx-2) && (X2 > ndimx-2))
    Touch = 0;
  if ((Y1 < 1) && (Y2 < 1))
    Touch = 0;
  if ((Y1 > ndimy-2) && (Y2 > ndimy-2))
    Touch = 0;

  return Touch;
}

extern void ExploreLine(int X1, int Y1, int X2, int Y2,
			float **nfs, float lowlim, int ndimx, int ndimy,
			SatTrackSetRec *PSatTrackSetProp)
{
  int Direction,BX,IX1,IY1,IX2,IY2;

  IX1 = X1;
  IX2 = X2;
  IY1 = Y1;
  IY2 = Y2;
  if (fabs(X2-X1)>=fabs(Y2-Y1))
    BX = 1;
  else
    BX = 0;

  for (Direction=0;Direction<2;Direction++) {
    int Done,Length;

    Done = 0;
    Length = 0;
    while (!Done) {
      float BestScore;
      int Tweak,BestTweak,BestX,BestY;
      
      BestScore = 0;
      for (Tweak=-2;Tweak<3;Tweak++) {
	float Score;
	int NX,NY,Size;
	
	if (Length == 0) Size = 25;
	else if (Length == 1) Size = 18;
	else if (Length == 2) Size = 10;

	TryExtension(X1,Y1,X2,Y2,BX,nfs,ndimx,ndimy,lowlim,Tweak,Direction,
		     Size,&Score,&NX,&NY);
	if (Score > BestScore) {
	  BestScore = Score;
	  BestTweak = Tweak;
	  BestX = NX;
	  BestY = NY;	  
	}
      }
      if (BestScore > 0) {
	if (Direction == 0) {
	  X1 = BestX;
	  Y1 = BestY;
	  if (BX) {
	    if ((X1 == 0)||(Y1==0)||(Y1==ndimy-1))
	      Done = 1;
	  }
	  else {
	    if ((Y1 == 0)||(X1==0)||(X1==ndimx-1))
	      Done = 1;
	  }
	}
	else {
	  X2 = BestX;
	  Y2 = BestY;
	  if (BX) {
	    if ((X2 == ndimx-1)||(Y2==0)||(Y2==ndimy-1))
	      Done = 1;
	  }
	  else {
	    if ((Y2 == ndimy-1)||(X2==0)||(X2==ndimx-1))
	      Done = 1;
	  }
	}
      }
      else {
	Length++;
	if (Length == 3)
	  Done = 1;
      }
    }
  }
  if (1) {
    int Edge,Leng;

    Edge = TouchingEdges(ndimx,ndimy,X1,X2,Y1,Y2);
    Leng = sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1));
    if (Leng > 80) {
      float Width;
      int OK;
      
      OK = 0;
      if (Edge == 1) OK = 1;
      OptimizeLine(nfs,ndimx,ndimy,lowlim,&X1,&X2,&Y1,&Y2);
      Width = CalculateLineWidth(nfs,ndimx,ndimy,lowlim,X1,Y1,X2,Y2);
      Leng = sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1));
      /*      if ((Leng > 105)&&(Width > 17)) OK = 1;*/
      if (OK) {
	if (!LineAlreadyExist(X1,Y1,X2,Y2,PSatTrackSetProp)) {
	  fprintf(stderr,"%i %i %i %i\n",X1,Y1,X2,Y2);
	  AddSatTrackProp(X1,Y1,X2,Y2,Width,PSatTrackSetProp);
	}
      }
    }
  }
}

float ran3(long *idum);

extern void IterateOverLineSearches(int *X1, int *Y1, int *X2, int *Y2, 
				    int Num2, 
				    float **nfs, int ndimx, int ndimy,
				    float lowlim, SatTrackSetRec *PSatTrackSetProp)
{
  int NPoint1,NPoint2,*ind1,*ind2,I;
  int TotTrial,Ind1[25],Ind2[75];

  NPoint1 = 0;
  for (I=0;I<25;I++)
    if (nfs[Y1[I]][X1[I]] > lowlim) {
      Ind1[NPoint1] = I;
      NPoint1++;
    }
  NPoint2 = 0;
  for (I=0;I<Num2;I++)
    if (nfs[Y2[I]][X2[I]] > lowlim) {
      Ind2[NPoint2] = I;
      NPoint2++;
    }

  TotTrial = NPoint1*NPoint2;
  if (TotTrial > 400) {
    TotTrial = 400;
   
    ind1 = (int *)calloc(TotTrial,sizeof(int));
    ind2 = (int *)calloc(TotTrial,sizeof(int));
    
    for (I=0;I<TotTrial;I++) {
      ind1[I] = (int)(ran3(&Globalidum)*NPoint1);
      ind2[I] = (int)(ran3(&Globalidum)*NPoint2);
    }
  }
  else {
    ind1 = (int *)calloc(TotTrial,sizeof(int));
    ind2 = (int *)calloc(TotTrial,sizeof(int));
    
    for (I=0;I<TotTrial;I++) {
      ind1[I] = (I % NPoint1);
      ind2[I] = (I / NPoint1);
    }
  }
  for (I=0;I<TotTrial;I++) {
    int CX1,CX2,CY1,CY2,t;
    
    CX1 = X1[Ind1[ind1[I]]];
    CY1 = Y1[Ind1[ind1[I]]];
    CX2 = X2[Ind2[ind2[I]]];
    CY2 = Y2[Ind2[ind2[I]]];
    if (fabs(CX1-CX2)<fabs(CY1-CY2)) {
      if (CY1 > CY2) {
	t = CY1; CY1 = CY2; CY2 = t;
      }
    }
    else {
      if (CX1 > CX2) {
	t = CX1; CX1 = CX2; CX2 = t;
      }
    }
    if (LineSearch(CX1,CY1,CX2,CY2,nfs,lowlim)) {      
      OptimizeLine(nfs,ndimx,ndimy,lowlim,&CX1,&CX2,&CY1,&CY2);
      if (!LineAlreadyExist(CX1,CY1,CX2,CY2,PSatTrackSetProp)) {
	ExploreLine(CX1,CY1,CX2,CY2,nfs,lowlim,ndimx,ndimy,PSatTrackSetProp);
      }
    }
  }
  free(ind1);
  free(ind2);
}

extern void IdentifyLines(float **fs, int dimx, int dimy, SatTrackSetRec *PSatTrackSetProp)
{
  float **nfs,Tot,*f,med,sig,lowlim;
  int *ind;
  int ndimx,ndimy,cdimx,cdimy,Num2,Low,High;
  int I,J,K,upv,I0,J0;
  int *X1,*X2,*Y1,*Y2;

  ndimx = dimx / 4;
  ndimy = dimy / 4;
  allocFloatArray(&nfs,ndimx,ndimy);
  for (I=0;I<ndimx;I++)
    for (J=0;J<ndimy;J++) {
      Tot = 0.0;
      for (I0=0;I0<4;I0++)
	for (J0=0;J0<4;J0++)
	  Tot += fs[J*4+J0][I*4+I0];
      nfs[J][I] = Tot;
    }
  WriteFloatFits("A",nfs,ndimx,ndimy);
  f = (float *)calloc(ndimx*ndimy,sizeof(float));
  for (I=0;I<ndimx;I++)
    for (J=0;J<ndimy;J++) {
      f[J*ndimx+I] = nfs[J][I];
    }
  ind = (int *)calloc(ndimx*ndimy,sizeof(int));
  hpsort(ndimx*ndimy,f,ind);
  med = f[ind[ndimx*ndimy/2]];
  sig = sqrt(f[ind[ndimx*ndimy/2]]);
  upv = ndimx*ndimy*0.8;
  lowlim = f[ind[upv]];
  lowlim = med+sig*4;
  fprintf(stderr,"Median = %g, Sigma = %g\n",med,sig);
  fprintf(stderr,"Lower Limit = %g\n",lowlim);

  cdimx = ndimx/25;
  cdimy = ndimy/25;
  X1 = (int *)calloc(25,sizeof(int));
  Y1 = (int *)calloc(25,sizeof(int));
  X2 = (int *)calloc(75,sizeof(int));
  Y2 = (int *)calloc(75,sizeof(int));
  PSatTrackSetProp->Num = 0;
  for (I=0;I<cdimx;I++)
    for (J=0;J<cdimy;J++) {            
      if (J*25+50<ndimy) {
	for (K=0;K<25;K++) {
	  X1[K] = I*25+K;
	  Y1[K] = J*25+24;
	}
	if (I>0) Low = I*25-25;
	else Low = I*25;
	if (I<cdimx-1) High = I*25+50;
	else High = I*25+25;
	Num2 = High-Low;
	for (K=0;K<Num2;K++) {
	  X2[K] = Low+K;
	  Y2[K] = J*25+50;
	}	
	IterateOverLineSearches(X1,Y1,X2,Y2,Num2,nfs,ndimx,ndimy,
				lowlim,PSatTrackSetProp);
      }
      if (I*25+50<ndimx) {
	for (K=0;K<25;K++) {
	  X1[K] = I*25+24;
	  Y1[K] = J*25+K;
	}
	if (J>0) Low = J*25-25;
	else Low = J*25;
	if (J<cdimy-1) High = J*25+50;
	else High = J*25+25;
	Num2 = High-Low;
	for (K=0;K<Num2;K++) {
	  X2[K] = I*25+50;
	  Y2[K] = Low+K;
	}
	IterateOverLineSearches(X1,Y1,X2,Y2,Num2,nfs,ndimx,ndimy,
				lowlim,PSatTrackSetProp);
      }
    }
  free(X1);
  free(Y1);
  free(X2);
  free(Y2);
  free(f);
  free(ind);
  freeFloatArray(nfs,ndimx,ndimy);
}

extern void Extrap(int *PX1, int *PY1, int X2, int Y2, int NX1)
{
  *PY1 = MakeInteger(*PY1 + (NX1 - *PX1) * (*PY1 - Y2) / (*PX1 - X2)+0.5);
  *PX1 = NX1;
}

extern void removev(char *PathName)
{
  if (remove(PathName) == -1) {
    perror(PathName);
    exit(1);
  }
}

extern int DoesFileExist(char *FN)
{
  int NotType;
  FILE *fin;

  NotType = 0;
  fin = fopen(FN,"r");
  if (fin != NULL) {
    fclose(fin);
  }
  else NotType = 1;
  
  return !NotType;
}

extern void OutputSatTrackSetProp(char *FN, int dimx, int dimy, 
				  SatTrackSetRec *PSatTrackSetProp, int WhExt, int Index)
{
  int I;
  char NFN[200],BaseFN[200];
  FILE *fout1,*fout2;

  if (PSatTrackSetProp->Num==0)
    return;
  strcpy(BaseFN,FN);
  BaseFN[strlen(BaseFN)-5] = 0;
  sprintf(NFN,"%s_SCI_%i.reg",BaseFN,WhExt);
  if (!DoesFileExist(NFN)) {
    fout1 = fopen(NFN,"w");
    fprintf(fout1,"global color=green font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n");
  }
  else {
    fout1 = fopen(NFN,"a");
  }
  sprintf(NFN,"%s.mask",BaseFN);
  fout2 = fopen(NFN,"a");

  for (I=0;I<PSatTrackSetProp->Num;I++) {
    SatTrackRec *PSatTrackProp;
    int X[4],Y[4];

    PSatTrackProp = &PSatTrackSetProp->PSatTrackProp[I];
    PSatTrackProp->X1 = PSatTrackProp->X1*4 + 2;
    PSatTrackProp->Y1 = PSatTrackProp->Y1*4 + 2;
    PSatTrackProp->X2 = PSatTrackProp->X2*4 + 2;
    PSatTrackProp->Y2 = PSatTrackProp->Y2*4 + 2;
    PSatTrackProp->Width *= 4;

    if (PSatTrackProp->X1 < 3)
      Extrap(&PSatTrackProp->X1,&PSatTrackProp->Y1,
	     PSatTrackProp->X2,PSatTrackProp->Y2,-1);
    if (PSatTrackProp->X2 < 3)
      Extrap(&PSatTrackProp->X2,&PSatTrackProp->Y2,
	     PSatTrackProp->X1,PSatTrackProp->Y1,-1);
    if (PSatTrackProp->Y1 < 3)
      Extrap(&PSatTrackProp->Y1,&PSatTrackProp->X1,
	     PSatTrackProp->Y2,PSatTrackProp->X2,-1);
    if (PSatTrackProp->Y2 < 3)
      Extrap(&PSatTrackProp->Y2,&PSatTrackProp->X2,
	     PSatTrackProp->Y1,PSatTrackProp->X1,-1);
    if (PSatTrackProp->X1 > dimx - 7)
      Extrap(&PSatTrackProp->X1,&PSatTrackProp->Y1,
	     PSatTrackProp->X2,PSatTrackProp->Y2,dimx);
    if (PSatTrackProp->X2 > dimx - 7)
      Extrap(&PSatTrackProp->X2,&PSatTrackProp->Y2,
	     PSatTrackProp->X1,PSatTrackProp->Y1,dimx);
    if (PSatTrackProp->Y1 > dimy - 7)
      Extrap(&PSatTrackProp->Y1,&PSatTrackProp->X1,
	     PSatTrackProp->Y2,PSatTrackProp->X2,dimy);
    if (PSatTrackProp->Y2 > dimy - 7)
      Extrap(&PSatTrackProp->Y2,&PSatTrackProp->X2,
	     PSatTrackProp->Y1,PSatTrackProp->X1,dimy);    

    if (fabs(PSatTrackProp->X2-PSatTrackProp->X1)<fabs(PSatTrackProp->Y2-PSatTrackProp->Y1)) {            
      WidthLXY(PSatTrackProp->X1,PSatTrackProp->Y1,PSatTrackProp->X2,PSatTrackProp->Y2,dimx,dimy,
	       &X[0],&Y[0],PSatTrackProp->Width);
      WidthHXY(PSatTrackProp->X1,PSatTrackProp->Y1,PSatTrackProp->X2,PSatTrackProp->Y2,dimx,dimy,
	       &X[1],&Y[1],PSatTrackProp->Width);
      WidthHXY(PSatTrackProp->X2,PSatTrackProp->Y2,PSatTrackProp->X1,PSatTrackProp->Y1,dimx,dimy,
	       &X[2],&Y[2],PSatTrackProp->Width);
      WidthLXY(PSatTrackProp->X2,PSatTrackProp->Y2,PSatTrackProp->X1,PSatTrackProp->Y1,dimx,dimy,
	       &X[3],&Y[3],PSatTrackProp->Width);
    }
    else {
      WidthLXY(PSatTrackProp->Y1,PSatTrackProp->X1,PSatTrackProp->Y2,PSatTrackProp->X2,dimy,dimx,
	       &Y[0],&X[0],PSatTrackProp->Width);
      WidthHXY(PSatTrackProp->Y1,PSatTrackProp->X1,PSatTrackProp->Y2,PSatTrackProp->X2,dimy,dimx,
	       &Y[1],&X[1],PSatTrackProp->Width);
      WidthHXY(PSatTrackProp->Y2,PSatTrackProp->X2,PSatTrackProp->Y1,PSatTrackProp->X1,dimy,dimx,
	       &Y[2],&X[2],PSatTrackProp->Width);
      WidthLXY(PSatTrackProp->Y2,PSatTrackProp->X2,PSatTrackProp->Y1,PSatTrackProp->X1,dimy,dimx,
	       &Y[3],&X[3],PSatTrackProp->Width);
    }
    fprintf(fout1,"physical;polygon(%i,%i,%i,%i,%i,%i,%i,%i)\n",
	    X[0]+1,Y[0]+1,X[1]+1,Y[1]+1,X[2]+1,Y[2]+1,X[3]+1,Y[3]+1);
    fprintf(fout2,"%s_inmask%i.fits 4 %i %i %i %i %i %i %i %i\n",
	    BaseFN,WhExt,X[0],Y[0],X[1],Y[1],X[2],Y[2],X[3],Y[3]);	    
  }
  fclose(fout1);
  fclose(fout2);
}

extern void DestructSatTrackSetProp(SatTrackSetRec *PSatTrackSetProp)
{
  if (PSatTrackSetProp->Num)
    free(PSatTrackSetProp->PSatTrackProp);
}

main(int argc, char *argv[])
{
  SatTrackSetRec SatTrackSetProp;
  int dimx,dimy,Index;
  float **fs;

  Globalidum = -1;
  for (Index=0;Index<2;Index++) {
    SReadFits(argv[1],&fs,&dimx,&dimy,1+Index*3);
    IdentifyLines(fs,dimx,dimy,&SatTrackSetProp);
    OutputSatTrackSetProp(argv[1],dimx,dimy,&SatTrackSetProp,Index+1,Index);
    freeFloatArray(fs,dimx,dimy);
    DestructSatTrackSetProp(&SatTrackSetProp);
  }
}
