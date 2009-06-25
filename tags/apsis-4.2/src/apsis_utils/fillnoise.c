#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "longnam.h"
#include "arrays.h"

long Globalidum;

#define BLOCKSIZE 1e7

float ran3(long *idum);

extern void printerror( int status)
{

  if (status) {
    fits_report_error(stderr, status); 
    
    exit( status ); 
  }
  return;
}

extern void SReadFits(char *FN, int LX, int HX, int LY, int HY,
		      float ***ptfs, int *dimx, int *dimy,
		      int WhImage, int *pcrpix1, int *pcrpix2)
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
  char comment[500],TFN[200];

  sprintf(TFN,"%s[%i:%i,%i:%i]",FN,LX+1,HX,LY+1,HY);
  fits_open_file(&afptr, TFN, READONLY, &status);

  HDU = WhImage+1;

  if (fits_movabs_hdu(afptr,HDU,&hdutype,&status))
    printerror( status );
  if (fits_read_key(afptr,TINT,"CRPIX1",pcrpix1,comment,&status))
    printerror( status );
  if (fits_read_key(afptr,TINT,"CRPIX2",pcrpix2,comment,&status))
    printerror( status );

  fits_get_img_size(afptr, 2, anaxes, &status);
  *dimx = anaxes[0];
  *dimy = anaxes[1]; 
  fprintf(stderr,"%i %i\n",*dimx,*dimy);
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

extern void WriteFloatFits(char *FN, float **tfs, int dimx, int dimy,
			   int crpix1, int crpix2)
{
  int status,naxis,J;
  float *ffs;
  fitsfile *fptr;
  long fpixel,naxes[2];
  char OFN[200];

  status = 0;
  sprintf(OFN,"%s",FN);
  if (fits_open_image(&fptr,OFN,READWRITE,&status))
    printerror( status );

  naxis = 2;
  fpixel = 1;
  naxes[0] = dimx;
  naxes[1] = dimy;

  if (fits_update_key(fptr,TINT,"NAXIS1",&dimx,"",&status))
    printerror( status );
  if (fits_update_key(fptr,TINT,"NAXIS2",&dimy,"",&status))
    printerror( status );

  if (fits_update_key(fptr,TINT,"CRPIX1",&crpix1,"",&status))
    printerror( status );
  if (fits_update_key(fptr,TINT,"CRPIX2",&crpix2,"",&status))
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

extern void AddNoise(char *FN, char *WFN, char *SLX, char *SLY, char *SHX, char *SHY)
{
  float **fs,**wfs;
  int dimx,dimy,I,J,crpix1,crpix2,t1,t2,LX,LY,HX,HY;
  char BaseFN[100];

  sscanf(SLX,"%i",&LX); 
  sscanf(SHX,"%i",&HX);
  sscanf(SLY,"%i",&LY);
  sscanf(SHY,"%i",&HY);
  fprintf(stderr,"Read in %s.\n",FN);
  SReadFits(FN,LX,HX,LY,HY,&fs,&dimx,&dimy,0,&crpix1,&crpix2);
  fprintf(stderr,"Read in %s.\n",WFN);
  SReadFits(WFN,LX,HX,LY,HY,&wfs,&dimx,&dimy,0,&t1,&t2);
  fprintf(stderr,"Adding noise.\n");
  for (I=0;I<dimx;I++) {
    if (I%100==0) fprintf(stderr,"%i\n",I);
    for (J=0;J<dimy;J++)
      if (wfs[J][I] < 1.1e5) 
	fs[J][I] = ran3(&Globalidum)*0.003;
  }
  fprintf(stderr,"Writing out file.\n");
  WriteFloatFits(FN,fs,dimx,dimy,crpix1,crpix2);
}

main(int argc, char *argv[])
{
  Globalidum = -1;

  AddNoise(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
}
