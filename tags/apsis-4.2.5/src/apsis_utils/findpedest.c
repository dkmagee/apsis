#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "longnam.h"
#include "arrays.h"

#define BLOCKSIZE 1e7

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

extern void SexUpSexDefaults()
{
  FILE *fout;
  char BaseTempFN[200];

  sprintf(BaseTempFN,"default.sex");
  fout = fopen(BaseTempFN,"w");
  fprintf(fout,"CATALOG_NAME    test.cat\n");
  fprintf(fout,"CATALOG_TYPE    ASCII_HEAD\n");
  fprintf(fout,"PARAMETERS_NAME default.param\n");
  fprintf(fout,"DETECT_TYPE     CCD\n");
  fprintf(fout,"FLAG_IMAGE      flag.fits\n");
  fprintf(fout,"DETECT_MINAREA  30\n");
  fprintf(fout,"DETECT_THRESH   200.\n");
  fprintf(fout,"ANALYSIS_THRESH 200.\n");
  fprintf(fout,"FILTER          Y\n");
  fprintf(fout,"FILTER_NAME     default.conv\n");
  fprintf(fout,"DEBLEND_NTHRESH 32\n");
  fprintf(fout,"DEBLEND_MINCONT 0.005\n");
  fprintf(fout,"CLEAN           Y\n");
  fprintf(fout,"CLEAN_PARAM     1.0\n");
  fprintf(fout,"MASK_TYPE       CORRECT\n");
  fprintf(fout,"PHOT_APERTURES  5\n");
  fprintf(fout,"PHOT_AUTOPARAMS 2.5, 3.5\n");
  fprintf(fout,"SATUR_LEVEL     50000.0\n");
  fprintf(fout,"MAG_ZEROPOINT   0.0\n");
  fprintf(fout,"MAG_GAMMA       4.0\n");
  fprintf(fout,"GAIN            0.0\n");
  fprintf(fout,"PIXEL_SCALE     1.0\n");
  fprintf(fout,"SEEING_FWHM     1.2\n");
  fprintf(fout,"STARNNW_NAME    default.nnw\n");
  fprintf(fout,"BACK_SIZE       128\n");
  fprintf(fout,"BACK_FILTERSIZE 8\n");
  fprintf(fout,"BACKPHOTO_TYPE  GLOBAL\n");
  fprintf(fout,"CHECKIMAGE_TYPE BACKGROUND\n");
  fprintf(fout,"CHECKIMAGE_NAME check.fits\n");
  fprintf(fout,"MEMORY_OBJSTACK 2000\n");
  fprintf(fout,"MEMORY_PIXSTACK 200000\n");
  fprintf(fout,"MEMORY_BUFSIZE  1024\n");
  fprintf(fout,"VERBOSE_TYPE    NORMAL\n");
  fclose(fout);
  sprintf(BaseTempFN,"default.param");
  fout = fopen(BaseTempFN,"w");
  fprintf(fout,"NUMBER\n");
  fprintf(fout,"MAG_ISO\n");
  fclose(fout);
  sprintf(BaseTempFN,"default.conv");
  fout = fopen(BaseTempFN,"w");
  fprintf(fout,"CONV NORM\n");
  fprintf(fout,"1 2 1\n");
  fprintf(fout,"2 4 2\n");
  fprintf(fout,"1 2 1\n");
  fclose(fout);
}

extern void EstimateBackground(float **fs, int dimx, int dimy,
			       float ***pbfs)
{
  char CMD[200];
  char BaseTempFN[200];
  int tdimx,tdimy;

  SexUpSexDefaults();
  WriteFloatFits("t.fits",fs,dimx,dimy);
  
  sprintf(CMD,"sex t.fits\n");
  system(CMD);
  sprintf(BaseTempFN,"check.fits");
  SReadFits(BaseTempFN,pbfs,&tdimx,&tdimy,0);
  remove("check.fits");
  remove("t.fits");
  remove("test.cat");
  remove("default.sex");
  remove("default.param");
  remove("default.conv");
}

long Globalidum;

#define BLOCKSIZE 1e7

#define PIX_SORT(a,b) {if (a>b) {float temp=a;a=b;b=temp;}}

float opt_med9(float *p)
{
  PIX_SORT(p[1], p[2]); PIX_SORT(p[4], p[5]); PIX_SORT(p[7], p[8]);
  PIX_SORT(p[0], p[1]); PIX_SORT(p[3], p[4]); PIX_SORT(p[6], p[7]);
  PIX_SORT(p[1], p[2]); PIX_SORT(p[4], p[5]); PIX_SORT(p[7], p[8]);
  PIX_SORT(p[0], p[3]); PIX_SORT(p[5], p[8]); PIX_SORT(p[4], p[7]);
  PIX_SORT(p[3], p[6]); PIX_SORT(p[1], p[4]); PIX_SORT(p[2], p[5]);
  PIX_SORT(p[4], p[7]); PIX_SORT(p[4], p[2]); PIX_SORT(p[6], p[4]);
  PIX_SORT(p[4], p[2]); return(p[4]);
}

float opt_med25(float *p)
{

  PIX_SORT(p[0], p[1]); PIX_SORT(p[3], p[4]); PIX_SORT(p[2], p[4]);
  PIX_SORT(p[2], p[3]); PIX_SORT(p[6], p[7]); PIX_SORT(p[5], p[7]);
  PIX_SORT(p[5], p[6]); PIX_SORT(p[9], p[10]); PIX_SORT(p[8], p[10]);
  PIX_SORT(p[8], p[9]); PIX_SORT(p[12], p[13]); PIX_SORT(p[11], p[13]);
  PIX_SORT(p[11], p[12]); PIX_SORT(p[15], p[16]); PIX_SORT(p[14], p[16]);
  PIX_SORT(p[14], p[15]); PIX_SORT(p[18], p[19]); PIX_SORT(p[17], p[19]);
  PIX_SORT(p[17], p[18]); PIX_SORT(p[21], p[22]); PIX_SORT(p[20], p[22]);
  PIX_SORT(p[20], p[21]); PIX_SORT(p[23], p[24]); PIX_SORT(p[2], p[5]);
  PIX_SORT(p[3], p[6]); PIX_SORT(p[0], p[6]); PIX_SORT(p[0], p[3]);
  PIX_SORT(p[4], p[7]); PIX_SORT(p[1], p[7]); PIX_SORT(p[1], p[4]);
  PIX_SORT(p[11], p[14]); PIX_SORT(p[8], p[14]); PIX_SORT(p[8], p[11]);
  PIX_SORT(p[12], p[15]); PIX_SORT(p[9], p[15]); PIX_SORT(p[9], p[12]);
  PIX_SORT(p[13], p[16]); PIX_SORT(p[10], p[16]); PIX_SORT(p[10], p[13]);
  PIX_SORT(p[20], p[23]); PIX_SORT(p[17], p[23]); PIX_SORT(p[17], p[20]);
  PIX_SORT(p[21], p[24]); PIX_SORT(p[18], p[24]); PIX_SORT(p[18], p[21]);
  PIX_SORT(p[19], p[22]); PIX_SORT(p[8], p[17]); PIX_SORT(p[9], p[18]);
  PIX_SORT(p[0], p[18]); PIX_SORT(p[0], p[9]); PIX_SORT(p[10], p[19]);
  PIX_SORT(p[1], p[19]); PIX_SORT(p[1], p[10]); PIX_SORT(p[11], p[20]);
  PIX_SORT(p[2], p[20]); PIX_SORT(p[2], p[11]); PIX_SORT(p[12], p[21]);
  PIX_SORT(p[3], p[21]); PIX_SORT(p[3], p[12]); PIX_SORT(p[13], p[22]);
  PIX_SORT(p[4], p[22]); PIX_SORT(p[4], p[13]); PIX_SORT(p[14], p[23]);
  PIX_SORT(p[5], p[23]); PIX_SORT(p[5], p[14]); PIX_SORT(p[15], p[24]);
  PIX_SORT(p[6], p[24]); PIX_SORT(p[6], p[15]); PIX_SORT(p[7], p[16]);
  PIX_SORT(p[7], p[19]); PIX_SORT(p[13], p[21]); PIX_SORT(p[15], p[23]);
  PIX_SORT(p[7], p[13]); PIX_SORT(p[7], p[15]); PIX_SORT(p[1], p[9]);
  PIX_SORT(p[3], p[11]); PIX_SORT(p[5], p[17]); PIX_SORT(p[11], p[17]);
  PIX_SORT(p[9], p[17]); PIX_SORT(p[4], p[10]); PIX_SORT(p[6], p[12]);
  PIX_SORT(p[7], p[14]); PIX_SORT(p[4], p[6]); PIX_SORT(p[4], p[7]);
  PIX_SORT(p[12], p[14]); PIX_SORT(p[10], p[14]); PIX_SORT(p[6], p[7]);
  PIX_SORT(p[10], p[12]); PIX_SORT(p[6], p[10]); PIX_SORT(p[6], p[17]);
  PIX_SORT(p[12], p[17]); PIX_SORT(p[7], p[17]); PIX_SORT(p[7], p[10]);
  PIX_SORT(p[12], p[18]); PIX_SORT(p[7], p[12]); PIX_SORT(p[10], p[18]);
  PIX_SORT(p[12], p[20]); PIX_SORT(p[10], p[20]); PIX_SORT(p[10], p[12]);

  return (p[12]);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float selectk(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i,ir,j,l,mid;
  float a,temp;
  
  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	  }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  }
      arr[l]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef SWAP

float OvalShape[9][2] = {{-21,25},{-1,25},{20,8},{28,-18},{17,-28},
			 {4,-27},{-8,-20},{-29,10},{-21,25}};

extern int MakeInteger(float F)
{
  int V;
  if (F >= 0) V = (int)F;
  else V = (int)F - 1;
  return V;
}

extern void OvalScore(float **fs, int dimx, int dimy, int X, int Y, 
		      float scale, int ndens, int xcut, int ycut,
		      float clip, float rms,
		      float *pmed, float *plowquart, float *pmode)
{
  int Tot,nring,ndiam,J,I,CX,CY,OK,*ind,mod[100],max,wh,whmax;
  float *f,ringscale,V,dX,dY,AvgX,AvgY;

  ndiam = (int)sqrt(ndens / 8);
  if (ndiam < 1) ndiam = 1;
  nring = ndens / ndiam;

  f = (float *)calloc(ndiam*nring,sizeof(float));
  ind = (int *)calloc(ndiam*nring,sizeof(int));
  Tot = 0;
  if (ndens > 1) {
    for (I=0;I<100;I++)
      mod[I] = 0;
    max = 0;
    whmax = 0;
  }
    
  for (J=0;J<ndiam;J++) {
    ringscale = (0.8 + 0.5*(J-ndiam/2.+0.5)/ndiam);
    for (I=0;I<nring;I++) {
      OK = 1;
      V = ((float)I*8 / nring + (float)J/ndiam);
      if (V > 8)
	V -= 8;
      dX = ((int)V+1-V)*OvalShape[(int)V][0] + 
	(V-(int)V)*OvalShape[(int)V+1][0];
      dY = ((int)V+1-V)*OvalShape[(int)V][1] + 
	(V-(int)V)*OvalShape[(int)V+1][1];

      if ((xcut==1)&&(dX>0)||(xcut==2)&&(dX<0))
	OK = 0;
      if ((ycut==1)&&(dY>0)||(ycut==2)&&(dY<0))
	OK = 0;
      if (OK) {
	CX = MakeInteger(X + dX * scale * ringscale+0.5);
	CY = MakeInteger(Y + dY * scale * ringscale+0.5);
	if (clip < 1)
	  fs[CY][CX] = 1e4;
	if ((CX >= 0) && (CX < dimx) && (CY >= 0) && (CY < dimy)) {
	  if (fs[CY][CX] < clip) {
	    if (ndens > 1) {
	      wh = MakeInteger(fs[CY][CX] * 5 / rms + 10.5);
	      if (wh > 99) wh = 99;
	      if (wh >= 0) {
		mod[wh]++;
		if (mod[wh] > max) {
		  max = mod[wh];
		  whmax = wh;
		}
	      }
	    }
	    f[Tot] = fs[CY][CX];
	    Tot++;
	  }
	}
      }
    }
  }  
  hpsort(Tot,f,ind);
  V = (Tot-1)/2;
  *pmed = ((int)V+1-V)*f[ind[(int)V]] + (V-(int)V)*f[ind[(int)V]];
  V = (Tot-1)/4;
  *plowquart = ((int)V+1-V)*f[ind[(int)V]] + (V-(int)V)*f[ind[(int)V]];  
  if (ndens > 1) {
    *pmode = rms*(whmax-10)/5.0;
  }
  free(f);
  free(ind);
}

extern float InnerOverSearch(float **nfs, int ndimx, int ndimy,
			     int X, int Y)
{
  float Sum,Radius;
  int CX,CY,DX,DY;
  float NDX,NDY;

  Sum = 0.0;

  for (CX=X-10;CX<X+10;CX++)
    for (CY=Y-9;CY<Y+9;CY++) {
      DX = CX-X;
      DY = CY-Y;
      NDX = DX*0.642 + DY*0.766;
      NDY = -DX*0.766 + DY*0.642;
      Radius = (NDX*NDX + (float)NDY*NDY/9);
      if ((CX >= 0) && (CX < ndimx) && (CY >= 0) && (CY < ndimy) 
	  && (Radius < 36)) {
	Sum += nfs[CY][CX];
      }
    }
  return Sum;
}

extern void ClearNeigh(float **nfs, int ndimx, int ndimy,
		       int X, int Y)
{
  int CX,CY;

  for (CX=X-50;CX<X+50;CX++)
    for (CY=Y-50;CY<Y+50;CY++) {
      if ((CX >= 0) && (CX < ndimx) && (CY >= 0) && (CY < ndimy)) {
	nfs[CY][CX] = 0;
      }
    }
}

typedef struct {
  float X,Y;
} OvalRec;

typedef struct {
  int Num;
  OvalRec *POvalProp;
} OvalSetRec;

extern void SearchForOvalPlateau(float **nfs, int ndimx, int ndimy,
				 float rms, OvalSetRec *POvalSetProp)
{
  int X,Y;
  float **medfs,**difffs,med,lowquart,**goodfs,**derivfs,**modefs;
  float medx0,medx1,medy0,medy1,lowq,divf,rms3,rms4,mode;

  rms3 = rms * 1.4 / 5.6;
  rms4 = rms * 7 / 5.6;
  allocFloatArray(&medfs,ndimx,ndimy);
  allocFloatArray(&difffs,ndimx,ndimy); 
  allocFloatArray(&derivfs,ndimx,ndimy); 
  allocFloatArray(&modefs,ndimx,ndimy);
  allocFloatArray(&goodfs,ndimx,ndimy);
  fprintf(stderr,"Starting first pass.\n");
  for (X=0;X<ndimx;X++)
    for (Y=0;Y<ndimy;Y++) {
      OvalScore(nfs,ndimx,ndimy,X,Y,1.0,16,0,0,1e10,rms,&med,&lowquart,&mode);
      medfs[Y][X] = med;
      modefs[Y][X] = mode;
      divf = 1;
      if (med > rms4) divf = pow(med/rms4,0.25);
      difffs[Y][X] = (lowquart - med)/divf;
      if ((medfs[Y][X] > rms3) && (difffs[Y][X] > -rms4))
	goodfs[Y][X] = 1;
    }
  /*WriteFloatFits("B",medfs,ndimx,ndimy);  
  WriteFloatFits("C",difffs,ndimx,ndimy);
  WriteFloatFits("M",modefs,ndimx,ndimy);*/

  fprintf(stderr,"First pass done.\n");
  fprintf(stderr,"Starting second pass.\n");
  for (X=0;X<ndimx;X++)
    for (Y=0;Y<ndimy;Y++) {
      if (goodfs[Y][X] > 0) {
	OvalScore(nfs,ndimx,ndimy,X,Y,1.0,256,0,0,medfs[Y][X]+rms*3,rms,
		  &med,&lowquart,&mode);
	divf = 1;
	if (med > rms4) divf = pow(med/rms4,0.25);
	if ((med > rms3) && ((lowquart-med)/divf > -rms4) &&
	    (mode > rms3)) {
	  goodfs[Y][X] = 1;
	  
	  OvalScore(nfs,ndimx,ndimy,X,Y,1.0,64,1,0,1e7,rms,&medx0,&lowq,
		    &mode);
	  OvalScore(nfs,ndimx,ndimy,X,Y,1.0,64,2,0,1e7,rms,&medx1,&lowq,
		    &mode);
	  OvalScore(nfs,ndimx,ndimy,X,Y,1.0,64,0,1,1e7,rms,&medy0,&lowq,
		    &mode);
	  OvalScore(nfs,ndimx,ndimy,X,Y,1.0,64,0,2,1e7,rms,&medy1,&lowq,
		    &mode);
	  derivfs[Y][X] = sqrt((medx0 - medx1)*(medx0 - medx1) +
			       (medy0 - medy1)*(medy0 - medy1))/med;
	  if (derivfs[Y][X] > 0.55)
	    goodfs[Y][X] = 0;
	}
	else
	  goodfs[Y][X] = 0;
      }
    }
  /*  OvalScore(goodfs,ndimx,ndimy,545,513,1.0,256,0,0,-1,rms,&medy1,&lowq,
      &mode);*/

  /*  WriteFloatFits("E",goodfs,ndimx,ndimy);*/
  POvalSetProp->Num = 0;
  if (1) {
    int Done;
    
    Done = 0;
    while (!Done) {
      int BestX,BestY;
      float Best,Cur;
      
      Best = 0;
      for (X=0;X<ndimx;X++)
	for (Y=0;Y<ndimy;Y++) {
	  if (goodfs[Y][X] > 0) {
	    Cur = InnerOverSearch(goodfs,ndimx,ndimy,X,Y);
	    if (Cur > Best) {
	      Best = Cur;
	      BestX = X;
	      BestY = Y;
	    }
	  }
	}
      if (Best > 40) {
	OvalRec *POvalProp;
	if (POvalSetProp->Num==0)
	  POvalSetProp->POvalProp = (OvalRec *)
	    calloc((POvalSetProp->Num+1),sizeof(OvalRec));
	else
	  POvalSetProp->POvalProp = (OvalRec *)
	    realloc(POvalSetProp->POvalProp,
		    (POvalSetProp->Num+1)*sizeof(OvalRec));
	POvalProp = &POvalSetProp->POvalProp[POvalSetProp->Num];
	POvalProp->X = BestX;
	POvalProp->Y = BestY;

	POvalSetProp->Num++;
	fprintf(stderr,"%g %i %i\n",Best,BestX,BestY);
	ClearNeigh(goodfs,ndimx,ndimy,BestX,BestY);
      }
      else Done = 1;
    }
  }
}

extern void IdentifyPedest(char *FN, float **fs, int dimx, int dimy,
			   OvalSetRec *POvalSetProp)
{
  float *f,**nfs,avg,rms,med,lowv;
  int ndimx,ndimy,I,J,I0,J0,Tot;
  float **bfs;

  f = (float *)calloc(81,sizeof(float));  
  ndimx = dimx / 3;
  ndimy = dimy / 3;
  EstimateBackground(fs,dimx,dimy,&bfs);
  avg = 0;
  for (I=0;I<dimx;I++)
    for (J=0;J<dimy;J++) {
      avg += bfs[J][I];
      fs[J][I] -= bfs[J][I];
    }
  avg /= (dimx*dimy);
  rms = sqrt(avg);
  freeFloatArray(bfs,dimx,dimy);

  allocFloatArray(&nfs,ndimx,ndimy);
  for (I=0;I<ndimx;I++)
    for (J=0;J<ndimy;J++) {
      for (I0=0;I0<3;I0++)
	for (J0=0;J0<3;J0++)
	  f[I0*3+J0] = fs[J*3+J0][I*3+I0];
      nfs[J][I] = opt_med9(f);
    }
  fprintf(stderr,"rms=%.3f\n",rms);
  /*  WriteFloatFits("A",nfs,ndimx,ndimy);  */
  /*  SearchForSquarePlateau(nfs,ndimx,ndimy,9,rms);*/
  SearchForOvalPlateau(nfs,ndimx,ndimy,rms,POvalSetProp);
  freeFloatArray(nfs,ndimx,ndimy);
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

extern void OutputOvalSetProp(char *FN, OvalSetRec *POvalSetProp,
			      int WhExt, int Index)
{
  int I,X,Y;
  char NFN[200],BaseFN[200];
  FILE *fout1,*fout2;

  if (POvalSetProp->Num==0)
    return;
  strcpy(BaseFN,FN);
  BaseFN[strlen(BaseFN)-5] = 0;
  sprintf(NFN,"%s_SCI_%i.reg",BaseFN,WhExt);
  if (!DoesFileExist(NFN)) {
    fout1 = fopen(NFN,"a");
    fprintf(fout1,"global color=green font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n");
  }
  else {
    fout1 = fopen(NFN,"w");    
  }
  sprintf(NFN,"%s.mask",BaseFN);
  fout2 = fopen(NFN,"a");

  for (I=0;I<POvalSetProp->Num;I++) {
    OvalRec *POvalProp;
    int J;

    POvalProp = &POvalSetProp->POvalProp[I];
    fprintf(fout1,"physical;polygon(");
    fprintf(fout2,"%s_inmask%i.fits 8",BaseFN,WhExt);
    for (J=0;J<8;J++) {
      X = MakeInteger((POvalProp->X + OvalShape[J][0]*1.5)*3+0.5);
      Y = MakeInteger((POvalProp->Y + OvalShape[J][1]*1.5)*3+0.5);
      fprintf(fout1,"%i,%i",X,Y);
      if (J<7)
	fprintf(fout1,",");
      else
	fprintf(fout1,")\n");
      fprintf(fout2," %i %i",X,Y);
    }
    fprintf(fout2,"\n");
  }
  fclose(fout1);
  fclose(fout2);
}

main(int argc, char *argv[])
{
  int dimx,dimy,Index;
  float **fs;

  Globalidum = -1;
  for (Index=0;Index<2;Index++) {
    OvalSetRec OvalSetProp;

    SReadFits(argv[1],&fs,&dimx,&dimy,1+Index*3);
    IdentifyPedest(argv[1],fs,dimx,dimy,&OvalSetProp);
    OutputOvalSetProp(argv[1],&OvalSetProp,Index+1,Index);
  }
}
