#define NFPRINTF(w,x) {fprintf(w, "\33[1M> %s\n\33[1A",x);}

#include <stdio.h>
#include <stdlib.h>
#include "arrays.h"
#include "longnam.h"
#include "fitsio.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;


extern void printerror( int status)
{

  if (status) {
    fits_report_error(stderr, status); 
    
    exit( status ); 
  }
  return;
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

typedef char FileStr[130];

typedef struct {
  int NExp,NHigh,Detail;
  float NLow;
  FileStr *PDz;
  FileStr *PCx;
} FileListRec;

void ReadFileList(char *FN, FileListRec *PFileListProp)
{
  FILE *fin;
  int I;

  fin = fopen(FN,"r");
  fscanf(fin,"%i %g %i%*[^\n]",&PFileListProp->NExp,
	 &PFileListProp->NLow,&PFileListProp->NHigh,
	 &PFileListProp->Detail);
  fgetc(fin);

  PFileListProp->PDz = (FileStr *)
    calloc(PFileListProp->NExp,sizeof(FileStr));
  PFileListProp->PCx = (FileStr *)
    calloc(PFileListProp->NExp,sizeof(FileStr));

  for (I=0;I<PFileListProp->NExp;I++) {
    fscanf(fin,"%s %[^\n]",PFileListProp->PDz[I],PFileListProp->PCx[I]);
    fgetc(fin);
  }
  fclose(fin);
}

main(int argc, char *argv[])
{
  FileListRec FileListProp;
  fitsfile **pfin,*fptr;
  int status = 0;
  long firstpixel;
  int hdutype,HDU;
  long anaxes[] = {0,0};
  long firstpix = 1, npixels;
  int I,J,Wh,K;
  int keysexist,morekeys;
  float **scidata,*oarray,*outline,*apix,*PExpTime;
  int **maskdata,Tot,N1,N2;
  char Card[90],Comment[100];

  fprintf(stdout,"Reading Input file '%s'.\n",argv[1]);
  ReadFileList(argv[1],&FileListProp);
  pfin = (fitsfile **)calloc(FileListProp.NExp*2,sizeof(fitsfile *));
  PExpTime = (float *)calloc(FileListProp.NExp,sizeof(float));

  fprintf(stdout,"Median Filtering...  NFile = %i, NLow = %.2f, NHigh = %i\n\n",
	  FileListProp.NExp,FileListProp.NLow,FileListProp.NHigh);

  for (I=0;I<FileListProp.NExp;I++) {
    fits_open_file(&pfin[2*I], FileListProp.PDz[I], READONLY, &status);
    if (fits_movabs_hdu(pfin[2*I],1,&hdutype,&status))
      printerror( status );
    if (fits_read_key(pfin[2*I],TFLOAT,"EXPTIME",&PExpTime[I],
		      Comment,&status)) {
      printerror( status );
    }
    if (I==0) {
      fits_get_img_size(pfin[0], 2, anaxes, &status);
    }

    fprintf(stdout,"%s %s %.1f %.3f\n",FileListProp.PDz[I],
	    FileListProp.PCx[I],PExpTime[I],PExpTime[I]/(PExpTime[0]+1e-15));
    if (fits_close_file(pfin[I*2],&status))
      printerror( status );      
  }

  if (DoesFileExist(argv[2]))
    remove(argv[2]);
  if (fits_create_file(&fptr,argv[2],&status))
    printerror( status );
  
  fprintf(stdout,"\nGenerating File '%s'\n",argv[2]);

  HDU = 1;
  if (fits_movabs_hdu(fptr,HDU,&hdutype,&status))
    printerror( status );
  if (fits_create_img(fptr,FLOAT_IMG,2,anaxes,&status))
    printerror( status );

  N1 = anaxes[0]*20;
  N2 = anaxes[0]*anaxes[1] / N1;
  oarray = (float *)calloc(FileListProp.NExp+1,sizeof(float));
  outline = (float *)calloc(N1,sizeof(float));
  allocFloatArray(&scidata,N1,FileListProp.NExp);
  allocIntArray(&maskdata,N1,FileListProp.NExp);
  apix = (float *)calloc(N1,sizeof(float));
  for (J=0;J<N2;J++) {
    int NL;
    char STR1[100];

    if (J%10==0) {
      sprintf(STR1,"%i%% Done.",J*100/N2);
      NFPRINTF(stdout,STR1);
    }
    NL = anaxes[0]*anaxes[1] - N1*J;
    if (NL > N1)
      NL = N1;
    for (I=0;I<FileListProp.NExp;I++) {
      firstpix = 1+N1*J;
      npixels = NL;

      fits_open_file(&pfin[I*2], FileListProp.PDz[I], READONLY, &status);
      if (fits_movabs_hdu(pfin[I*2],1,&hdutype,&status))
	printerror( status );
      fits_read_img(pfin[I*2], TFLOAT, firstpix, npixels,
		    NULL, apix, NULL, &status);
      for (K=0;K<NL;K++)
	scidata[I][K] = apix[K];
      if (fits_close_file(pfin[I*2],&status))
	printerror( status );      

      fits_open_file(&pfin[I*2+1], FileListProp.PCx[I], READONLY, &status);
      if (fits_movabs_hdu(pfin[I*2+1],1,&hdutype,&status))
	printerror( status );
      fits_read_img(pfin[I*2+1], TFLOAT, firstpix, npixels,
		    NULL, apix, NULL, &status);      
      for (K=0;K<NL;K++)
	maskdata[I][K] = (int)(apix[K]+0.5);
      if (fits_close_file(pfin[I*2+1],&status))
	printerror( status );      
    }
    for (K=0;K<NL;K++) {
      Tot = 0;
      for (I=0;I<FileListProp.NExp;I++) { 
	if (maskdata[I][K]) {
	  oarray[Tot+1] = scidata[I][K]*PExpTime[0]/(PExpTime[I]+1e-15);
	  Tot++;
	}
      }      
      if (Tot > 0) {
	if (FileListProp.Detail) {
	  float WhF;
	  
	  WhF = Tot*(float)FileListProp.NLow/FileListProp.NExp-0.5;
	  outline[K] = selectk((int)WhF+1,Tot,oarray) * (1 + (int)WhF - WhF);
	  if ((int)(WhF)+2 <= Tot)
	    outline[K] += selectk((int)WhF+2,Tot,oarray) * (WhF - (int)WhF);
	}
	else {
	  Wh = (int)(Tot*(float)FileListProp.NLow/FileListProp.NExp);
	  outline[K] = selectk(Wh+1,Tot,oarray);
	}
      }
      else
	outline[K] = 0.0;
    }    
    firstpixel = 1+J*N1;
    if (fits_write_img(fptr,TFLOAT,firstpixel,NL,outline,&status))
      printerror( status );
  }

  fits_open_file(&pfin[0], FileListProp.PDz[0], READONLY, &status);
  if (fits_movabs_hdu(pfin[0],1,&hdutype,&status))
    printerror( status );
  if (fits_get_hdrspace(pfin[0],&keysexist,&morekeys,&status))
    printerror( status );

  for (J=0;J<keysexist;J++) {
    if (fits_read_record(pfin[0],J+1,Card,&status))
      printerror( status );
    
    if (fits_write_record(fptr,Card,&status))
      printerror( status );
  }
  if (fits_close_file(pfin[0],&status))
    printerror( status );

  if (fits_close_file(fptr,&status))
    printerror( status );
}
