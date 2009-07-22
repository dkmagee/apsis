#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "longnam.h"
#include "superalign.h"

float ran3(long *idum);

typedef struct {
  char FN[200];
  int DimX,DimY;
  float CX,CY,Angle;
  short **mask;
} MaskRec;

typedef struct {
  int NMask;
  MaskRec *PMaskProp;
} MaskListRec;

typedef struct {
  int ID;
  int LI;
  int NMI;
  float *PX,*PY;
  int NVert;
} PolygonRec;

extern void DestructPolygonProp(PolygonRec *PPolygonProp)
{
  free(PPolygonProp->PX);
  free(PPolygonProp->PY);
}

extern int WithinPolygon(PolygonRec *PPolygonProp,
			 float X, float Y)
{
  int I,Inside;
  float OAng,NAng,DiffAng,Rot;

  Rot = 0.0;
  for (I=1;I<PPolygonProp->NVert+1;I++) {
    if (I==1) {
      OAng = atan2(PPolygonProp->PY[0]-Y,PPolygonProp->PX[0]-X);
    }
    if (I<PPolygonProp->NVert)
      NAng = atan2(PPolygonProp->PY[I]-Y,PPolygonProp->PX[I]-X);
    else
      NAng = atan2(PPolygonProp->PY[0]-Y,PPolygonProp->PX[0]-X);

    DiffAng = NAng-OAng;
    if (DiffAng < 0) DiffAng += 2*M_PI;
    if (DiffAng > M_PI)
      Rot += (DiffAng - 2*M_PI);
    else
      Rot += DiffAng;
    OAng = NAng;
  }
  if ((Rot > M_PI) || (Rot < -M_PI))
    Inside = 1;
  else 
    Inside = 0;

  return Inside;
}

float X[2][12] = {{188,2538,4234,4200,4176,1883,87,134},
		  {307,1289,2637,3196,4295,4259,4235,3287,1995,1213,186,246}};
float Y[2][12] = {{2081,2177,2231,1085,153,104,55,1096},
		  {4094,4154,4224,4249,4292,3160,2286,2258,2213,2182,2134,3167}};

extern void SetUpPolygonProp(PolygonRec *PPolygonProp, MaskRec *PMaskProp,
			     int Wh)
{
  int I,CX,CY,TCX,TCY;
  float CosAng,SinAng;

  CosAng = cos(PMaskProp->Angle*M_PI/180);
  SinAng = sin(PMaskProp->Angle*M_PI/180);

  if (Wh==0)
    PPolygonProp->NVert = 8;
  else
    PPolygonProp->NVert = 12;

  TCX = 2178;
  TCY = 2184;
  PPolygonProp->PX = (float *)calloc(PPolygonProp->NVert,sizeof(float));
  PPolygonProp->PY = (float *)calloc(PPolygonProp->NVert,sizeof(float));
  for (I=0;I<PPolygonProp->NVert;I++) {
    PPolygonProp->PX[I] = PMaskProp->CX + 
      (X[Wh][I]-TCX)*CosAng/20 - (Y[Wh][I]-TCY)*SinAng/20;
    PPolygonProp->PY[I] = PMaskProp->CY +
      (X[Wh][I]-TCX)*SinAng/20 + (Y[Wh][I]-TCY)*CosAng/20;
  }
}

extern void printerror( int status)
{

  if (status) {
    fits_report_error(stderr, status); 
    
    exit( status ); 
  }
  return;
}

extern void WriteShortFits(char *FN, short **tfs, int dimx, int dimy)
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

extern void GenerateMask(MaskRec *PMaskProp)
{
  int I,J,OK;
  PolygonRec PolygonProp1,PolygonProp2;
  
  SetUpPolygonProp(&PolygonProp1,PMaskProp,0);
  SetUpPolygonProp(&PolygonProp2,PMaskProp,1);
  allocShortArray(&PMaskProp->mask,PMaskProp->DimX,PMaskProp->DimY);
  for (I=0;I<PMaskProp->DimX;I++)
    for (J=0;J<PMaskProp->DimY;J++) {
      OK = WithinPolygon(&PolygonProp1,I+0.5,J+0.5);
      if (!OK)
	OK = WithinPolygon(&PolygonProp2,I+0.5,J+0.5);
  
      PMaskProp->mask[J][I] = OK;
    }
}

extern void AddToStack(short **stack, MaskRec *PMaskProp)
{
  short **mask;
  int DimX,DimY,I,J;

  mask = PMaskProp->mask;
  DimX = PMaskProp->DimX;
  DimY = PMaskProp->DimY;
  for (I=0;I<DimX;I++)
    for (J=0;J<DimY;J++)
      stack[J][I] += mask[J][I];
}

extern void RemoveFromStack(short **stack, MaskRec *PMaskProp)
{
  short **mask;
  int DimX,DimY,I,J;

  mask = PMaskProp->mask;
  DimX = PMaskProp->DimX;
  DimY = PMaskProp->DimY;
  for (I=0;I<DimX;I++)
    for (J=0;J<DimY;J++)
      stack[J][I] -= mask[J][I];
}

extern float ScoreStack(short **stack, int DimX, int DimY, int Min)
{
  int I,J,Score,Max;
  float FScore;

  Score = 0;
  for (I=0;I<DimX;I++)
    for (J=0;J<DimY;J++) {
      if (stack[J][I] > Min)
	Max = Min;
      else
	Max = stack[J][I];
      Score += Max;
    }
  FScore = (float)Score / (DimX*DimY);
  return FScore;
}

extern void FindSolution(MaskListRec *PMaskListProp, int *PUse, int Min)
{
  short **stack;
  int DimX,DimY,I;
  float BScore;
  long Globalidum;
  
  Globalidum = -5;
  DimX = PMaskListProp->PMaskProp[0].DimX;
  DimY = PMaskListProp->PMaskProp[0].DimY;
  allocShortArray(&stack,DimX,DimY);

  BScore = 0;
  for (I=0;I<6000;I++) {
    MaskRec *PMaskProp,*PMaskPropS;
    float NScore,NScoreS;
    int Wh,WhS;

    Wh = (int)(ran3(&Globalidum)*PMaskListProp->NMask);
    WhS = (int)(ran3(&Globalidum)*PMaskListProp->NMask);
    PMaskProp = &PMaskListProp->PMaskProp[Wh];
    PMaskPropS = &PMaskListProp->PMaskProp[WhS];
    /*    WriteShortFits("TOT",stack,PMaskProp->DimX,PMaskProp->DimY);
    WriteShortFits("A",PMaskProp->mask,PMaskProp->DimX,PMaskProp->DimY);
    WriteShortFits("S",PMaskPropS->mask,PMaskProp->DimX,PMaskProp->DimY);*/
    if (!PUse[Wh]) {      
      AddToStack(stack,PMaskProp);
      NScore = ScoreStack(stack,DimX,DimY,Min);
      RemoveFromStack(stack,PMaskProp);
      if (PUse[WhS]) {
	RemoveFromStack(stack,PMaskPropS);
	NScoreS = ScoreStack(stack,DimX,DimY,Min);
	AddToStack(stack,PMaskPropS);
	if (NScoreS > BScore) {
	  BScore = NScoreS;
	  PUse[Wh] = 1;
	  PUse[WhS] = 0;
	  AddToStack(stack,PMaskProp);
	  RemoveFromStack(stack,PMaskPropS);
	}
	else if (NScore > BScore) {
	  BScore = NScore;
	  PUse[Wh] = 1;
	  AddToStack(stack,PMaskProp);
	}
	else if (NScoreS > BScore-1e-7) {
	  BScore = NScoreS;
	  PUse[WhS] = 0;
	  RemoveFromStack(stack,PMaskPropS);
	}
      }
      else if (NScore > BScore) {
	BScore = NScore;
	AddToStack(stack,PMaskProp);
	PUse[Wh] = 1;
      }      
    }      
  }
}

extern void ReadMaskList(char *FN, MaskListRec *PMaskListProp, 
			 int DimX, int DimY)
{
  FILE *fin;
  int I;

  fin = fopen(FN,"r");
  fscanf(fin,"%i%[^\n]",&PMaskListProp->NMask);
  fgetc(fin);

  PMaskListProp->NMask /= 2;
  PMaskListProp->PMaskProp = (MaskRec *)
    calloc(PMaskListProp->NMask,sizeof(MaskRec));

  for (I=0;I<PMaskListProp->NMask;I++) {
    MaskRec *PMaskProp;
    int ch;
    
    PMaskProp = &PMaskListProp->PMaskProp[I];
    fscanf(fin,"%s %g %g %g%*[^\n]",
	   PMaskProp->FN,&PMaskProp->CX,&PMaskProp->CY,&PMaskProp->Angle);
    fgetc(fin);
    ch = 0;
    while (PMaskProp->FN[ch] != '[')
      ch++;
    PMaskProp->FN[ch] = 0;
    fscanf(fin,"%*[^\n]");
    fgetc(fin);
    PMaskProp->DimX = DimX;
    PMaskProp->DimY = DimY;
    GenerateMask(PMaskProp);
  }
}

extern void OutputSolution(char *FN, MaskListRec *PMaskListProp,
			   int *PWh)
{
  FILE *fout;
  int I;

  fout = fopen(FN,"w");
  for (I=0;I<PMaskListProp->NMask;I++) {
    MaskRec *PMaskProp;

    PMaskProp = &PMaskListProp->PMaskProp[I];
    if (PWh[I])
      fprintf(fout,"%s[sci,1]\n%s[sci,2]\n",PMaskProp->FN,PMaskProp->FN);
  }
  fclose(fout);
}

main(int argc, char *argv[])
{
  MaskListRec MaskListProp;
  int *PUse,DimX,DimY,Min;

  sscanf(argv[3],"%i",&DimX);
  sscanf(argv[4],"%i",&DimY);
  sscanf(argv[5],"%i",&Min);
  ReadMaskList(argv[1],&MaskListProp,DimX,DimY);
  PUse = (int *)calloc(MaskListProp.NMask,sizeof(int));
  FindSolution(&MaskListProp,PUse,Min);
  OutputSolution(argv[2],&MaskListProp,PUse);
}
