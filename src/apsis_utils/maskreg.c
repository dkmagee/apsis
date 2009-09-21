#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"
#include "longnam.h"
#include "arrays.h"

extern void printerror( int status)
{

  if (status) {
    fits_report_error(stderr, status); 
    
    exit( status ); 
  }
  return;
}

typedef struct {
  int NPoint;
  float *PX,*PY;
} MaskRegRec;

typedef struct {
  char Name[100];
  int NRegion;
  MaskRegRec *PMaskRegProp;
} ImageMaskRec;

typedef struct {
  int NImage;
  ImageMaskRec *PImageMaskProp;
} MaskRec;

typedef char LineStr[100];

void ReadFile(char *FN, MaskRec *PMaskProp)
{
  FILE *fin;  
  int Tot,NReal,I;
  char STR[200];

  Tot = 0;
  fin = fopen(FN,"r");
  while (fscanf(fin,"%[^\n]",STR) != EOF) {
    fgetc(fin);
    Tot++;
  }
  fclose(fin);

  NReal = 0;
  PMaskProp->PImageMaskProp = (ImageMaskRec *)
    calloc(Tot,sizeof(ImageMaskRec));
  NReal = -1;
  fin = fopen(FN,"r");
  while (fscanf(fin,"%s",STR) != EOF) {
    ImageMaskRec *PImageMaskProp;
    MaskRegRec *PMaskRegProp;

    if ((NReal<0)||
	(strcmp(PMaskProp->PImageMaskProp[NReal].Name,STR))) {
      NReal++;
      strcpy(PMaskProp->PImageMaskProp[NReal].Name,STR);
    }
    fgetc(fin);
    PImageMaskProp = &PMaskProp->PImageMaskProp[NReal];
    if (PImageMaskProp->NRegion == 0)
      PImageMaskProp->PMaskRegProp = (MaskRegRec *)
	calloc(1,sizeof(MaskRegRec));
    else
      PImageMaskProp->PMaskRegProp = (MaskRegRec *)
	realloc(PImageMaskProp->PMaskRegProp,
		(PImageMaskProp->NRegion+1)*sizeof(MaskRegRec));

    PMaskRegProp = &PImageMaskProp->PMaskRegProp[PImageMaskProp->NRegion];
    fscanf(fin,"%i",&PMaskRegProp->NPoint);
    PMaskRegProp->PX = (float *)
      calloc(PMaskRegProp->NPoint,sizeof(float));
    PMaskRegProp->PY = (float *)
      calloc(PMaskRegProp->NPoint,sizeof(float));

    for (I=0;I<PMaskRegProp->NPoint;I++)
      fscanf(fin,"%g %g",&PMaskRegProp->PX[I],&PMaskRegProp->PY[I]);
    fscanf(fin,"%*[^\n]");
    fgetc(fin);
    PImageMaskProp->NRegion++;
  }
  PMaskProp->NImage = NReal+1;
  fclose(fin);  
}

extern int WithinPolygon(MaskRegRec *PMaskRegProp,
			 float X, float Y)
{
  int I,Inside;
  float OAng,NAng,DiffAng,Rot;

  Rot = 0.0;
  for (I=1;I<PMaskRegProp->NPoint+1;I++) {
    if (I==1) {
      OAng = atan2(PMaskRegProp->PY[0]-Y,PMaskRegProp->PX[0]-X);
    }
    if (I<PMaskRegProp->NPoint)
      NAng = atan2(PMaskRegProp->PY[I]-Y,PMaskRegProp->PX[I]-X);
    else
      NAng = atan2(PMaskRegProp->PY[0]-Y,PMaskRegProp->PX[0]-X);

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

extern void MaskOutRegion(MaskRegRec *PMaskRegProp, int *ifs,
			  int dimx, int dimy)
{
  int LX,HX,LY,HY,X,Y,I;

  for (I=0;I<PMaskRegProp->NPoint;I++) {
    if ((I==0)||(PMaskRegProp->PX[I]<LX))
      LX = PMaskRegProp->PX[I]-1;
    if ((I==0)||(PMaskRegProp->PX[I]>HX))
      HX = PMaskRegProp->PX[I]+1;
    if ((I==0)||(PMaskRegProp->PY[I]<LY))
      LY = PMaskRegProp->PY[I]-1;
    if ((I==0)||(PMaskRegProp->PY[I]>HY))
      HY = PMaskRegProp->PY[I]+1;
  }
  for (X=LX;X<HX;X++)
    for (Y=LY;Y<HY;Y++) {
      if ((X>=0)&&(Y>=0)&&(X<dimx)&&(Y<dimy)&&
	  WithinPolygon(PMaskRegProp,(float)X,(float)Y))
	ifs[X+Y*dimx] = 0;
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

extern void ModifyMask(ImageMaskRec *PImageMaskProp, char *Path)
{
  fitsfile *fptr;
  char NFN[200];
  int status,nval,anynul,hdutype;
  long fpixel,naxes[2],nelements;
  char comment1[500],comment2[500];
  int *nfs,I;

  status = 0;
  sprintf(NFN,"%s%s",Path,PImageMaskProp->Name);

  if (!DoesFileExist(NFN)) {
    fprintf(stderr,"Skipping masks to %s...\n",NFN);
    return;
  }
  fprintf(stderr,"Applying masks to %s...\n",NFN);

  fits_open_file(&fptr, NFN, READWRITE, &status);

  if (fits_movabs_hdu(fptr,1,&hdutype,&status))
    printerror( status );

  nval = 0;
  fpixel = 1;

  if (fits_read_key(fptr,TINT,"NAXIS1",&naxes[0],comment1,&status))
    printerror( status );
  if (fits_read_key(fptr,TINT,"NAXIS2",&naxes[1],comment2,&status))
    printerror( status );
  nelements = naxes[0]*naxes[1];
  nfs = (int *)calloc(naxes[0]*naxes[1],sizeof(int));
  if (fits_read_img(fptr,TINT,fpixel,nelements,&nval,nfs,
		    &anynul,&status))
    printerror( status );

  for (I=0;I<PImageMaskProp->NRegion;I++)
    MaskOutRegion(&PImageMaskProp->PMaskRegProp[I],nfs,naxes[0],naxes[1]);

  if (fits_write_img(fptr,TINT,fpixel,nelements,nfs,&status))
    printerror( status );

  if (fits_close_file(fptr,&status))
    printerror( status );  

  free(nfs);
}

extern void ModifyMasks(MaskRec *PMaskProp, char *Path)
{
  int I;

  for (I=0;I<PMaskProp->NImage;I++) {
    ModifyMask(&PMaskProp->PImageMaskProp[I],Path);
  }
}

main(int argc, char *argv[])
{
  MaskRec MaskProp;
  char NStr[100];

  fprintf(stderr,"Applying masks...\n");
  ReadFile(argv[1],&MaskProp);
  if (argv[2][strlen(argv[2])-1] != '/')
    sprintf(NStr,"%s/",argv[2]);
  else
    sprintf(NStr,"%s",argv[2]);

  ModifyMasks(&MaskProp,NStr);
}
