#include "parser.h"
#include <math.h>
#include <stdio.h>
#ifdef USEMPATROL
#include <mpatrol.h>
#else
#include <malloc.h>
#endif
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

static double Expr ();

extern char *strsub(STR1,Str,Index,Length)
char *STR1;
char *Str;
int Index;
int Length;
{
  int I;
  for (I=0; I<Length; I++) 
    STR1[I] = Str[Index-1+I];
  STR1[Length] = '\0';
  return STR1;
}

static void NextP(int TotalC)
{  /*Next P*/
  int l,i;
  
  if (PCurListVar->ProdNew) {
    l = strlen(PCurListVar->formulaN);
    if (PCurListVar->p)
      PCurListVar->formulaN[l] = PCurListVar->Ch;
    for (i=0;i<TotalC-1;i++)
      PCurListVar->formulaN[l+i+1] = PCurListVar->formula[PCurListVar->p+i];
    PCurListVar->formulaN[l+i+1] = 0;
  }
  PCurListVar->p += TotalC-1;
  do {
    PCurListVar->p++;
    if (PCurListVar->p <= strlen(PCurListVar->formula))
      PCurListVar->Ch = PCurListVar->formula[PCurListVar->p - 1];
    else
      PCurListVar->Ch = '\015';
  } while (PCurListVar->Ch == ' ');
}  /*NextP*/

static void CopyP(int TotalC, char *NStr)
{ 
  int l,i;
  
  if (PCurListVar->ProdNew) {
    int tc;

    l = strlen(PCurListVar->formulaN);
    tc = strlen(NStr);
    for (i=0;i<tc;i++)
      PCurListVar->formulaN[l+i] = NStr[i];
    PCurListVar->formulaN[l+i] = 0;
  }
  PCurListVar->p += TotalC-1;
  do {
    PCurListVar->p++;
    if (PCurListVar->p <= strlen(PCurListVar->formula))
      PCurListVar->Ch = PCurListVar->formula[PCurListVar->p - 1];
    else
      PCurListVar->Ch = '\015';
  } while (PCurListVar->Ch == ' ');
}  /*NextP*/

static void ProcessAsNumber(f)
double *f;
{
  int code, start;
  char STR1[256];

  start = PCurListVar->p;
  do {
    NextP(1);
  } while (isdigit(PCurListVar->Ch)||
	  (PCurListVar->Ch=='.'));
  if (PCurListVar->Ch == '.') {
    do {
      NextP(1);
    } while (isdigit(PCurListVar->Ch)||
		(PCurListVar->Ch=='.'));
  }
  if (PCurListVar->Ch == 'E') {
    NextP(1);
    do {
      NextP(1);
    } while (isdigit(PCurListVar->Ch)||
		(PCurListVar->Ch=='.'));
  }
  code = (sscanf(strsub(STR1, PCurListVar->formula,
	  start, PCurListVar->p - start),
	  "%lg", f) == 0);
}  /*ProcessAsNumber*/

static void ProcessAsNewExpr(f)
double *f;
{
  NextP(1);
  *f = Expr();
  if (PCurListVar->ErrorType)
    return;
  if (PCurListVar->Ch == ')')
    NextP(1);
  else {
    PCurListVar->ErrorType = 1;
  }
}  /*ProcessAsNewExpr*/

static double Fact(i)
int i;
{
  if (i > 0)
    return (i * Fact(i - 1));
  else
    return 1.0;
}  /*Fact*/


static double Fct ();

static double S_Fact();

void FunctionReplace(FunctionReplRec FunctionReplProp, 
		     int *PpBeg, int *pLoc)
{
  int i,j;

  for (i=0;i<strlen(FunctionReplProp.NewName);i++) {
    PCurListVar->formulaN[*PpBeg] = FunctionReplProp.NewName[i];
    (*PpBeg)++;
  }
  PCurListVar->formulaN[*PpBeg] = '(';
  (*PpBeg)++;
  for (i=0;i<FunctionReplProp.NArg;i++) {
    if (FunctionReplProp.UseMap[i]==1) {
      int Wh;

      Wh = FunctionReplProp.Map[i];
      for (j=pLoc[Wh];j<pLoc[Wh+1]-1;j++) {
	PCurListVar->formulaN[*PpBeg] = PCurListVar->formula[j-1];
	(*PpBeg)++;
      }	   	    
    }
    else if (FunctionReplProp.UseMap[i]==2) {
      int Wh;

      Wh = FunctionReplProp.Map[i];
      for (j=pLoc[Wh];(j<pLoc[Wh+1]-1)&&(PCurListVar->formula[j-1]!='-');j++) {
	PCurListVar->formulaN[*PpBeg] = PCurListVar->formula[j-1];
	(*PpBeg)++;
      }
      PCurListVar->formulaN[*PpBeg] = '-';
      (*PpBeg)++;
      for (j=0;j<strlen(FunctionReplProp.Str[i]);j++) {
	PCurListVar->formulaN[*PpBeg] = FunctionReplProp.Str[i][j];
	(*PpBeg)++;
      }	   
    }
    else {
      for (j=0;j<strlen(FunctionReplProp.Str[i]);j++) {
	PCurListVar->formulaN[*PpBeg] = FunctionReplProp.Str[i][j];
	(*PpBeg)++;
      }	   
    }
    if (i==FunctionReplProp.NArg-1) {
      PCurListVar->formulaN[*PpBeg] = ')';
      (*PpBeg)++;	    	    
    }
    else {	
      PCurListVar->formulaN[*PpBeg] = ',';
      (*PpBeg)++;	    
    }
  }	
  PCurListVar->formulaN[*PpBeg] = 0;
}


static void ProcessAsStandardFunction(double *f)
{
  int index,i;
  char STR1[256];
  int FORLIM;
  double Parm[10];
  int pLoc[10];
  int pBeg;
  BrParmRec TBrParameter[10];
  FunctionReplRec TFunctionReplProp;

  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "ABS")) {
    NextP(3);
    *f = S_Fact();
    *f = fabs(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 4), "SQRT")) {
    NextP(4);
    *f = S_Fact();
    *f = sqrt(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "SQR")) {
    NextP(3);
    *f = S_Fact();
    *f *= *f;
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "SIN")) {
    NextP(3);
    *f = S_Fact();
    *f = sin(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "COS")) {
    NextP(3);
    *f = S_Fact();
    *f = cos(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 6), "ARCTAN")) {
    NextP(6);
    *f = S_Fact();
    *f = atan(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "LN")) {
    NextP(2);
    *f = S_Fact();
    *f = log(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "LOG")) {
    NextP(3);
    *f = S_Fact();
    *f = log(*f) / log(10.0);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 3), "EXP")) {
    NextP(3);
    *f = S_Fact();
    *f = exp(*f);
    return;
  }
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 4), "FACT")) {
    NextP(4);
    *f = S_Fact();
    *f = Fact((int)((long)(*f)));
    return;
  }
  PCurListVar->ErrorType = 0;
  FORLIM = PCurListVar->Numchar;
  for (index = 0; index < FORLIM; index++) {
    int L,L2;

    L = strlen(PCurListVar->charList[index]);
    L2 = strlen(PCurListVar->formula);
    if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p,L),
		PCurListVar->charList[index])) {
      if ((L+PCurListVar->p-1 >= L2) || 
	  (!isalnum(PCurListVar->formula[PCurListVar->p+L-1]))) {
	*f = PCurListVar->ValList[index];
	NextP(strlen(PCurListVar->charList[index]));
	return;
      }
    }
  }
  FORLIM = PCurListVar->NumFunc;
  for (index = 0; index < FORLIM; index++) {
    if ((!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p,
			strlen(PCurListVar->FuncList[index])),
		 PCurListVar->FuncList[index])) &&
        (PCurListVar->formula[PCurListVar->p + 
			     strlen(PCurListVar->FuncList[index]) - 1]
	 == '(')) {
      pBeg = strlen(PCurListVar->formulaN);
      PCurListVar->ErrorType = 0;
      NextP(strlen(PCurListVar->FuncList[index]));
      if (PCurListVar->Ch != '(') {
	PCurListVar->ErrorType = 3;
        return;
      }
      NextP(1);
      TFunctionReplProp.ReplFunc = 0;
      if (PCurListVar->ProdNew) {
	FunctionReplProp.ReplFunc = 0;
	for (i = 0; i < PCurListVar->FuncNArg[index]; i++)
	  BrParameter[i] = NULL;
	*f = (*(PCurListVar->FuncProc)[index])(0,PCurListVar->IDArg[index]);
	for (i = 0; i < PCurListVar->FuncNArg[index]; i++)
	  TBrParameter[i] = BrParameter[i];
	TFunctionReplProp = FunctionReplProp;
      }      
      PARSEERROR[0] = 0;
      for (i = 0; i < PCurListVar->FuncNArg[index]; i++) {
	if (i > 0) {
	  if (PCurListVar->Ch != ',') {
	    PCurListVar->ErrorType = 4;
            return;
	  }
	  NextP(1);
	}
	pLoc[i] = PCurListVar->p;
	if ((PCurListVar->ProdNew)&&(TBrParameter[i])) {
	  int Adv;
	  char NStr[100];

	  Parm[i] = (*(TBrParameter[i]))(PCurListVar->formula,PCurListVar->p,
					 NStr,&Adv);
	  if (Adv) CopyP(Adv,NStr);
	  else {
	    Parm[i] = S_Fact();
	    if (PCurListVar->ErrorType)
	      return;
	  }
	}
	else {
	  Parm[i] = S_Fact();
	  if (PCurListVar->ErrorType)
	    return;
	}
      }
      if (PCurListVar->Ch != ')') {
	PCurListVar->ErrorType = 5;
        return;
      }
      NextP(1);
      pLoc[i] = PCurListVar->p;
      for (i = 0; i < PCurListVar->FuncNArg[index]; i++)
	Parameter[i] = Parm[i];
      *f = (*(PCurListVar->FuncProc)[index])(1,PCurListVar->IDArg[index]);
      if (TFunctionReplProp.ReplFunc) {	       
	FunctionReplace(TFunctionReplProp,&pBeg,pLoc);
      }
      if (PARSEERROR[0] != 0) {
	PCurListVar->ErrorType = 6;
	PCurListVar->p = pLoc[WhParseEntry];
        return;
      }
      return;
    }
  }
  PCurListVar->ErrorType = 2;
}  /*ProcessAsStandardFunction*/


static double Fct()
{
  double f;

  if (isdigit(PCurListVar->Ch)||(PCurListVar->Ch=='.')) {
    ProcessAsNumber(&f);
    return f;
  }
  if (PCurListVar->Ch == '(')
    ProcessAsNewExpr(&f);
  else
    ProcessAsStandardFunction(&f);
  return f;
}  /*Fct*/


static double S_Fact()
{
  if (PCurListVar->Ch == '-') {
    NextP(1);
    return (-Fct());
  } else if (PCurListVar->Ch == '+') {
    NextP(1);
    return (Fct());
  } else
    return (Fct());
}  /*S_Fact*/


static double Term()
{
  double t, t1;

  t = S_Fact();
  if (PCurListVar->ErrorType) return 1.0;
  while (PCurListVar->Ch == '^') {
    NextP(1);
    t1 = S_Fact();
    if (t < -1e-8 && (((long)floor(t1 + 0.5)) & 1) == 1) {
      t = -exp(log(fabs(t)) * t1);
      continue;
    }
    if (fabs(t) != 0)
      t = exp(log(fabs(t)) * t1);
    else if (t1 != 0)
      t = 0.0;
    else
      t = 1.0;
  }
  return t;
}  /*Term*/


static double SmplExpr()
{
  double s;
  char Operator;

  s = Term();
  if (PCurListVar->ErrorType) return 1.0;
  while (PCurListVar->Ch == '/' || PCurListVar->Ch == '*') {
    Operator = PCurListVar->Ch;
    NextP(1);
    switch (Operator) {

    case '*':
      s *= Term();
      break;

    case '/':
      s /= Term();
      break;
    }
  }
  return s;
}


static double NSExpr()
{
  double ns,tv;
  char Operator;

  ns = SmplExpr();
  if (PCurListVar->ErrorType) return 1.0;
  while (PCurListVar->Ch == '-' || PCurListVar->Ch == '+') {
    Operator = PCurListVar->Ch;
    NextP(1);
    switch (Operator) {

    case '+':
      tv = SmplExpr();
      ns += tv;
      break;

    case '-':
      tv = SmplExpr();
      ns -= tv;
      break;
    }
  }
  return ns;
}


static double B3Expr()
{
  double BE3;
  char STR1[256];

  BE3 = NSExpr();
  if (PCurListVar->ErrorType) return 1.0;
  while (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "&&")) {
    NextP(2);
    BE3 = ((long)floor(BE3 + 0.5)) & ((long)floor(NSExpr() + 0.5));
  }
  return BE3;
}


static double B2Expr()
{
  double BE2;
  char STR1[256];

  BE2 = B3Expr();
  if (PCurListVar->ErrorType) return 1.0;
  while (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "||")) {
    NextP(2);
    BE2 = ((long)floor(BE2 + 0.5)) | ((long)floor(B3Expr() + 0.5));
  }
  return BE2;
}


static double Expr()
{
  double e;
  char STR1[256];

  e = B2Expr();
  if (PCurListVar->ErrorType) return 1.0;
  if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "==")) {
    NextP(2);
    if (e == B2Expr())
      return 1.0;
    else
      return 0.0;
  } else if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "!=")) {
    NextP(2);
    if (e != B2Expr())
      return 1.0;
    else
      return 0.0;
  } else if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), "<=")) {
    NextP(2);
    if (e <= B2Expr())
      return 1.0;
    else
      return 0.0;
  } else if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 2), ">=")) {
    NextP(2);
    if (e >= B2Expr())
      return 1.0;
    else
      return 0.0;
  } else if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 1), "<")) {
    NextP(1);
    if (e < B2Expr())
      return 1.0;
    else
      return 0.0;
  } else if (!strcmp(strsub(STR1, PCurListVar->formula, PCurListVar->p, 1), ">")) {
    NextP(1);
    if (e > B2Expr())
      return 1.0;
    else
      return 0.0;
  } else
    return e;
}

static void Eval(char *Formula_, double *Value, int *Error, int Handle, 
		 int ProduceNew)
{
  int FORLIM;
  char STR1[256];

  if (strlen(Formula_) > sizeof(MaxCompStr)) {
    fprintf(stderr,"Eval inputstring too long.\n");
    exit(1);
  }
  strcpy(PCurListVar->formula, Formula_);
  PCurListVar->ProdNew = ProduceNew;
  FORLIM = strlen(PCurListVar->formula);
  /*  for (i = 0; i < FORLIM; i++)
      PCurListVar->formula[i] = toupper(PCurListVar->formula[i]);*/
  if (PCurListVar->formula[0] == '.')
    sprintf(PCurListVar->formula, "0%s", strcpy(STR1, PCurListVar->formula));
  PCurListVar->p = 0;
  PCurListVar->ErrorType = 0;
  NextP(1);
  *Value = Expr();

  if (!PCurListVar->ErrorType)
    *Error = false;
  else
    *Error = true;
  if (*Error) {
    int i,j;    
    if (Handle==0)
      fprintf(stderr,"Error in %s at position %i\n",
	      PCurListVar->formula,PCurListVar->p);

    switch (PCurListVar->ErrorType) {
    case 1:
      sprintf(PARSEERROR,"')' expected.");
      break;
    case 2:
      sprintf(PARSEERROR,"Expression not found.");
      break;
    case 3:
      sprintf(PARSEERROR,"'(' expected.");
      break;
    case 4:
      sprintf(PARSEERROR,"',' expected.");
      break;
    case 5:
      sprintf(PARSEERROR,"')' expected.  Expression has too many parameters");
      break;
    }
    if (Handle == 0) {
      j = 0;
      i = i / j;
      exit(1);
    }
  }
  PCurListVar->Breakpoint = PCurListVar->p;
}  /*Eval*/


extern double ComputeFormula(int *p2, char *Strg, int *Err, int Handle, 
			     int ProduceNew)
{
  double r;

  Eval(Strg, &r, Err,Handle,ProduceNew);
  *p2 = PCurListVar->p;
  return r;
}  /*ComputeFormula*/

void SetUpVar(Ch)
char *Ch;
{
  if (PCurListVar->Numchar < MAXCHAR) {
    PCurListVar->Numchar++;
    strcpy(PCurListVar->charList[PCurListVar->Numchar - 1], Ch);
  }
}


void SetVar(char *Ch, double Value)
{
  int i, FORLIM;

  FORLIM = PCurListVar->Numchar;
  for (i = 0; i < FORLIM; i++) {
    if (!strcmp(PCurListVar->charList[i], Ch))
      PCurListVar->ValList[i] = Value;
  }
}

extern void SetUpFunc(char *Ch, double (*Loc) (int, int), int IDArg, int NArg)
{
  if (PCurListVar->NumFunc < MAXFUNC) {
    PCurListVar->NumFunc++;
    strcpy(PCurListVar->FuncList[PCurListVar->NumFunc - 1], Ch);
	PCurListVar->FuncProc[PCurListVar->NumFunc-1] = Loc;
    PCurListVar->FuncNArg[PCurListVar->NumFunc-1] = NArg;
    PCurListVar->IDArg[PCurListVar->NumFunc-1] = IDArg;
  }
}

extern void InitVar()
{
  ListVarType *PPrevListVar;

  PPrevListVar = PCurListVar;
  PCurListVar = (ListVarType *)calloc(1,sizeof(ListVarType));
  PCurListVar->prev = (struct ListVarType *)PPrevListVar;
  PCurListVar->Numchar = 0;
  PCurListVar->NumFunc = 0;
}

extern void DestVar()
{
  ListVarType *PPrevListVar;

  PPrevListVar = (ListVarType *)PCurListVar->prev;
  free(PCurListVar);
  PCurListVar = PPrevListVar;
}

extern void Parser_init()
{
  PCurListVar = NULL;
}
