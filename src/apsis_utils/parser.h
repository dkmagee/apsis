typedef double (*BrParmRec) (char *, int, char *, int *); 

BrParmRec BrParameter[10];

double Parameter[10]; /* Parameters that are to be passed to function */

char PARSEERROR[300];
int WhParseEntry;

#define MAXCHAR 30
#define MAXFUNC 50

typedef char MaxCompStr[1000];

typedef struct {
  int ReplFunc;
  char NewName[20];
  int NArg;
  int UseMap[10];
  int Map[10];
  MaxCompStr Str[10];
} FunctionReplRec;

FunctionReplRec FunctionReplProp;

typedef struct {
  struct ListVarType *prev;
  int Numchar;
  int NumFunc;
  char charList[MAXCHAR][10];
  double ValList[MAXCHAR];
  char FuncList[MAXFUNC][10];
  double (*FuncProc[MAXFUNC]) (int, int);
  int FuncNArg[MAXFUNC];
  int IDArg[MAXFUNC];
  int p;
  char Ch;
  MaxCompStr formula;
  int ProdNew;
  MaxCompStr formulaN;
  int Breakpoint;
  int ErrorType;
} ListVarType;

ListVarType *PCurListVar;

#ifndef true
# define true    1
# define false   0
#endif

extern double ComputeFormula (int *p2, char *Strg, int *Err,
			      int Handle, int ProduceNew);

extern void SetUpVar (char *Ch);

extern void SetVar (char *Ch, double Value);

extern void SetUpFunc(char *Ch, double (*Loc) (int, int), int IDArg, int NArg);

extern void InitVar();

extern void DestVar();

extern char *strsub(char *STR1, char *Str, int Index, int Length);

extern void Parser_init();

/* End. */
