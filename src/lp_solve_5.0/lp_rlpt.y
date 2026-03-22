/* ========================================================================= */
/* NAME  : lp_rlpt.y                                                         */
/* ========================================================================= */


%token VAR CONS INTCONS VARIABLECOLON INF FR SEC_INT SEC_SEC SEC_SOS SOSTYPE SIGN RE_OPLE RE_OPGE MINIMISE MAXIMISE SUBJECTTO BOUNDS END UNDEFINED


%{
#include <string.h>
#include <ctype.h>

#include "lpkit.h"
#include "yacc_read.h"

static char Last_var[NAMELEN], Last_var0[NAMELEN];
static REAL f, f0;
static int x;
static int Sign;
static int isign, isign0;      /* internal_sign variable to make sure nothing goes wrong */
		/* with lookahead */
static int make_neg;   /* is true after the relational operator is seen in order */
		/* to remember if lin_term stands before or after re_op */
static int Within_gen_decl = FALSE; /* TRUE when we are within an gen declaration */
static int Within_bin_decl = FALSE; /* TRUE when we are within an bin declaration */
static int Within_sec_decl = FALSE; /* TRUE when we are within an sec declaration */
static int Within_sos_decl = FALSE; /* TRUE when we are within an sos declaration */
static short SOStype; /* SOS type */
static int SOSNr;
static int weight; /* SOS weight */
static int SOSweight = 0; /* SOS weight */

static int HadConstraint;
static int HadVar;
static int Had_lineair_sum;

#define YY_FATAL_ERROR lex_fatal_error

/* let's please C++ users */
#ifdef __cplusplus
extern "C" {
#endif

static int wrap(void)
{
  return(1);
}

#ifdef __cplusplus
};
#endif

#define yywrap wrap
#define yyerror read_error

#include "lp_rlpt.h"

%}



%start inputfile
%%

EMPTY: /* EMPTY */
                ;

inputfile	:
{
  isign = 0;
  make_neg = 0;
  Sign = 0;
  HadConstraint = FALSE;
  HadVar = FALSE;
}
		  objective_function
		  x_constraints
		  bounds
                  int_sec_sos_declarations
                  end
		;

/* start objective_function */

/*

 objective_function: MAXIMISE of | MINIMISE of;
 of:                 real_of | VARIABLECOLON real_of;
 real_of:            lineair_sum;

*/

objective_function:
                  MAXIMISE
                  of
{
  set_obj_dir(TRUE);
}
                | MINIMISE
                  of
{
  set_obj_dir(FALSE);
}
                ;

of              : real_of
                | VARIABLECOLON
{
  if(!add_constraint_name(Last_var))
    YYABORT;
  /* HadConstraint = TRUE; */
}
                  real_of
                ;

real_of:          lineair_sum
{
  add_row();
  /* HadConstraint = FALSE; */
  HadVar = FALSE;
  isign = 0;
  make_neg = 0;
}
		;

/* end objective_function */



/* start constraints */

/*

 constraints:        EMPTY | x_constraints;
 x_constraints:      SUBJECTTO constraints;
 constraints:        constraint | constraints constraint;
 constraint:         real_constraint | VARIABLECOLON real_constraint;
 real_constraint:    lineair_sum RE_OP cons_term RHS_STORE;
 RE_OP:              RE_OPLE | RE_OPGE;
 cons_term:          x_SIGN REALCONS | INF;
 lineair_sum:        lineair_term | lineair_sum lineair_term;
 lineair_term:       x_SIGN VARIABLE VAR_STORE | x_SIGN REALCONS VARIABLE VAR_STORE;
 x_SIGN:             EMPTY | SIGN;
 VARIABLE:           VAR | FR;
 REALCONS:           INTCONS | CONS;

*/

x_constraints:    SUBJECTTO
{
  HadConstraint = TRUE;
}
                  constraints
{
  HadConstraint = FALSE;
}
                ;

constraints	: constraint
		| constraints
		  constraint
		;

constraint      : real_constraint
                | VARIABLECOLON
{
  if(!add_constraint_name(Last_var))
    YYABORT;
  /* HadConstraint = TRUE; */
}
                  real_constraint
                ;

real_constraint	: lineair_sum
		  RE_OP
{
  if(!store_re_op((char *) yytext, HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  make_neg = 1;
}
		  cons_term
                  RHS_STORE
{
  Had_lineair_sum = TRUE;
  add_row();
  /* HadConstraint = FALSE; */
  HadVar = FALSE;
  isign = 0;
  make_neg = 0;
  null_tmp_store(TRUE);
}
		;

lineair_sum	: lineair_term
		| lineair_sum
		  lineair_term
		;

lineair_term	: x_SIGN
                  VARIABLE
{
  f = 1.0;
}
                  VAR_STORE
		| x_SIGN
                  REALCONS
		  VARIABLE
                  VAR_STORE
                ;

RE_OP: RE_OPLE | RE_OPGE
                ;

cons_term:        x_SIGN
		  REALCONS
                | INF
{
  isign = Sign;
}
                ;

/* end constraints */



/* start bounds */

/*

 bounds:             EMPTY | BOUNDS x_bounds;
 x_bounds:           bound | x_bounds bound;
 bound:              VARIABLE VAR_STORE bound2 | cons_term RE_OP VARIABLE VAR_STORE RHS_STORE optionalbound;
 VARIABLE:           VAR | FR;
 bound2:             RE_OP cons_term RHS_STORE | FR RHS_STORE RHS_STORE;
 cons_term:          x_SIGN REALCONS | INF;
 optionalbound:      EMPTY | RE_OP cons_term RHS_STORE;

*/

bounds:           EMPTY
                | BOUNDS
                  x_bounds
                ;

x_bounds:         bound
                | x_bounds
                  bound
                ;

bound:            VARIABLE
{
  f = 1.0;
  isign = 0;
}
                  VAR_STORE
                  bound2
		| cons_term
{
  f0 = f;
  isign0 = isign;
}
		  RE_OP
{
  if(!store_re_op((char *) yytext, HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  make_neg = 0;
}
                  VARIABLE
{
  isign = 0;
  f = -1.0;
}
                  VAR_STORE
{
  isign = isign0;
  f = f0;
}
                  RHS_STORE
{
  if(!store_bounds(TRUE))
    YYABORT;
}
                  optionalbound
{
  /* HadConstraint = FALSE; */
  HadVar = FALSE;
  isign = 0;
  make_neg = 0;
  null_tmp_store(TRUE);
}
                ;

bound2:           RE_OP
{
  if(!store_re_op((char *) yytext, HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  make_neg = 1;
}
		  cons_term
                  RHS_STORE
{
  if(!store_bounds(TRUE))
    YYABORT;
  /* HadConstraint = FALSE; */
  HadVar = FALSE;
  isign = 0;
  make_neg = 0;
  null_tmp_store(TRUE);
}
                | FR
{
  if(!store_re_op(">", HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  make_neg = 1;
  isign = 0;
  f = -DEF_INFINITE;
}
                  RHS_STORE
{
  if(!store_bounds(FALSE))
    YYABORT;

  if(!store_re_op("<", HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  f = DEF_INFINITE;
  isign = 0;
}
                  RHS_STORE
{
  if(!store_bounds(FALSE))
    YYABORT;
  /* HadConstraint = FALSE; */
  HadVar = FALSE;
  isign = 0;
  make_neg = 0;
  null_tmp_store(TRUE);
}
                ;


optionalbound:    EMPTY
                | RE_OP
{
  if(!store_re_op((*yytext == '<') ? ">" : (*yytext == '>') ? "<" : (char *) yytext, HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  make_neg = 0;
  isign = 0;
}
                  cons_term
                  RHS_STORE
{
  if(!store_bounds(TRUE))
    YYABORT;
}
                ;

/* end bounds */



REALCONS: INTCONS | CONS
                ;

RHS_STORE:        EMPTY
{
  if (    (isign || !make_neg)
      && !(isign && !make_neg)) /* but not both! */
    f = -f;
  if(!rhs_store(f, HadConstraint, HadVar, Had_lineair_sum))
    YYABORT;
  isign = 0;
}
                ;

x_SIGN:           EMPTY
{
  isign = 0;
}
                | SIGN
{
  isign = Sign;
}
                ;

VAR_STORE:        EMPTY
{
  if (    (isign || make_neg)
      && !(isign && make_neg)) /* but not both! */
    f = -f;
  if(!var_store(Last_var, f, HadConstraint, HadVar, Had_lineair_sum)) {
    yyerror("var_store failed");
    YYABORT;
  }
  /* HadConstraint |= HadVar; */
  HadVar = TRUE;
  isign = 0;
}
                ;


/* start int_sec_sos_declarations */

int_sec_sos_declarations: x_int_declarations
                          x_sec_declarations
                          x_sos_declarations
                ;

VARIABLES:        ONEVARIABLE
                | VARIABLES
                  ONEVARIABLE
                ;

SOSVARIABLES:     ONESOSVARIABLE
                | SOSVARIABLES
                  ONESOSVARIABLE
                ;

ONEVARIABLE:      VARIABLE
{
  storevarandweight(Last_var);
}
                ;

ONESOSVARIABLE:   VARIABLE
{
  SOSNr++;
  weight = SOSNr;
  storevarandweight(Last_var);
  set_sos_weight(weight, 2);
}
                | VARIABLECOLON
{
  storevarandweight(Last_var);
}
                  INTCONS
{
  weight = (int) (f + .1);
  set_sos_weight(weight, 2);
}
                ;

x_int_declarations:
                  EMPTY
                | int_declarations
                ;

int_declarations: int_declaration
                | int_declarations
                  int_declaration
                ;

int_declaration: SEC_INT
{
  check_int_sec_sos_decl(Within_gen_decl ? 1 : Within_bin_decl ? 2 : 0, 0, 0);
}
                 VARIABLES
                ;

x_sec_declarations:
                 EMPTY
                | sec_declarations
                ;

sec_declarations: sec_declaration
                | sec_declarations
                  sec_declaration
                ;

sec_declaration: SEC_SEC
{
  check_int_sec_sos_decl(0, 1, 0);
}
                 VARIABLES
                ;


x_sos_declarations:
                  EMPTY
                | sos_declarations
                ;

sos_declarations: sos_declaration
                | sos_declarations
                  sos_declaration
                ;

sos_declaration: SEC_SOS
                 x_single_sos_declarations
                ;

x_single_sos_declarations:
                  EMPTY
                | single_sos_declarations
                ;

single_sos_declarations:
                  single_sos_declaration
                | single_sos_declarations
                  single_sos_declaration
                ;

single_sos_declaration:
                 SOSTYPE
{
  char buf[16];

  check_int_sec_sos_decl(0, 0, 1);
  SOSweight++;
  sprintf(buf, "SOS%d", SOSweight);
  storevarandweight(buf);
  SOStype = Last_var[1] - '0';
  set_sos_type(SOStype);
  check_int_sec_sos_decl(0, 0, 2);
  weight = 0;
  SOSNr = 0;
}
                 SOSVARIABLES
{
  set_sos_weight(SOSweight, 1);
}
                ;

/* end int_sec_sos_declarations */


VARIABLE: VAR | FR
                ;

/* start end */

end:              EMPTY
                | END
                ;

/* end end */

%%

static void yy_delete_allocated_memory(void)
{
  /* free memory allocated by flex. Otherwise some memory is not freed.
     This is a bit tricky. There is not much documentation about this, but a lot of
     reports of memory that keeps allocated */

  /* If you get errors on this function call, just comment it. This will only result
     in some memory that is not being freed. */

# if defined YY_CURRENT_BUFFER
    /* flex defines the macro YY_CURRENT_BUFFER, so you should only get here if lp_rlpt.h is
       generated by flex */
    /* lex doesn't define this macro and thus should not come here, but lex doesn't has
       this memory leak also ...*/

    yy_delete_buffer(YY_CURRENT_BUFFER); /* comment this line if you have problems with it */
    yy_init = 1; /* make sure that the next time memory is allocated again */
    yy_start = 0;
# endif
}

static int parse(void)
{
  return(yyparse());
}

lprec * __WINAPI read_lpt(FILE *filename, int verbose, char *lp_name)
{
  yyin = filename;
  return(yacc_read(verbose, lp_name, &yylineno, parse, yy_delete_allocated_memory));
}

lprec * __WINAPI read_LPT(char *filename, int verbose, char *lp_name)
{
  FILE *fpin;
  lprec *lp = NULL;

  if((fpin = fopen(filename, "r")) != NULL) {
    lp = read_lpt(fpin, verbose, lp_name);
    fclose(fpin);
  }
  return(lp);
}
