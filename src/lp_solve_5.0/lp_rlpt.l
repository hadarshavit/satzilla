WS   [ \t\r\n]+
LT   [A-Za-z!\"#$%&\(\)/,;\?@_`'\{\}|~]
KR   [A-Za-z0-9!\"#$%&\(\)/,\.;\?@_`'\{\}|~]
DI   [0-9]
INNM {DI}+
NM   {DI}*\.?{DI}+([Ee][-+]?{DI}+)?
VR   {LT}{KR}*
S_OP [-+]+
OPLE [<]?=?
OPGE [>]=?
INF  [ \r\n\t]*[-+][Ii][Nn][Ff]([Ii][Nn][Ii][Tt]([Ee]|[Yy]))?
FR   [Ff][Rr][Ee][Ee]
MIN  ^[mM][iI][nN](([iI][mM][iI][zZ][eE])?|([iI][mM][uU][mM])?)
MAX  ^[mM][aA][xX](([iI][mM][iI][zZ][eE])?|([iI][mM][uU][mM])?)
ST   ^([sS][uU][bB][jJ][eE][cC][tT][ ][tT][oO])|([Ss][Uu][Cc][Hh][ ][Tt][Hh][Aa][Tt])|([sS][tT])|([sS][.][tT][.])|([sS][tT][.])
BD   ^[bB][oO][uU][nN][dD][sS]?
INT  ^([Ii][Nn][Tt]([Ee][Gg][Ee][Rr][Ss]?)?)|([Gg][Ee][Nn]([Ee][Rr][Aa][Ll][Ss]?)?)|([Bb][Ii][Nn]([Aa][Rr]([Yy]|([Ii][Ee][Ss])))?)
SEC  ^([Ss][Ee][Mm][Ii]-[Cc][Oo][Nn][Tt][Ii][Nn][Uu][Oo][Uu][Ss])|([Ss][Ee][Mm][Ii])|([Ss][Ee][Mm][Ii][Ss])
SOS  ^[Ss][Oo][Ss]
SOST [Ss][12]::
LBL  {VR}:
END  ^[Ee][Nn][Dd]

%start LINECOMMENT

%%

<INITIAL>"\\" {
  BEGIN LINECOMMENT;
} /* begin skip LINECOMMENT */

<LINECOMMENT>\n {
  BEGIN INITIAL;
} /* end skip LINECOMMENT */

<LINECOMMENT>\r {
  BEGIN INITIAL;
} /* end skip LINECOMMENT */

<LINECOMMENT>. {
}

<INITIAL>{WS} {
}

<INITIAL>{MIN} {
  /* strcpy(Last_var, (char *)yytext); */
  return(MINIMISE);
}

<INITIAL>{MAX} {
  /* strcpy(Last_var, (char *)yytext); */
  return(MAXIMISE);
}

<INITIAL>{ST} {
  /* strcpy(Last_var, (char *)yytext); */
  return(SUBJECTTO);
}

<INITIAL>{BD} {
  /* strcpy(Last_var, (char *)yytext); */
  return(BOUNDS);
}

<INITIAL>{END} {
  /* strcpy(Last_var, (char *)yytext); */
  return(END);
}

<INITIAL>{INNM} {
  f = atof((char *)yytext);
  return(INTCONS);
} /* f contains the last float */

<INITIAL>{NM} {
  f = atof((char *)yytext);
  return(CONS);
} /* f contains the last float */

<INITIAL>{INF} {
  char *ptr;

  f = DEF_INFINITE;
  Sign = 0;
  ptr = (char *)yytext;
  while (isspace(*ptr)) ptr++;
  if(*ptr == '-')
    Sign = 1;
  return(INF);
} /* f contains the last float */

<INITIAL>{FR} {
  strcpy(Last_var0, Last_var);
  strcpy(Last_var, (char *)yytext);
  return(FR);
}

<INITIAL>{S_OP} {
  Sign = 0;
  for(x = 0; x < yyleng; x++)
    if(yytext[x] == '-' || yytext[x] == '+')
      Sign = (Sign == (yytext[x] == '+'));
  return (SIGN);
  /* Sign is TRUE if the sign-string
     represents a '-'. Otherwise Sign
     is FALSE */
}

<INITIAL>{INT} {
  Within_gen_decl = Within_bin_decl = Within_sec_decl = Within_sos_decl = FALSE;
  if((toupper(*yytext) == 'G') || (toupper(*yytext) == 'I'))
    Within_gen_decl = TRUE;
  else
    Within_bin_decl = TRUE;
  return(SEC_INT);
}

<INITIAL>{SEC} {
  Within_gen_decl = Within_bin_decl = Within_sec_decl = Within_sos_decl = FALSE;
  Within_sec_decl = TRUE;
  return(SEC_SEC);
}

<INITIAL>{SOS} {
  Within_gen_decl = Within_bin_decl = Within_sec_decl = Within_sos_decl = FALSE;
  Within_sos_decl = TRUE;
  return(SEC_SOS);
}

<INITIAL>{SOST} {
  strcpy(Last_var, (char *)yytext);
  return(SOSTYPE);
}

<INITIAL>{LBL} {
  strcpy(Last_var, (char *)yytext);
  Last_var[strlen(Last_var) - 1] = 0;
  return(VARIABLECOLON);
}

<INITIAL>{VR} {
  strcpy(Last_var, (char *)yytext);
  return(VAR);
}

<INITIAL>{OPLE} {
  return(RE_OPLE);
}

<INITIAL>{OPGE} {
  return(RE_OPGE);
}

<INITIAL>. {
  report(NULL, CRITICAL, "LEX ERROR : %s lineno %d\n", yytext, yylineno);
  return(UNDEFINED);
}

%%
