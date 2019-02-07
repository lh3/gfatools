%{
// to compile: flex gfa.l && gcc -O2 lex.yy.c -o gfa-chk
#include <stdio.h>
int lineno;
%}

sgn          \+|\-
id           [!-~]+
num          [0-9]+
int          {sgn}?{num}
float        {sgn}?({num}?"."{num}([eE]{sgn}?{num})?)|({num}[eE]{sgn}?{num})

key          [A-Za-z][A-Za-z0-9_]
tag_i        {key}:i:{int}
tag_f        {key}:f:{float}
tag_A        {key}:A:[ -~]
tag_Z        {key}:Z:[ -~]+
tag_Bi       {key}:B:[CsSiI](,{int})+
tag_Bf       {key}:B:f(,{float})+
tag          {tag_i}|{tag_f}|{tag_A}|{tag_Z}|{tag_Bi}|{tag_Bf}

seq          [!-~]+|\*
cigar        ([0-9]+[MIDNS])+
overlap      {num}(:{num})?|{num}:|:{num}|{cigar}

header       H(\t{tag})*
segment      S\t{id}\t{seq}(\t{tag})*
link         L\t{id}\t{sgn}\t{id}\t{sgn}\t{overlap}(\t{tag})*

%%

^[^HSL].*\n    { ++lineno; return 1; }
^{header}\n    { ++lineno; return 2; }
^{segment}\n   { ++lineno; return 3; }
^{link}\n      { ++lineno; return 4; }
^.*\n          { ++lineno; fprintf(stderr, "Parse error at line %d:\n%s", lineno, yytext); return -1; }

%%

int yywrap(void) { return 1; }

int main(void)
{
	while (yylex() != 0);
	return 0;
}
