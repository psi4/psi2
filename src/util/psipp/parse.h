#ifndef BISON_Y_TAB_H
# define BISON_Y_TAB_H

#ifndef YYSTYPE
typedef union {
  int ival;
  char *sval;
  } yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	T_FORTRAN	257
# define	T_INCLUDE	258
# define	T_IF	259
# define	T_INACTIVE_IF	260
# define	T_ELIF	261
# define	T_IFDEF	262
# define	T_IFNDEF	263
# define	T_ELSE	264
# define	T_ENDIF	265
# define	T_DEFINE	266
# define	T_UNDEF	267
# define	T_EXPR_DONE	268
# define	T_OR	269
# define	T_AND	270
# define	T_EQ	271
# define	T_NOT_EQ	272
# define	T_GREA_EQ	273
# define	T_LESS_EQ	274
# define	T_NAME	275
# define	T_DEFINED	276
# define	T_STREQ	277


extern YYSTYPE yylval;

#endif /* not BISON_Y_TAB_H */
