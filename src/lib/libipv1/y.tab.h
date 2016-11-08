#ifndef BISON_Y_TAB_H
# define BISON_Y_TAB_H

#ifndef YYSTYPE
typedef union {
  ip_value_t *val;
  char *str;
  } yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	T_STRING	257


extern YYSTYPE yylval;

#endif /* not BISON_Y_TAB_H */
