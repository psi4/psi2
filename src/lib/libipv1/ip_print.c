
/* $Log: ip_print.c,v $
/* Revision 1.4  1994/06/02 02:22:28  seidl
/* using new tmpl now...change .global to .gbl and .local to .lcl
/*
 * Revision 1.1.1.1  1994/05/02  17:05:52  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 1.3  1991/07/30  03:28:45  seidl
 * add rcs log and id
 * */

#include <stdio.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"

#include "ip_print.gbl"
#include "ip_print.lcl"

#define N_INDENT 2

static char rcs_id[] = "$Id: ip_print.c,v 1.4 1994/06/02 02:22:28 seidl Exp $";

GLOBAL_FUNCTION VOID
ip_print_keyword(fp,st)
FILE *fp;
ip_keyword_tree_t *st;
{
  if (st->up) ip_print_keyword(fp,st->up);
  fprintf(fp,"%s:",st->keyword);
  }

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
GLOBAL_FUNCTION VOID
ip_print_tree(fp,tree)
FILE *fp;
ip_keyword_tree_t *tree;
{
  if (!tree) tree = ip_tree;

  ip_print_tree_(fp,tree,0);
  }


/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out.  Indent is used to record how deep in the tree we
 * are, so we know how far to indent things. */
LOCAL_FUNCTION VOID
ip_print_tree_(fp,tree,indent)
FILE *fp;
ip_keyword_tree_t *tree;
int indent;
{
  ip_keyword_tree_t *I;

  I=tree;
  do {
    if (I->value && I->down) {
      ip_warn("ip_print_tree: tree has both value and subtrees - can't print");
      }

    if (!(I->value || I->down)) {
      ip_warn("ip_print_tree: tree has neither value nor subtrees-impossible");
      }

    if (!I->keyword) {
      ip_warn("ip_print_tree: tree has no keyword - impossible");
      }

    ip_indent(fp,indent);
    if (ip_special_characters(I->keyword)) {
      fprintf(fp,"\"%s\"",I->keyword);
      }
    else {
      fprintf(fp,"%s",I->keyword);
      }

    if (I->down) {
      fprintf(fp,": (\n");
      ip_print_tree_(fp,I->down,indent + N_INDENT);
      ip_indent(fp,indent + N_INDENT);
      fprintf(fp,")\n");
      }
    else {
      fprintf(fp," = ");
      ip_print_value_(fp,I->value,indent + strlen(I->keyword) + 3);
      fprintf(fp,"\n");
      }

    } while ((I = I->across) != tree);

  }

LOCAL_FUNCTION VOID
ip_indent(fp,n)
FILE *fp;
int n;
{
  int i;

  for (i=0; i<n; i++) fprintf(fp," ");
  }

GLOBAL_FUNCTION VOID
ip_print_value(fp,value)
FILE *fp;
ip_value_t *value;
{
  if (!value) return;
  ip_print_value_(fp,value,0);
  fprintf(fp,"\n");
  }

LOCAL_FUNCTION VOID
ip_print_value_(fp,value,indent)
FILE *fp;
ip_value_t *value;
int indent;
{
  if (value->type == IP_SCALAR) {
    if (ip_special_characters(value->v.scalar)) {
      fprintf(fp,"\"%s\"",value->v.scalar);
      }
    else {
      fprintf(fp,"%s",value->v.scalar);
      }
    }
  else if (value->type == IP_ARRAY) {
    ip_print_array_(fp,value->v.array,indent);
    }
  else {
    ip_warn("ip_print_value_: bad value type");
    }
  }

LOCAL_FUNCTION VOID
ip_print_array_(fp,array,indent)
FILE *fp;
ip_array_t *array;
int indent;
{
  int i;
  fprintf(fp,"(");
  for (i=0; i<array->n; i++) {
    ip_print_value_(fp,array->values[i],indent);
    if (i != array->n - 1) fprintf(fp," ");
    }
  fprintf(fp,")");
  }

LOCAL_FUNCTION int
ip_special_characters(keyword)
char *keyword;
{
  char *ch=keyword;

  if (!keyword) return 0;
  while (*ch) {
    if (!(  (*ch >= 'a' && *ch <= 'z')
          ||(*ch >= 'A' && *ch <= 'Z')
          ||(*ch >= '0' && *ch <= '9')
          ||(*ch == '+')
          ||(*ch == '-')
          ||(*ch == '.')
          ||(*ch == '_'))) return 1;

    ch++;
    }

  return 0;
  }
