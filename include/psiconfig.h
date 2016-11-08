
#ifndef _psiconfig_h
#define _psiconfig_h

typedef int psi_int;
typedef unsigned int psi_uint;
typedef double psi_real;

#ifdef NO_VOID_PTR
typedef char* void_ptr;
#else
typedef void* void_ptr;
#endif

#endif /* _psiconfig_h */
