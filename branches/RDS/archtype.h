
        /* File  archtype.h */

//#ifndef ARCHTYPE
//#include "rcsid.h"

//GENERATE_RCSID(archtype_rcsid, 
//	"$Id: archtype.h,v 1.3 1994/06/07 18:18:19 pmontgom Exp $");

#ifdef CRAY
#define ARCH_CRAY 1
#define ARCHTYPE "CRAY"
#define VECTOR_REGISTER_LENGTH 128
#else
#define ARCH_CRAY 0
#endif


#ifdef sun
#define ARCH_SPARC 1
#define ARCHTYPE "SPARC"
#else
#define ARCH_SPARC 0
#endif


#ifdef __sgi
#define ARCH_SGI 1
#define ARCHTYPE "SGI"
#else
#define ARCH_SGI 0
#endif

//#ifndef ARCHTYPE
//   #error in archtype.h -- unknown architecture
//#endif

//#endif /* ARCHTYPE */
