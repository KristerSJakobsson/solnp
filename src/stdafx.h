#ifndef CPP_SOLNP_STDAFX_H
#define CPP_SOLNP_STDAFX_H

/* For the local tests and when imported in C++, we fetch and compile dlib from source
 * The main reason for this is that dlib will detect BLAS and LAPACK on compilation, which reduces rounding errors.
 * The drawback is that we will compile a lot of unneeded code.
*/

#ifndef USE_PRECOMPILED_DLIB
    #include "../library/dlib/matrix.h"
#else
    #include "dlib/matrix.h"
#endif

#endif //CPP_SOLNP_STDAFX_H
