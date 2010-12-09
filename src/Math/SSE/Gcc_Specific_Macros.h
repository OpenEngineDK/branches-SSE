#ifndef _GCC_SPECIFIC_MACROS_H
#define _GCC_SPECIFIC_MACROS_H

#define FORCE_INLINE inline __attribute__((always_inline))
#define ALIGNED(x) __attribute__((aligned(x)))

#endif
