// Define SSE_NEEDED if its required to make something work
// Define SSE2_IMPROVE if SSE2 would improve performence, but it would work on SSE1

#ifndef _OE_SSE_STATUS_H
#define _OE_SSE_STATUS_H

#ifndef __SSE__
    #ifdef SSE_NEEDED
        #error NOTE: __SSE__ not enabled (Try to pass -msse to compiler)
    #endif //SSE_NEEDED
#endif //__SSE__

#ifndef __SSE2__
    #ifdef SSE2_NEEDED
        #error NOTE: __SSE2__ not enabled (Try to pass -msse2 to compiler)
    #else
        #ifdef SSE2_IMPROVE
            #warning NOTE:  __SSE2__ not enabled (Passing -msse2 to compiler might peform better)
        #endif
    #endif
#endif //__SSE2__

#ifndef __SSE3__
    #ifdef SSE3_NEEDED
        #error NOTE: __SSE3__ not enabled (Try to pass -msse3 to compiler)
    #else
        #ifdef SSE3_IMPROVE
            #warning NOTE:  __SSE3__ not enabled (Passing -msse3 to compiler might peform better)
        #endif
    #endif
#endif //__SSE3__

#ifndef __SSE4__
    #ifdef SSE4_NEEDED
        #error NOTE: __SSE4__ not enabled (Try to pass -msse4 to compiler)
    #else
        #ifdef SSE4_IMPROVE
            #warning NOTE:  __SSE4__ not enabled (Passing -msse4 to compiler might peform better)
        #endif
    #endif
#endif //__SSE4__

#ifndef __SSE5__
    #ifdef SSE5_NEEDED
        #error NOTE: __SSE5__ not enabled (Try to pass -msse5 to compiler)
    #else
        #ifdef SSE5_IMPROVE
            #warning NOTE:  __SSE5__ not enabled (Passing -msse5 to compiler might peform better)
        #endif
    #endif
#endif //__SSE5__

#endif //_OE_SSE_STATUS_H
