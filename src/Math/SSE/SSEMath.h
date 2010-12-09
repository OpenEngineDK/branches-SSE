// SSE Math Libary
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS)
//
// This program is free software; It is covered by the GNU General
// Public License version 2 or any later version.
// See the GNU General Public License for more details (see LICENSE).
// --------------------------------------------------------------------

#ifndef _OE_SSE_MATH_H_
#define _OE_SSE_MATH_H_

#ifdef __SSE__

#include <Math/SSE/float4.h>

// Header for SSE
#include <xmmintrin.h>
#include <pmmintrin.h>

namespace OpenEngine {
namespace Math {
namespace SSE {

   //Getting a boolean
   FORCE_INLINE void Truth(const float4 m, bool *All)
   {
        float p[4];
        _mm_store_ps((float*)p, m.getData());
        for(int x=0; x<4; x++)
            All[x] = (int)p[x];
   }
   
   template<unsigned short N>
   FORCE_INLINE bool Truth(const float4 m)
   {
        bool Array[4];
        Truth(m, Array);
        for(int x=0; x<N; x++)
        {
            if(!(Array[x]))
                return false;
        }
        return true;
    
   }

   FORCE_INLINE float4 sqrt(float4 m)
   {
       return _mm_sqrt_ss(m.getData());
   }

   FORCE_INLINE float4 sum(float4 elm)
   {
         // Commented out, for testing
/*       #ifdef __SSE3__
           SSE::SSE_float4 m = {0xffffffff, 0xffffffff, 0xffffffff, 0x0};
           SSE::float4 a(_mm_and_ps(elm.getData(), m));
       	   a = _mm_hadd_ps(a.getData(), a.getData());
	       a = _mm_hadd_ps(a.getData(), a.getData());
	       return a;
	   #else
*/	       SSE::float4 m = elm;
	       SSE::float4 r = m;
           for(int x=0; x<3; x++)
           {
               m >>= 1;
	           r += m;
           }
	       return r;
	//   #endif
   }
   
   FORCE_INLINE float4 min(float4 elm)
   {
          SSE::float4 m = elm;
	      SSE::float4 r = m;
          for(int x=0; x<3; x++)
          {
             m <<= 1;
	         r = _mm_min_ps(r.getData(), m.getData());
          }
	      return r;
   }
   
   FORCE_INLINE float4 max(float4 elm)
   {
          SSE::float4 m = elm;
	      SSE::float4 r = m;
          for(int x=0; x<3; x++)
          {
             m <<= 1;
	         r = _mm_max_ps(r.getData(), m.getData());
          }
	      return r;
   }
   
   FORCE_INLINE float4 dotproduct(float4 elm, float4 velm)
   {
        /*
        #ifdef __SSE3__
            SSE::SSE_float4 a = {0xffffffff, 0xffffffff, 0xffffffff, 0x0};
            SSE::float4 m = _mm_and_ps(elm.getData()*velm.getData(), a);
	        m = _mm_hadd_ps(m.getData(), m.getData());
	        m = _mm_hadd_ps(m.getData(), m.getData());
            return m;
         #else
         */
	        SSE::float4 m = elm * velm;
	        SSE::float4 r = m;
            for(int x=0; x<3; x++)
            {
                m >>= 1;
	            r += m;
            }
	        return r;
         //#endif
   }

}  // NS SSE
}  // NS OpenEngine
}  // NS Math

#endif // __SSE__
#endif // _OE_SSE_MATH_H
