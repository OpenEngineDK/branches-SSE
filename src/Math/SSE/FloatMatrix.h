// Math matrix specialization
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS)
//
// This program is free software; It is covered by the GNU General
// Public License version 2 or any later version.
// See the GNU General Public License for more details (see LICENSE).
//--------------------------------------------------------------------

#ifndef _OE_FLOAT_MATRIX_H_
#define _OE_FLOAT_MATRIX_H_

#ifdef __SSE__

#include <Math/SSE/Matrix44.h>
//#include <Math/SSE/Matrix34.h>

namespace OpenEngine {
namespace Math {

class FloatMatrix : public SSE::Matrix44
{
   public:
   FORCE_INLINE FloatMatrix(const float f) : Matrix44(f)
   {}
};
/*
template<>
class Matrix<3,4,float> : Matrix34
{
}

template<>
class Matrix<3,3,float> : Matrix34 // Matrix34 on purpose, as float4 takes 4floats no matter what
{
}
*/
}  // NS OpenEngine
}  // NS Math
#endif // __SSE__
#endif // _OE_FLOAT_MATRIX_H_
