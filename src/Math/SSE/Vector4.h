// Math vector.
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS)
//
// This program is free software; It is covered by the GNU General
// Public License version 2 or any later version.
// See the GNU General Public License for more details (see LICENSE).
// --------------------------------------------------------------------

#warning MAKE IT COMPILE WHEN -msse3 isn't defined!

#ifndef _OE_VECTOR4_H_
#define _OE_VECTOR4_H_

#ifdef __SSE__

#include <Math/Exceptions.h>
#include <Math/SSE/float4.h>
#include <Math/SSE/SSEMath.h>
#include <Math/SSE/Gcc_Specific_Macros.h>

#include <string>
#include <sstream>
#include <ostream>
#include <math.h>
#include <boost/static_assert.hpp>

// Header for SSE
#include <xmmintrin.h>
#include <pmmintrin.h>

namespace OpenEngine {
namespace Math {

/**
 * Vector.
 *
 * @class Vector Vector.h Math/Vector.h
 * @param N Number of elements
 * @param T Type of elements
 */
template<>
class Vector<4, float> 
{
protected:
    // vector elements
    SSE::float4 elm;

public:
    /**
     * Create zero vector.
     * @code
     * Vector<3,int> v;   // [0, 0, 0]
     * @endcode
     */
    FORCE_INLINE Vector<4, float> ()
    : elm()
    {
    }
    
    /**
     * Create vector from scalar.
     * @code
     * Vector<3,int> v(7);   // [7, 7, 7]
     * @endcode
     *
     * @param s Scalar value in all indexes
     */
    FORCE_INLINE explicit Vector<4, float> (const float s)
    : elm(s)
    {
    }
    
    /**
     * Copy constructor.
     *
     * @param v Vector to copy
     */
    FORCE_INLINE Vector<4, float> (const Vector<4, float> & v)
    : elm(v.elm)
    {
    }
    
    FORCE_INLINE Vector<4, float> (const SSE::float4& v)
    : elm(v)
    {
    }
    
    FORCE_INLINE Vector<4, float> & operator=(const Vector<4, float> & other)
    {
        if (this != &other)
            elm = other.elm;
        return *this;
    }
    
    /**
     * Create vector from array.
     *
     * @param a Array to copy
     */
    FORCE_INLINE explicit Vector<4, float> (const float *a)
    : elm(a)
    {
    }

    /**
     * Constructor for a 4 element vector.
     */
    FORCE_INLINE Vector<4, float> (const float x, const float y, const float z, const float w)
    : elm(x,y,z,w)
    {   
    }
    
    /**
     * Index access to vector elements.
     * @code
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * v[1]                      // 2
     * @endcode
     *
     * @param i Index of element to return
     * @exception IndexOutOfBounds Index is out of bounds.
     * @return Element at index \a i
     */
    FORCE_INLINE float& operator[](const unsigned int i)
    {
#if OE_SAFE
        if (i >= 4)
            throw IndexOutOfBounds(i,0,4);
#endif
        return elm[i];
    }
    
    /**
     * Vector equality.
     */
    FORCE_INLINE bool operator==(const Vector<4, float> & v) const
    {
        if(SSE::Truth<4>(elm==v.elm))
            return false;
        return true;
    }
    
    /**
     * Vector inequality.
     */
    FORCE_INLINE bool operator!=(const Vector<4, float> & v) const
    {
        return !(*this == v);
    }
    
    /**
     * Scalar addition. 
     * This is a commutative operation.
     *
     * @code
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * v + 10                    // [11, 12, 13]
     * 10 + v                    // [11, 12, 13]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator+(const float s) const
    {
        return elm+SSE::float4(s);
    }
    
    /**
     * Vector addition.
     * @code
     * Vector<3,int> u(1,2,3);   // [1, 2, 3]
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * u + v                     // [2, 4, 6]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator+(const Vector<4, float> & v) const
    {
        return elm + v.elm;
    }
    
    /**
     * Scalar subtraction.
     *
     * @code
     * Vector<3,int> v(1,2,3);   // [11, 12, 13]
     * v - 10                    // [1, 2, 3]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator-(const float s) const
    {
        return elm-SSE::float4(s);
    }
    
    /**
     * Vector subtraction.
     * @code
     * Vector<3,int> u(1,2,3);   // [2, 4, 6]
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * u - v                     // [1, 2, 3]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator-(const Vector<4, float> & v) const
    {
        return elm-v.elm;
    }
    
    /**
     * Scalar multiplication.
     * @code
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * v * 10                    // [10, 20, 30]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator*(const float s) const
    {
        return elm*SSE::float4(s);
    }
    
    /**
     * Scalar division.
     * @code
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * v / 10                    // [0.5, 1.0, 1.5]
     * @endcode
     *
     * @todo What should the type of the returning vector be?
     *
     * @param s Scalar to divide by.
     * @exception DivisionByZero Cannot divide by a zero scalar.
     * @return Vector where all elements are divided by \a s
     */
    FORCE_INLINE const Vector<4, float>  operator/(const float s) const
    {
#if OE_SAFE
        if (s == 0)
            throw DivisionByZero();
#endif
        return elm/SSE::float4(s);
    }

    /**
     * Additive inverse.
     * @code
     * Vector<3,int> u(1,2,3);   // [ 1,  2,  3]
     * -u                        // [-1, -2, -3]
     * @endcode
     */
    FORCE_INLINE const Vector<4, float>  operator-() const
    {
        return -elm;
    }
    
    /**
     * Dot/scalar product
     * @code
     * Vector<3,int> u(1,2,3);   // [1, 2, 3]
     * Vector<3,int> v(1,2,3);   // [1, 2, 3]
     * u * v                     // 14
     * @endcode
     */
    FORCE_INLINE float operator*(const Vector<4, float> & v) const
    {
	      SSE::float4 m = elm * v.elm;
	      SSE::float4 r = m;
         for(int x=0; x<3; x++)
         {
            m >>= 1;
	         r += m;
         }
	      return r[0];
    }
    
    /**
     * Destructive scalar addition.
     */
    FORCE_INLINE void operator+=(const float s)
    {
        elm += SSE::float4(s);
    }
    
    /**
     * Destructive vector addition.
     */
    FORCE_INLINE void operator+=(const Vector<4, float> & v)
    {
        elm += v.elm;
    }
    
    /**
     * Destructive scalar multiplication.
     */
    FORCE_INLINE void operator*=(const float s)
    {
        elm *= SSE::float4(s);
    }
    
    /**
     * Destructive scalar subtraction.
     */
    FORCE_INLINE void operator-=(const float s)
    {
        elm -= SSE::float4(s);
    }
    
    /**
     * Destructive scalar division.
     *
     * @param s Scalar to divide by.
     * @exception DivisionByZero Cannot divide by a zero scalar.
     * @return Vector where all elements are divided by \a s
     */
    FORCE_INLINE void operator/=(const float s)
    {
#if OE_SAFE
        if (s == 0)
            throw DivisionByZero();
#endif
        elm /= SSE::float4(s);
    }
    
    /**
     * Is this the zero vector.
     *
     * @return True if all elements are zero, false otherwise.
     */
    FORCE_INLINE bool IsZero() const
    {
        if(SSE::Truth<4>(SSE::float4()==elm))
            return false;
        return true;
    }
    
    /**
     * Length/modulo of vector.
     *
     * @return Vector length
     */
    FORCE_INLINE float GetLength() const
    {
        // * is skalar/dot product, not member by member multiplikation
        // which returns a float, hence the float sqrt
        return sqrt((float) ((*this) * (*this)));
    }
    
    /**
     * Length/modulo squared of vector.
     *
     * @return Vector length squared
     */
    FORCE_INLINE float GetLengthSquared() const
    {
        // * is skalar/dot product, not member by member multiplikation
        // which returns a float
        return (float) ((*this) * (*this));
    }
    
    /**
     * Normalized vector.
     *
     * @return Normalized vector
     */
    FORCE_INLINE Vector<4, float>  GetNormalize() const
    {
        Vector<4, float>  v(*this);
        v.Normalize();
        return v;
    }
    
    /**
     * Destructively Normalize vector.
     *
     * @exception ArithmeticException Normalizing zero vector
     */
    FORCE_INLINE void Normalize()
    {
        float norm = GetLength();
#if OE_SAFE
        if (norm == 0)
            throw ArithmeticException("Can not normalize the zero vector.");
#endif
      #warning note: CHECK THE A FLOAT USING "!=" COMPARSION
        for(unsigned int i=0; i<4; i++)
            if(elm[i]!= 0) //Checking a float using "!=" ??
	            elm[i]=(((float)elm[i])/norm);
    }
    
    /**
     * Maximum element value.
     */
    FORCE_INLINE float Max() const
    {
	      SSE::float4 m = elm;
	      SSE::float4 r = m;
         for(int x=0; x<3; x++)
         {
            m <<= 1;
            #warning NOTE: IMPLEMENT THIS!
	         r = _mm_max_ps(r.getData(), m.getData()); //IMPLEMENT IN SSE::float4
         }
	      return r[0];
    }
    
    /**
     * Minimum element value.
     */
    FORCE_INLINE float Min() const
    {
	      SSE::float4 m = elm;
	      SSE::float4 r = m;
         for(int x=0; x<3; x++)
         {
            m <<= 1;
            #warning NOTE: IMPLEMENT THIS!
	         r = _mm_min_ps(r.getData(), m.getData()); //IMPLEMENT IN SSE::float4
         }
	      return r[0];
    }
    
    /**
     * Index of maximum element value.
     */
    FORCE_INLINE unsigned int MaxIndex() const
    {
      unsigned int m = 0;
        for (unsigned int i=1; i<4; i++)
            if (elm.get(i) > elm.get(m)) m = i;
        return m;
    }
    
    /**
     * Index of minimum element value.
     */
    FORCE_INLINE unsigned int MinIndex() const
    {
        int m = 0;
        for (unsigned int i=1; i<4; i++)
            if (elm.get(i) < elm.get(m)) m = i;
        return m;
    }
    
    /**
     * Sum of elements.
     *
     * @return Sum of elements
     */
    FORCE_INLINE float Sum() const
    {
    /* // Alternative way
	      register SSE::float4 m = elm;
	      register SSE::float4 r = m;
         for(int x=0; x<3; x++)
         {
            m >>= 1;
	         r += m;
         }
	      return r[0];
   */
         SSE::float4 a = elm;
         #warning NOTE: IMPLEMENT THIS!
   	   a = _mm_hadd_ps(a.getData(), a.getData());  //IMPLEMENT IN SSE::float4
	      a = _mm_hadd_ps(a.getData(), a.getData());
	      return a[0];
    }
    /**
     * Get vector with integer entries.
     *
     * @return Vector of integers.
     */
    FORCE_INLINE Vector<4,int> ToInt() const
    {
        Vector<4,int> v;
        for(unsigned int i=0; i<4; i++)
            v[i] = (int)elm.get(i);
        return v;
    }
    
    /**
     * Get vector with floating point entries.
     *
     * @return Vector of floats.
     */
    FORCE_INLINE Vector<4, float>  ToFloat() const
    {
        Vector<4, float>  v(*this);
        return v;
    }
    
    /**
     * Get vector with double precision entries.
     *
     * @return Vector of doubles.
     */
    FORCE_INLINE Vector<4,double> ToDouble() const
    {
        Vector<4,double> v;
        for (unsigned int i=0; i<4; i++)
            v[i] = (double)elm.get(i);
        return v;
    }
    
    /**
     * Create array from vector.
     *
     * @param a Array to populate
     */
    FORCE_INLINE void ToArray(float *a) const
    {
        elm.ToArray(a);
    }
    
    /**
     * Returns a pointer directly into the vectors data core.
     *
     * @return The data pointer.
     */
    FORCE_INLINE float* ToArray() const
    {
        return (float*) &elm;
    }
    
    /**
     * String representation.
     * Ex. [1, 2, 3]
     *
     * @return Vector as string
     */
    FORCE_INLINE std::string ToString() const
    {
        return elm.ToString<4>();
    }
    
    /**
     * Get value of index (non modifiable).
     * 
     * @param i Index of element to return
     * @exception IndexOutOfBounds Index is out of bounds.
     * @return Element at index \a i
     */
    FORCE_INLINE float Get(const unsigned int i) const
    { 
#if OE_SAFE
        if (i >= 4)
            throw IndexOutOfBounds(i,0,4);
#endif
        return elm.get(i); 
    }
};

}  // NS OpenEngine
}  // NS Math

#endif // __SSE__
#endif // _OE_VECTOR4_H

