// Math matrix44
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS)
//
// This program is free software; It is covered by the GNU General
// Public License version 2 or any later version.
// See the GNU General Public License for more details (see LICENSE).
//--------------------------------------------------------------------

#ifndef _OE_MATRIX44_H_
#define _OE_MATRIX44_H_

#ifdef __SSE__

__attribute__((aligned(16))) const long _Sign_PNNP[4] = { 0x00000000, 0x80000000, 0x80000000, 0x00000000 };
#define _mm_ror_ps(vec, i) (((i)%4) ? (_mm_shuffle_ps(vec,vec, _MM_SHUFFLE((unsigned char)(i+3)%4,(unsigned char)(i+2)%4,(unsigned char)(i+1)%4,(unsigned char)(i+0)%4))) : (vec))

#include <string>
#include <sstream>
#include <boost/static_assert.hpp>
#include <Math/Vector.h>
#include <Math/SSE/float4.h>

// Header for SSE
#include <xmmintrin.h>
#include <pmmintrin.h>

#ifdef _OE_FLOAT4_INCLUDED

namespace OpenEngine {
namespace Math {

template<>
class Matrix<4,4,float>
{
private:
    // matrix elements
    SSE::float4 elm[4];

public:
    /**
     * Create identity matrix.
     * If the dimension is not NxN zero-rows will appear.
     *
     * @code
     * Matrix<3,3,int> m;   // [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
     * @endcode
     */
    Matrix<4,4,float>()
    {
        SSE::float4 Temp(1.0f, 0.0f, 0.0f, 0.0f);
        for(int x=0; x<4; x++)
        {
            elm[x] = Temp;
            Temp>>=1;
        }
    }
    /**
     * Copy constructor.
     *
     * @param m Matrix to copy
     */
    Matrix<4,4,float>(const Matrix<4,4,float>& m)
    {
        for(int x=0; x<4; x++)
            elm[x] = m.elm[x];
    }
    /**
     * Create matrix from scalar.
     * @code
     * Matrix<2,3,int> m(7);   // [(7, 7, 7), (7, 7, 7)]
     * @endcode
     *
     * @param s Scalar value in all indexes
     */
    explicit Matrix<4,4,float>(const float s)
    {
        SSE::float4 Temp(s);
        for(int x=0; x<4; x++)
            elm[x] = Temp;
    }
    /**
     * Create matrix from array.
     *
     * @param a Array to create fro
     */
    explicit Matrix<4,4,float>(const float* a)
    {
        for(int x=0; x<4; x++)
            elm[x] = SSE::float4(a+(4*x));
    }
    
    /**
     * Constructor for a 4x4 matrix.
     */
    Matrix<4,4,float>(const float a, const float b, const float c, const float d, 
           const float e, const float f, const float g, const float h,
           const float i, const float j, const float k, const float l, 
           const float n, const float m, const float o, const float p) {

        elm[0] = SSE::float4(a,b,c,d);
        elm[1] = SSE::float4(e,f,g,h);
        elm[2] = SSE::float4(i,j,k,l);
        elm[3] = SSE::float4(n,m,o,p);
    }
    
    /**
     * Index access to matrix elements.
     * @code
     * Matrix<2,2,int> m(1,2, 3,4);   // [(1, 2), (3, 4)]
     * m(1,0)                         // 3
     * @endcode
     *
     * @param i Row index
     * @param j Column index
     * @return Element at index \a (i,j)
     */
    float& operator()(const unsigned int i, const unsigned int j)
    {
#if OE_SAFE
        if (i >= 4)
            throw IndexOutOfBounds(i,0,4);
        if (j >= 4)
            throw IndexOutOfBounds(j,0,4);
#endif
        return elm[i][j];
    }
    /**
     * Matrix equality.
     * True if all index element are identical.
     */
    bool operator==(const Matrix<4,4,float>& m) const
    {
        for (unsigned int i=0; i<4; i++)
            for (unsigned int j=0; j<4; j++)
                if (elm[i].get(j) != m.elm[i].get(j))
                    return false;
        return true;
    }
    /**
     * Matrix inequality.
     * True if one or more index elements differ.
     */
    bool operator!=(const Matrix<4,4,float>& m) const
    {
        return !(*this == m);
    }
    /**
     * Matrix multiplication.
     * @code
     * Matrix<2,2,int> a(1,2,3,4);     // [(1,  2), ( 3,  4)]
     * a * a;                          // [(7, 10), (15, 22)]
     * @endcode
     */
    Matrix<4,4,float> operator*(const Matrix<4,4,float> m)
    {
        Matrix<4,4,float> r;
        for (unsigned int i=0; i<4; i++) 
            for (unsigned int j=0; j<4; j++)
            {
                float s = 0;
                for (unsigned int t=0; t<4; t++)
                    s += elm[i].get(t) * m.elm[t].get(j);
                r.elm[i][j] = s;
            }
        return r;
    }
    
    /**
     * Matrix-Vector multiplication of 4*4 matrices.
     */
    const Vector<4,float> operator*(const Vector<4,float> v)
    {
        Vector<4,float> r;
        r[0] = elm[0].get(0) * v.Get(0) + elm[0].get(1) * v.Get(1) + elm[0].get(2) * v.Get(2) + elm[0].get(3) * v.Get(3);
        r[1] = elm[1].get(0) * v.Get(0) + elm[1].get(1) * v.Get(1) + elm[1].get(2) * v.Get(2) + elm[1].get(3) * v.Get(3);
        r[2] = elm[2].get(0) * v.Get(0) + elm[2].get(1) * v.Get(1) + elm[2].get(2) * v.Get(2) + elm[2].get(3) * v.Get(3);
        r[3] = elm[3].get(0) * v.Get(0) + elm[3].get(1) * v.Get(1) + elm[3].get(2) * v.Get(2) + elm[3].get(3) * v.Get(3);
        return r;
    }
    /**
     * Scalar mult
     */
    const Matrix<4,4,float> operator*(const float s)
    {
        Matrix<4,4,float> r;
        SSE::float4 Temp(s);
        for(int x=0; x<4; x++)
            r.elm[x] = Temp*elm[x];
        return r;
    }
    /**
     * Destructive matrix addition.
     */
    void operator+=(Matrix<4,4,float> m)
    {
        for(int x=0; x<4; x++)
            elm[x] += m.elm[x];
    }
    /**
     * Matrix addition.
     */
    const Matrix<4,4,float> operator+(Matrix<4,4,float> m)
    {
        Matrix<4,4,float> r;
        for(int x=0; x<4; x++)
            r.elm[x] = elm[x] + m.elm[x];
        return r;
    }
    /**
     * Get matrix row vector.
     * @code
     * Matrix<2,2,int> m(1,2, 3,4);   // [(1, 2), (3, 4)]
     * m[1]                           // [3, 4]
     * @endcode
     *
     * @see GetRow()
     * @param i Row index
     * @return Row vector
     */
    Vector<4,float> operator[](const unsigned int i) {
        return this->GetRow(i);
    }
    /**
     * Set matrix row vector.
     *
     * @param i Row index
     * @param r Row vector
     */
    void SetRow(const unsigned int i, const Vector<4,float> r)
    {
        r.ToArray((float*)&elm[i]);
    }
    /**
     * Get matrix row vector.
     *
     * @param i Row index
     * @return Row vector
     */
    Vector<4,float> GetRow(const unsigned int i)
    {
        return Vector<4,float>((float*)&elm[i]);
    }
    /**
     * Get matrix column vector.
     *
     * @param j Column index
     * @return Column vector
     */
    Vector<4,float> GetColumn(const unsigned int j)
    {
        Vector<4,float> v;
        for (int i=0; i<4; i++)
            v[i] = elm[i].get(j);
        return v;
    }
    /**
     * Matrix trace.
     * Only defined for NxN matrices.
     *
     * @return Sum of elements in the main diagonal
     */
    float Trace()
    {
        return (elm[0] + (elm[1]<<1) + (elm[2]<<2) + (elm[3]<<3)).get(0);
        /* SECOND SOLUTION
        float t = 0;
        for (unsigned int i=0; i<4; i++)
            t += elm[i].get(i);
        return t;
        */
    }
    /**
     * Transpose matrix.
     * Note that this is a destructive operation and only works on
     * square matrices.
     */
    void Transpose()
    {
        float tmp;
        for (unsigned int i=1; i<4; i++)
            for (unsigned int j=0; j<i; j++)
            {
                tmp = elm[i].get(j);
                elm[i][j] = elm[j].get(i);
                elm[j][i] = tmp;
            }
    }

    /**
     * Returns the determinant of a matrix.
     * Copied from:
     *     http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm     
     * look at ?    
     *     http://www.google.com/codesearch/p?hl=en&sa=N&cd=3&ct=rc#xrFFuTdR3w0/tulip-2.0.2/library/tulip/include/tulip/cxx/Matrix.cxx&q=matrix%20inverse%20lang:c%2B%2B&l=13
     */
    float GetDeterminant() const
    {
    /*
        float value;
        value =
                elm[0].get(3) * elm[1].get(2) * elm[2].get(1) * elm[3].get(0)-elm[0].get(2) * elm[1].get(3) * elm[2].get(1) * elm[3].get(0)-elm[0].get(3) * elm[1].get(1) * elm[2].get(2) * elm[3].get(0)+elm[0].get(1) * elm[1].get(3)    * elm[2].get(2) * elm[3].get(0)+
                elm[0].get(2) * elm[1].get(1) * elm[2].get(3) * elm[3].get(0)-elm[0].get(1) * elm[1].get(2) * elm[2].get(3) * elm[3].get(0)-elm[0].get(3) * elm[1].get(2) * elm[2].get(0) * elm[3].get(1)+elm[0].get(2) * elm[1].get(3)    * elm[2].get(0) * elm[3].get(1)+
                elm[0].get(3) * elm[1].get(0) * elm[2].get(2) * elm[3].get(1)-elm[0].get(0) * elm[1].get(3) * elm[2].get(2) * elm[3].get(1)-elm[0].get(2) * elm[1].get(0) * elm[2].get(3) * elm[3].get(1)+elm[0].get(0) * elm[1].get(2)    * elm[2].get(3) * elm[3].get(1)+
                elm[0].get(3) * elm[1].get(1) * elm[2].get(0) * elm[3].get(2)-elm[0].get(1) * elm[1].get(3) * elm[2].get(0) * elm[3].get(2)-elm[0].get(3) * elm[1].get(0) * elm[2].get(1) * elm[3].get(2)+elm[0].get(0) * elm[1].get(3)    * elm[2].get(1) * elm[3].get(2)+
                elm[0].get(1) * elm[1].get(0) * elm[2].get(3) * elm[3].get(2)-elm[0].get(0) * elm[1].get(1) * elm[2].get(3) * elm[3].get(2)-elm[0].get(2) * elm[1].get(1) * elm[2].get(0) * elm[3].get(3)+elm[0].get(1) * elm[1].get(2)    * elm[2].get(0) * elm[3].get(3)+
                elm[0].get(2) * elm[1].get(0) * elm[2].get(1) * elm[3].get(3)-elm[0].get(0) * elm[1].get(2) * elm[2].get(1) * elm[3].get(3)-elm[0].get(1) * elm[1].get(0) * elm[2].get(2) * elm[3].get(3)+elm[0].get(0) * elm[1].get(1)    * elm[2].get(2) * elm[3].get(3);
            return value;
*/

    __m128 Va,Vb,Vc;
    __m128 r1,r2,r3,t1,t2,sum;
    __m128 Det;

    // First, Let's calculate the first four minterms of the first line
    t1 = elm[4].getData();
    t2 = _mm_ror_ps(elm[3].getData(),1); 
    Vc = _mm_mul_ps(t2,_mm_ror_ps(t1,0));                   // V3'·V4
    Va = _mm_mul_ps(t2,_mm_ror_ps(t1,2));                   // V3'·V4"
    Vb = _mm_mul_ps(t2,_mm_ror_ps(t1,3));                   // V3'·V4^

    r1 = _mm_sub_ps(_mm_ror_ps(Va,1),_mm_ror_ps(Vc,2));     // V3"·V4^ - V3^·V4"
    r2 = _mm_sub_ps(_mm_ror_ps(Vb,2),_mm_ror_ps(Vb,0));     // V3^·V4' - V3'·V4^
    r3 = _mm_sub_ps(_mm_ror_ps(Va,0),_mm_ror_ps(Vc,1));     // V3'·V4" - V3"·V4'

    Va = _mm_ror_ps(elm[2].getData(),1);
    sum = _mm_mul_ps(Va,r1);
    Vb = _mm_ror_ps(Va,1);
    sum = _mm_add_ps(sum,_mm_mul_ps(Vb,r2));
    Vc = _mm_ror_ps(Vb,1);
    sum = _mm_add_ps(sum,_mm_mul_ps(Vc,r3));

    // Now we can calculate the determinant:
    Det = _mm_mul_ps(sum,elm[1].getData());
    Det = _mm_add_ps(Det,_mm_movehl_ps(Det,Det));
    Det = _mm_sub_ss(Det,_mm_shuffle_ps(Det,Det,1));
    return *(float*)&Det;
    }

    Matrix<4,4,float> GetInverse()
    {
        /*
        float d = GetDeterminant();
        Matrix<4,4,float> b;
        
            b(0,0) = elm[1].get(2)*elm[2].get(3)*elm[3].get(1) - elm[1].get(3)*elm[2].get(2)*elm[3].get(1) + elm[1].get(3)*elm[2].get(1)*elm[3].get(2) - elm[1].get(1)*elm[2].get(3)*elm[3].get(2) - elm[1].get(2)*elm[2].get(1)*elm[3].get(3) + elm[1].get(1)*elm[2].get(2)*elm[3].get(3);
            b(0,1) = elm[0].get(3)*elm[2].get(2)*elm[3].get(1) - elm[0].get(2)*elm[2].get(3)*elm[3].get(1) - elm[0].get(3)*elm[2].get(1)*elm[3].get(2) + elm[0].get(1)*elm[2].get(3)*elm[3].get(2) + elm[0].get(2)*elm[2].get(1)*elm[3].get(3) - elm[0].get(1)*elm[2].get(2)*elm[3].get(3);
            b(0,2) = elm[0].get(2)*elm[1].get(3)*elm[3].get(1) - elm[0].get(3)*elm[1].get(2)*elm[3].get(1) + elm[0].get(3)*elm[1].get(1)*elm[3].get(2) - elm[0].get(1)*elm[1].get(3)*elm[3].get(2) - elm[0].get(2)*elm[1].get(1)*elm[3].get(3) + elm[0].get(1)*elm[1].get(2)*elm[3].get(3);
            b(0,3) = elm[0].get(3)*elm[1].get(2)*elm[2].get(1) - elm[0].get(2)*elm[1].get(3)*elm[2].get(1) - elm[0].get(3)*elm[1].get(1)*elm[2].get(2) + elm[0].get(1)*elm[1].get(3)*elm[2].get(2) + elm[0].get(2)*elm[1].get(1)*elm[2].get(3) - elm[0].get(1)*elm[1].get(2)*elm[2].get(3);
            b(1,0) = elm[1].get(3)*elm[2].get(2)*elm[3].get(0) - elm[1].get(2)*elm[2].get(3)*elm[3].get(0) - elm[1].get(3)*elm[2].get(0)*elm[3].get(2) + elm[1].get(0)*elm[2].get(3)*elm[3].get(2) + elm[1].get(2)*elm[2].get(0)*elm[3].get(3) - elm[1].get(0)*elm[2].get(2)*elm[3].get(3);
            b(1,1) = elm[0].get(2)*elm[2].get(3)*elm[3].get(0) - elm[0].get(3)*elm[2].get(2)*elm[3].get(0) + elm[0].get(3)*elm[2].get(0)*elm[3].get(2) - elm[0].get(0)*elm[2].get(3)*elm[3].get(2) - elm[0].get(2)*elm[2].get(0)*elm[3].get(3) + elm[0].get(0)*elm[2].get(2)*elm[3].get(3);
            b(1,2) = elm[0].get(3)*elm[1].get(2)*elm[3].get(0) - elm[0].get(2)*elm[1].get(3)*elm[3].get(0) - elm[0].get(3)*elm[1].get(0)*elm[3].get(2) + elm[0].get(0)*elm[1].get(3)*elm[3].get(2) + elm[0].get(2)*elm[1].get(0)*elm[3].get(3) - elm[0].get(0)*elm[1].get(2)*elm[3].get(3);
            b(1,3) = elm[0].get(2)*elm[1].get(3)*elm[2].get(0) - elm[0].get(3)*elm[1].get(2)*elm[2].get(0) + elm[0].get(3)*elm[1].get(0)*elm[2].get(2) - elm[0].get(0)*elm[1].get(3)*elm[2].get(2) - elm[0].get(2)*elm[1].get(0)*elm[2].get(3) + elm[0].get(0)*elm[1].get(2)*elm[2].get(3);
            b(2,0) = elm[1].get(1)*elm[2].get(3)*elm[3].get(0) - elm[1].get(3)*elm[2].get(1)*elm[3].get(0) + elm[1].get(3)*elm[2].get(0)*elm[3].get(1) - elm[1].get(0)*elm[2].get(3)*elm[3].get(1) - elm[1].get(1)*elm[2].get(0)*elm[3].get(3) + elm[1].get(0)*elm[2].get(1)*elm[3].get(3);
            b(2,1) = elm[0].get(3)*elm[2].get(1)*elm[3].get(0) - elm[0].get(1)*elm[2].get(3)*elm[3].get(0) - elm[0].get(3)*elm[2].get(0)*elm[3].get(1) + elm[0].get(0)*elm[2].get(3)*elm[3].get(1) + elm[0].get(1)*elm[2].get(0)*elm[3].get(3) - elm[0].get(0)*elm[2].get(1)*elm[3].get(3);
            b(2,2) = elm[0].get(1)*elm[1].get(3)*elm[3].get(0) - elm[0].get(3)*elm[1].get(1)*elm[3].get(0) + elm[0].get(3)*elm[1].get(0)*elm[3].get(1) - elm[0].get(0)*elm[1].get(3)*elm[3].get(1) - elm[0].get(1)*elm[1].get(0)*elm[3].get(3) + elm[0].get(0)*elm[1].get(1)*elm[3].get(3);
            b(2,3) = elm[0].get(3)*elm[1].get(1)*elm[2].get(0) - elm[0].get(1)*elm[1].get(3)*elm[2].get(0) - elm[0].get(3)*elm[1].get(0)*elm[2].get(1) + elm[0].get(0)*elm[1].get(3)*elm[2].get(1) + elm[0].get(1)*elm[1].get(0)*elm[2].get(3) - elm[0].get(0)*elm[1].get(1)*elm[2].get(3);
            b(3,0) = elm[1].get(2)*elm[2].get(1)*elm[3].get(0) - elm[1].get(1)*elm[2].get(2)*elm[3].get(0) - elm[1].get(2)*elm[2].get(0)*elm[3].get(1) + elm[1].get(0)*elm[2].get(2)*elm[3].get(1) + elm[1].get(1)*elm[2].get(0)*elm[3].get(2) - elm[1].get(0)*elm[2].get(1)*elm[3].get(2);
            b(3,1) = elm[0].get(1)*elm[2].get(2)*elm[3].get(0) - elm[0].get(2)*elm[2].get(1)*elm[3].get(0) + elm[0].get(2)*elm[2].get(0)*elm[3].get(1) - elm[0].get(0)*elm[2].get(2)*elm[3].get(1) - elm[0].get(1)*elm[2].get(0)*elm[3].get(2) + elm[0].get(0)*elm[2].get(1)*elm[3].get(2);
            b(3,2) = elm[0].get(2)*elm[1].get(1)*elm[3].get(0) - elm[0].get(1)*elm[1].get(2)*elm[3].get(0) - elm[0].get(2)*elm[1].get(0)*elm[3].get(1) + elm[0].get(0)*elm[1].get(2)*elm[3].get(1) + elm[0].get(1)*elm[1].get(0)*elm[3].get(2) - elm[0].get(0)*elm[1].get(1)*elm[3].get(2);
            b(3,3) = elm[0].get(1)*elm[1].get(2)*elm[2].get(0) - elm[0].get(2)*elm[1].get(1)*elm[2].get(0) + elm[0].get(2)*elm[1].get(0)*elm[2].get(1) - elm[0].get(0)*elm[1].get(2)*elm[2].get(1) - elm[0].get(1)*elm[1].get(0)*elm[2].get(2) + elm[0].get(0)*elm[1].get(1)*elm[2].get(2);
            
        return b * (1.0/d);
        */

    __m128 A = _mm_movelh_ps(elm[0].getData(), elm[1].getData()),    // the four sub-matrices 
            B = _mm_movehl_ps(elm[1].getData(), elm[0].getData()),
            C = _mm_movelh_ps(elm[2].getData(), elm[3].getData()),
            D = _mm_movehl_ps(elm[3].getData(), elm[2].getData());
    __m128 iA, iB, iC, iD,					// partial inverse of the sub-matrices
            DC, AB;
    __m128 dA, dB, dC, dD;                 // determinant of the sub-matrices
    __m128 det, d, d1, d2;
    __m128 rd;

    //  AB = A# * B
    AB = _mm_mul_ps(_mm_shuffle_ps(A,A,0x0F), B);
    AB -= (__m128)_mm_mul_ps(_mm_shuffle_ps(A,A,0xA5), _mm_shuffle_ps(B,B,0x4E));
    //  DC = D# * C
    DC = _mm_mul_ps(_mm_shuffle_ps(D,D,0x0F), C);
    DC -= (__m128)_mm_mul_ps(_mm_shuffle_ps(D,D,0xA5), _mm_shuffle_ps(C,C,0x4E));

    //  dA = |A|
    dA = _mm_mul_ps(_mm_shuffle_ps(A, A, 0x5F),A);
    dA = _mm_sub_ss(dA, _mm_movehl_ps(dA,dA));
    //  dB = |B|
    dB = _mm_mul_ps(_mm_shuffle_ps(B, B, 0x5F),B);
    dB = _mm_sub_ss(dB, _mm_movehl_ps(dB,dB));

    //  dC = |C|
    dC = _mm_mul_ps(_mm_shuffle_ps(C, C, 0x5F),C);
    dC = _mm_sub_ss(dC, _mm_movehl_ps(dC,dC));
    //  dD = |D|
    dD = _mm_mul_ps(_mm_shuffle_ps(D, D, 0x5F),D);
    dD = _mm_sub_ss(dD, _mm_movehl_ps(dD,dD));

    //  d = trace(AB*DC) = trace(A#*B*D#*C)
    d = _mm_mul_ps(_mm_shuffle_ps(DC,DC,0xD8),AB);

    //  iD = C*A#*B
    iD = _mm_mul_ps(_mm_shuffle_ps(C,C,0xA0), _mm_movelh_ps(AB,AB));
    iD += (__m128)_mm_mul_ps(_mm_shuffle_ps(C,C,0xF5), _mm_movehl_ps(AB,AB));
    //  iA = B*D#*C
    iA = _mm_mul_ps(_mm_shuffle_ps(B,B,0xA0), _mm_movelh_ps(DC,DC));
    iA += (__m128)_mm_mul_ps(_mm_shuffle_ps(B,B,0xF5), _mm_movehl_ps(DC,DC));

    //  d = trace(AB*DC) = trace(A#*B*D#*C) [continue]
    d = _mm_add_ps(d, _mm_movehl_ps(d, d));
    d = _mm_add_ss(d, _mm_shuffle_ps(d, d, 1));
    d1 = dA*dD;
    d2 = dB*dC;

    //  iD = D*|A| - C*A#*B
    iD = D*_mm_shuffle_ps(dA,dA,0) - iD;

    //  iA = A*|D| - B*D#*C;
    iA = A*_mm_shuffle_ps(dD,dD,0) - iA;

    //  det = |A|*|D| + |B|*|C| - trace(A#*B*D#*C)
    det = d1+d2-d;
    rd = (__m128)(SSE::float4((float)1.0f).getData()/det);
    #ifdef ZERO_SINGULAR
        rd = _mm_and_ps(_mm_cmpneq_ss(det,_mm_setzero_ps()), rd);
    #endif

    //  iB = D * (A#B)# = D*B#*A
    iB = _mm_mul_ps(D, _mm_shuffle_ps(AB,AB,0x33));
    iB -= (__m128)_mm_mul_ps(_mm_shuffle_ps(D,D,0xB1), _mm_shuffle_ps(AB,AB,0x66));
    //  iC = A * (D#C)# = A*C#*D
    iC = _mm_mul_ps(A, _mm_shuffle_ps(DC,DC,0x33));
    iC -= (__m128)_mm_mul_ps(_mm_shuffle_ps(A,A,0xB1), _mm_shuffle_ps(DC,DC,0x66));

    rd = _mm_shuffle_ps(rd,rd,0);
    rd = _mm_xor_ps(rd,(*(__m128*)&_Sign_PNNP));

    //  iB = C*|B| - D*B#*A
    iB = C*_mm_shuffle_ps(dB,dB,0) - iB;

    //  iC = B*|C| - A*C#*D;
    iC = B*_mm_shuffle_ps(dC,dC,0) - iC;

    //  iX = iX / det
    iA *= rd;
    iB *= rd;
    iC *= rd;
    iD *= rd;

    Matrix<4,4,float> b;
    b.elm[0] = _mm_shuffle_ps(iA,iB,0x77);
    b.elm[1] = _mm_shuffle_ps(iA,iB,0x22);
    b.elm[2] = _mm_shuffle_ps(iC,iD,0x77);
    b.elm[3] = _mm_shuffle_ps(iC,iD,0x22);

    return b;
    }
    

    /**
     * Get a matrix expanded by one column and one row.
     * The column and row will consist of zero elements and a one in
     * the diagonal entry.
     *
     * @code
     * Matrix<2,2,int> m(1,2, 3,4);  // [(1,2), (3,4)]
     * m.GetExpanded()               // [(1,2,0), (3,4,0), (0,0,1)]
     * @endcode
     *
     * @return Reduced matrix.
     */
    Matrix<5,5,float> GetExpanded()
    {
        Matrix<5,5,float> m;
        for (unsigned int i=0; i < 4; i++)
            for (unsigned int j=0; j < 4; j++)
                m(i,j) = elm[i].get(j);
        return m;
    }
    /**
     * Get a matrix reduced by one column and one row.
     *
     * @code
     * Matrix<3,3,int> m(1,2,3, 4,5,6, 7,8,9);  // [(1,2,3), (4,5,6), (7,8,9)]
     * m.GetReduced()                           // [(1,2), (4,5)]
     * @endcode
     *
     * @return Reduced matrix.
     */
    Matrix<3,3,float> GetReduced() {
        Matrix<3,3,float> m;
        for (unsigned int i=0; i < 3; i++)
            for (unsigned int j=0; j < 3; j++)
                m(i,j) = elm[i].get(j);
        return m;
    }
    /**
     * Create array of matrix.
     *
     * @param a Array to populate
     */
    void ToArray(float *a) const
    {
        for(int x=0; x<4; x++)
            elm[x].ToArray<4>(a+(x*4));
    }
    /**
     * String representation.
     * Ex. [(1, 2), (3, 4)]
     *
     * @return Matrix as string
     */
    std::string ToString() const
    {
        std::ostringstream out;
        out << "[";
        for (unsigned int i=0; i<3; i++) {
            out << "(";
            for (unsigned int j=0; j<3; j++)
                out << elm[i].get(j) << ", ";
            out << elm[i].get(3) << "), ";
        }
        out << "(";
        for (unsigned int j=0; j<3; j++)
            out << elm[3].get(j) << ", ";
        out << elm[3].get(3) << ")]";
        return out.str();
    }

}; // Matrix
}  // NS OpenEngine
}  // NS Math

#endif // __SSE__
#endif // _OE_FLOAT4_INCLUDED
#endif // _OE_MATRIX44_H_
