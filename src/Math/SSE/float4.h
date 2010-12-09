#ifndef _OE_FLOAT4_H
#define _OE_FLOAT4_H

#define SSE_NEEDED

// Include only if SSE is enabled (otherwise an error is thrown, then this is included).
// Possibly use CPUID()
#ifdef __SSE__
#define _OE_FLOAT4_INCLUDED

//For ToString
#include <string>
#include <sstream>

//Macros for FORCE_INLINE and ALIGNED(x)
#include "Gcc_Specific_Macros.h"

//Enable the float operators
#ifdef _OE_FLOAT4_FLOAT_OPERATORS
    #warning NOTE: "_OE_FLOAT4_FLOAT_OPERATORS" declared, float operators unlocked!
    #warning NOTE: ^^Its bad pratice to use floats directly with float4 (memory latencies)^^
#endif //_OE_FLOAT4_FLOAT_OPERATORS

// Header for SSE
#include <xmmintrin.h> 

namespace OpenEngine {
namespace Math {
namespace SSE {
    
    //Definintion (__m128 is the SSE registers type)
    typedef __m128 SSE_float4;

    //float4 class (wrapper around __m128, providing overloading, ect).
    //ALIGNED(16) class float4 //No need to align since __m128 forces it, and as thats the only member
    ALIGNED(x) class float4
    {
        //---------
        //Functions
        public:
            //Constructors
            FORCE_INLINE float4(); //Overvej at lave en tom default constructor (istedet for samtlige elementer = 0)
            #ifdef _OE_FLOAT4_FLOAT_OPERATORS
                FORCE_INLINE explicit float4(const float f);
                FORCE_INLINE explicit float4(const float *a);
                FORCE_INLINE explicit float4(const float x, const float y, const float z, const float w = 0.0f);
            #else
                FORCE_INLINE float4(const float f);
                FORCE_INLINE float4(const float *a);
                FORCE_INLINE float4(const float x, const float y, const float z, const float w = 0.0f);
            #endif  //_OE_FLOAT4_FLOAT_OPERATORS
            FORCE_INLINE float4(const float4 &m);
            FORCE_INLINE float4(const SSE_float4 m);
            FORCE_INLINE ~float4();
            
            //Debugging
            FORCE_INLINE void PrintString() const;
            template<unsigned N>
            FORCE_INLINE std::string ToString() const;
           
            //Getting Dataz
            FORCE_INLINE SSE_float4 getData() const;
            template<unsigned N>
            FORCE_INLINE void ToArray(float *p) const;
            FORCE_INLINE void ToArray(float& x, float& y, float& z, float& w) const;
                        
            //Operators (float scalar, avoid this, as its slow)
            //[] for access ((assert for values > 4)
            FORCE_INLINE float& operator[](const unsigned int i);
            FORCE_INLINE float get(const unsigned int i) const;
            //+ = + to each element (and for scalar)
            FORCE_INLINE float4 operator+(const float4 m) const;
            FORCE_INLINE void operator+=(const float4 m);
            //- = - to each element (and for scalar)
            FORCE_INLINE float4 operator-(const float4 m) const;
            FORCE_INLINE void operator-=(const float4 m);
            //* = * to each element (and for scalar)
            FORCE_INLINE float4 operator*(const float4 m) const;
            FORCE_INLINE void operator*=(const float4 m);
            // (/) = (/) to each element (and for scalar)
            FORCE_INLINE float4 operator/(const float4 m) const;
            FORCE_INLINE void operator/=(const float4 m);
            // (==) = (==) brug epsilon EPS fra OpenEngine
            FORCE_INLINE float4 operator==(const float4 m) const;
            FORCE_INLINE float4 ExactlyEqual(const float4 m) const;
            // (!=) = (!=)
            FORCE_INLINE float4 operator!=(const float4 m) const;
            FORCE_INLINE float4 ExactlyInEqual(const float4 m) const;
            // (>) = (>)
            FORCE_INLINE float4 operator>(const float4 m) const;
            // (<) = (<)
            FORCE_INLINE float4 operator<(const float4 m) const;
            // (>=) = (>=)
            FORCE_INLINE float4 operator>=(const float4 m) const;
            // (<=) = (<=)
            FORCE_INLINE float4 operator<=(const float4 m) const;
            // (!) = (!)
            //FORCE_INLINE float4 operator!() const;
            // (=) = (=)
            FORCE_INLINE float4& operator=(const float4& m);
            // (^) = (^)
            FORCE_INLINE float4 operator^(const float4 m) const;
            FORCE_INLINE void operator^=(const float4 m);
            // (&) = (&)
            FORCE_INLINE float4 operator&(const float4 m) const;
            FORCE_INLINE void operator&=(const float4 m);
            // (|) = (|)
            FORCE_INLINE float4 operator|(const float4 m) const;
            FORCE_INLINE void operator|=(const float4 m);
            // (~) = (~)
            //FORCE_INLINE float4 operator~();     
            //(>>) = (>>)
            FORCE_INLINE float4 operator>>(const int b) const;
            FORCE_INLINE void operator>>=(const int b);
            //(<<) = (<<)
            FORCE_INLINE float4 operator<<(const int b) const;
            FORCE_INLINE void operator<<=(const int b);
            // (-) = (-)
            FORCE_INLINE float4 operator-() const;
                        
            //if floating operators are included, get these too
            #ifdef _OE_FLOAT4_FLOAT_OPERATORS
                FORCE_INLINE float4 operator+(const float scalar) const;
                FORCE_INLINE void operator+=(const float scalar);
                FORCE_INLINE float4 operator-(const float scalar) const;
                FORCE_INLINE void operator-=(const float scalar);
                FORCE_INLINE float4 operator*(const float scalar) const;
                FORCE_INLINE void operator*=(const float scalar);
                FORCE_INLINE float4 operator/(const float scalar) const;
                FORCE_INLINE void operator/=(const float scalar);  
                FORCE_INLINE float4& operator=(const float m);
                FORCE_INLINE float4 operator^(const float m) const;
                FORCE_INLINE void operator^=(const float m);
                FORCE_INLINE float4 operator&(const float m) const;
                FORCE_INLINE void operator&=(const float m);   
                FORCE_INLINE float4 operator|(const float m) const;
                FORCE_INLINE void operator|=(const float m);  
                FORCE_INLINE float4 operator==(const float m) const;
                FORCE_INLINE float4 ExactlyEqual(const float m) const;
                FORCE_INLINE float4 operator!=(const float m) const; 
                FORCE_INLINE float4 ExactlyInEqual(const float m) const;  
            #endif //_OE_FLOAT4_FLOAT_OPERATORS
            
        private:
        //Internal Workers
            //Set
            FORCE_INLINE static SSE_float4 ZeroSetAll();
            FORCE_INLINE static SSE_float4 FloatSetAll(float f);
            FORCE_INLINE static SSE_float4 FloatSetEach(float x, float y, float z, float w);
            FORCE_INLINE static SSE_float4 FloatSetEachArray(const float *p);
            FORCE_INLINE static void Store(float *p, SSE_float4 m);
            //Arithmetic Operations
            FORCE_INLINE static SSE_float4 Add(SSE_float4 a, SSE_float4 b); //Adds
            FORCE_INLINE static SSE_float4 Sub(SSE_float4 a, SSE_float4 b); //Subtracts
            FORCE_INLINE static SSE_float4 Mul(SSE_float4 a, SSE_float4 b); //Multiplies
            FORCE_INLINE static SSE_float4 Div(SSE_float4 a, SSE_float4 b); //Divides
            //Logical Operations
            FORCE_INLINE static SSE_float4 AND(SSE_float4 a, SSE_float4 b); //Bitwise AND
            FORCE_INLINE static SSE_float4 NOT(SSE_float4 a, SSE_float4 b); //Logical NOT
            FORCE_INLINE static SSE_float4 OR(SSE_float4 a, SSE_float4 b);  //Bitwise OR
            FORCE_INLINE static SSE_float4 XOR(SSE_float4 a, SSE_float4 b); //Bitwise Exclusive OR
            //Shuffle
            FORCE_INLINE static SSE_float4 Float_Shift_Flip(SSE_float4 a);
            FORCE_INLINE static SSE_float4 Float_Shift_Left(SSE_float4 a);
            FORCE_INLINE static SSE_float4 Float_Shift_Right(SSE_float4 a);
            //Comparison
            FORCE_INLINE static SSE_float4 Equal(SSE_float4 a, SSE_float4 b);
            FORCE_INLINE static SSE_float4 LessThan(SSE_float4 a, SSE_float4 b);
            FORCE_INLINE static SSE_float4 LessThanOrEqual(SSE_float4 a, SSE_float4 b);
            FORCE_INLINE static SSE_float4 GreaterThan(SSE_float4 a, SSE_float4 b);
            FORCE_INLINE static SSE_float4 GreaterThanOrEqual(SSE_float4 a, SSE_float4 b);
            FORCE_INLINE static SSE_float4 InEqual(SSE_float4 a, SSE_float4 b);
        protected:
        //---------
        //Variables
        public:
        private:
            // Data Union
            ALIGNED(16) union
            {
                ALIGNED(16) SSE_float4 Data;
                float DataArray[4];
                struct
                {
                    float _x, _y, _z, _w;
                };
            };
        protected:    
    };
    
}  // NS OpenEngine
}  // NS Math
}  // NS SSE

//class implementation
#include "float4_impl.h"

//to overload << and >> for iostream
#ifdef _OE_FLOAT4_IOSTREAM_OVERLOAD
    #include "float4_iostream.h"
#endif //_OE_FLOAT4_IOSTREAM_OVERLOAD

#endif //__SSE__

//Check SSE status based upon defines
#include "SSE_Status.h"

#endif //_OE_FLOAT4_H
