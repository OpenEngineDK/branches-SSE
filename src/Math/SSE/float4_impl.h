#ifndef _OE_FLOAT4_IMPL_H
#define _OE_FLOAT4_IMPL_H

//GET EPS
#include <Math/Math.h>

// Streaming
#include <iostream>

namespace OpenEngine {
namespace Math {
namespace SSE {
    //------------//
    //Constructors//
    //------------//
    float4::float4()
    {
        Data = ZeroSetAll();
    }

    float4::float4(const float f)
    {
        Data = FloatSetAll(f);
    }
    
    float4::float4(const float *a)
    {
        Data = FloatSetEachArray(a);
    }
    
    float4::float4(const float x, const float y, const float z, const float w)
    {
        Data = FloatSetEach(x, y, z, w);
    }
    
    float4::float4(const float4 &m)
    {
        for (unsigned int i=0;i<4;i++)
            DataArray[i] = m.DataArray[i];
    }
    
    float4::float4(const SSE_float4 m)
    : Data(m)
    {
    }
    
    //-----------//
    //Destructors//
    //-----------//
    float4::~float4()
    {
    }
    
    //-----//
    //Debug//
    //-----//
    void float4::PrintString() const
    {
      for(int x=0; x<3; x++)
	     ::std::cout << DataArray[x] << ", ";
      ::std::cout << DataArray[3];
    }
    
    template<unsigned N>
    std::string float4::ToString() const
    {
        std::ostringstream out;
        out << "(";
        for (unsigned int i=0; i<N-1; i++)
        {
            out << DataArray[i] << ", ";
        }
        out << DataArray[N-1] << ")";
        return out.str();
    }
    
    SSE_float4 float4::getData() const
    {
        return Data;
    }
    
    template<unsigned N>
    void float4::ToArray(float *p) const
    {
        for(unsigned int x=0; x<N; x++)
            p[x] = DataArray[x];        
    }
    
    void float4::ToArray(float& x, float& y, float& z, float& w) const
    {
        x = _x;
        y = _y;
        z = _z;
        w = _w;
    }
        
    //---------//
    //Operators//
    //---------//
    float& float4::operator[](const unsigned int i)
    {
        return DataArray[i];
        //return reinterpret_cast<float*>(&Data)[i];
    }
    
    float float4::get(const unsigned int i) const
    {
        return DataArray[i];
    }
    
    //+ = + to each element (and for scalar)
    float4 float4::operator+(const float4 m) const
    {
        return Add(Data, m.Data);
    }

    void float4::operator+=(const float4 m)
    {    
        Data = Add(Data, m.Data);
    }

    //- = - to each element (and for scalar)
    float4 float4::operator-(const float4 m) const
    {
        return Sub(Data, m.Data);
    }
    
    void float4::operator-=(const float4 m)
    {
        Data = Sub(Data, m.Data);
    }

    //* = * to each element (and for scalar)
    float4 float4::operator*(const float4 m) const
    {
        return Mul(Data, m.Data);
    }

    void float4::operator*=(const float4 m)
    {
        Data = Mul(Data, m.Data);
    }

    // (/) = (/) to each element (and for scalar)
    float4 float4::operator/(const float4 m) const
    {
        return Div(Data, m.Data);
    }
    
    void float4::operator/=(const float4 m)
    {
        Data = Div(Data, m.Data);
    }

    // (==) = (==) brug epsilon EPS fra OpenEngine
    float4 float4::operator==(const float4 m) const
    {
        return (float4(Data-m.Data) < float4(EPS).Data);
    }
    
    float4 float4::ExactlyEqual(const float4 m) const
    {
        return Equal(Data, m.Data);
    }
    
    float4 float4::operator!=(const float4 m) const
    {
        return (float4(Data-m.Data) > float4(EPS).Data);
    }
    
    float4 float4::ExactlyInEqual(const float4 m) const
    {
        return InEqual(Data, m.Data);
    }
    
    // (>) = (>)
    float4 float4::operator>(const float4 m) const
    {
        return GreaterThan(Data, m.Data);
    }
    
    // (<) = (<)
    float4 float4::operator<(const float4 m) const
    {
        return LessThan(Data, m.Data);
    }
    
    // (>=) = (>=)
    float4 float4::operator>=(const float4 m) const
    {
        return GreaterThanOrEqual(Data, m.Data);
    }
    
    // (<=) = (<=)
    float4 float4::operator<=(const float4 m) const
    {
        return LessThanOrEqual(Data, m.Data);
    }
    
    // (=) = (=)
    float4& float4::operator=(const float4& m)
    {
        Data = m.Data;
        return *this;
    }

    // (^) = (^)
    float4 float4::operator^(const float4 m) const
    {
        return XOR(Data, m.Data);
    }
    
    void float4::operator^=(const float4 m)
    {
        Data = XOR(Data, m.Data);
    }
    
    // (&) = (&)
    float4 float4::operator&(const float4 m) const
    {
        return AND(Data, m.Data);
    }
    
    void float4::operator&=(const float4 m)
    {
        Data = AND(Data, m.Data);
    }
    
    // (|) = (|)
    float4 float4::operator|(const float4 m) const
    {
        return OR(Data, m.Data);
    }
    

    void float4::operator|=(const float4 m)
    {
        Data = OR(Data, m.Data);
    }
    
    //(>>) = (>>)
    float4 float4::operator>>(const int b) const
    {
        SSE_float4 Temp = Data;
        for(int x=0; x<b; x++)
            Temp = Float_Shift_Right(Temp);
        return Temp;
    }
    
    void float4::operator>>=(const int b)
    {
        for(int x=0; x<b; x++)
            Data = Float_Shift_Right(Data);
    }

    //(<<) = (<<)
    float4 float4::operator<<(const int b) const
    {
        SSE_float4 Temp = Data;
        for(int x=0; x<b; x++)
            Temp = Float_Shift_Left(Temp);
        return Temp;
    }
    
    void float4::operator<<=(const int b)
    {
        for(int x=0; x<b; x++)
            Data = Float_Shift_Right(Data);
    }

    float4 float4::operator-() const
    {
        return (Data * float4((float)-1.0f).getData());
    }
    
    #ifdef _OE_FLOAT4_FLOAT_OPERATORS
        float4 float4::operator+(const float scalar) const
        {
            float4 Temp(scalar);
            return Add(Data, Temp.Data);
        }
    
        void float4::operator+=(const float scalar)
        {
            float4 Temp(scalar);
            Data = Add(Data, Temp.Data);
        }
    
        float4 float4::operator-(const float scalar) const
        {
            float4 Temp(scalar);
            return Sub(Data, Temp.Data);
        }
    
        void float4::operator-=(const float scalar)
        {
            float4 Temp(scalar);
            Data = Sub(Data, Temp.Data);
        }
    
        float4 float4::operator*(const float scalar) const
        {
            float4 Temp(scalar);
            return Mul(Data, Temp.Data);
        }
    
        void float4::operator*=(const float scalar)
        {
            float4 Temp(scalar);
            Data = Mul(Data, Temp.Data);
        }
    
        float4 float4::operator/(const float scalar) const
        {
            float4 Temp(scalar);
            return Div(Data, Temp.Data);
        }
        
        void float4::operator/=(const float scalar)
        {
            float4 Temp(scalar);
            Data = Div(Data, Temp.Data);
        }
    
    
        float4& float4::operator=(const float m)
        {
            Data = FloatSetAll(m);
        }
    
        
        float4 float4::operator^(const float m) const
        {
            float4 Temp(m);
            return XOR(Data, Temp.Data);
        }
    
    
        void float4::operator^=(const float m)
        {
            float4 Temp(m);
            Data = XOR(Data, Temp.Data);
        }
    
        float4 float4::operator&(const float m) const
        {
            float4 Temp(m);
            return AND(Data, Temp.Data);
        }
    
        void float4::operator&=(const float m)
        {
            float4 Temp(m);
            Data = AND(Data, Temp.Data);
        }
    
        float4 float4::operator|(const float m) const
        {
            float4 Temp(m);
            return OR(Data, Temp.Data);
        }
        
        void float4::operator|=(const float m)
        {
            float4 Temp(m);
            Data = OR(Data, Temp.Data);
        } 
        
        bool float4::operator==(const float m) const
        {
            float4 Temp(m);
            if((Data-Temp.Data)>SSE_EPS)
               return true;
            else
               return false;            
        }

        bool float4::isExactly(const float m) const
        {
            float4 Temp(m);
            return Equal(Data, Temp.Data); 
        }
        
        bool float4::operator!=(const float m) const
        {
            float4 Temp(m);
            return !(Data == Temp.Data);
        }
        
        float4 float4::ExactlyInEqual(const float m) const
        {
            float4 Temp(m);
            return InEqual(Data, Temp.Data); 
        }
        
    #endif //_OE_FLOAT4_FLOAT_OPERATORS
    //----------------//
    //Internal Workers//
    //----------------//
    //Set//
    //---//
    /**
     * ZeroSet()
     * @brief "Clears the four single-precision, floating-point values."
     * @brief returns an SSE_float4 in the form of 0,0,0,0.
     * @return SSE_float4(0x0, 0x0, 0x0, 0x0) 
     */
    SSE_float4 float4::ZeroSetAll()
    {
        return _mm_setzero_ps();
    }
    /**
     * FloatSet(float f)
     * @brief "Sets the four single-precision, floating-point values to f."
     * @brief returns an SSE_float4 in the form of f,f,f,f.
     * @param f to be stored in SSE_float4
     * @return SSE_float4(f, f, f, f) 
     */
    SSE_float4 float4::FloatSetAll(float f)
    {
        return _mm_set1_ps(f);
    }
    /**
     * FloatSet(float x, float y, float z, float w)
     * @brief "Sets the four single-precision, floating-point values to the four inputs."
     * @brief returns an SSE_float4 in the form of x,y,z,w.
     * @param x to be stored in SSE_float4(1)
     * @param y to be stored in SSE_float4(2) (default val = 0)
     * @param z to be stored in SSE_float4(3) (default val = 0)
     * @param w to be stored in SSE_float4(4) (default val = 0)
     * @return SSE_float4(x, y, z, f) 
     */
    SSE_float4 float4::FloatSetEach(float x, float y, float z, float w)
    {
        return _mm_setr_ps(x, y, z, w);
    }
    /**
     * FloatSet(const float *p)
     * @brief "Loads four single-precision, floating-point values."
     * @brief returns an SSE_float4 in the form of p[0], p[1], p[2], p[3]
     * @param p[0] to be stored in SSE_float4(1)
     * @param p[1] to be stored in SSE_float4(2) (default val = 0)
     * @param p[2] to be stored in SSE_float4(3) (default val = 0)
     * @param p[3] to be stored in SSE_float4(4) (default val = 0)
     * @return SSE_float4(p[0], p[1], p[2], p[3]) 
     */    
    SSE_float4 float4::FloatSetEachArray(const float *p)
    {
        return _mm_load_ps(p);
    }
    
    void float4::Store(float *p, SSE_float4 m)
    {
        _mm_store_ps((float*)p, m);
    }
    
    /**
     * Add(SSE_float4 a, SSE_float4 b)
     * @breif "Adds the four single-precision, floating-point values of a and b."
     * @param a SSE_float4
     * @param b SSE_float4     
     * @return SSE_float4(a0+b0, a1+b1, a2+b2, a3+b3) 
     */
    SSE_float4 float4::Add(SSE_float4 a, SSE_float4 b)
    {
        return _mm_add_ps(a, b);
    }

    /**
     * Sub(SSE_float4 a, SSE_float4 b)
     * @breif "Subtracts the four single-precision, floating-point values of a and b."
     * @param a SSE_float4
     * @param b SSE_float4     
     * @return SSE_float4(a0-b0, a1-b1, a2-b2, a3-b3) 
     */
    SSE_float4 float4::Sub(SSE_float4 a, SSE_float4 b)
    {
        return _mm_sub_ps(a, b);
    }

    /**
     * Mul(SSE_float4 a, SSE_float4 b)
     * @breif "Multiplies the four single-precision, floating-point values of a and b."
     * @param a SSE_float4
     * @param b SSE_float4     
     * @return SSE_float4(a0*b0, a1*b1, a2*b2, a3*b3) 
     */
    SSE_float4 float4::Mul(SSE_float4 a, SSE_float4 b)
    {
        return _mm_mul_ps(a, b);
    }

    /**
     * Div(SSE_float4 a, SSE_float4 b)
     * @breif "Divides the four single-precision, floating-point values of a and b."
     * @param a SSE_float4
     * @param b SSE_float4     
     * @return SSE_float4(a0/b0, a1/b1, a2/b2, a3/b3) 
     */
    SSE_float4 float4::Div(SSE_float4 a, SSE_float4 b)
    {
        return _mm_div_ps(a, b);
    }

    /**
     * AND(SSE_float4 a, SSE_float4 b)
     * @breif "Computes the bitwise AND of the four single-precision, floating-point values of a and b."
     * @param a SSE_float4  
     * @param b SSE_float4 
     * @return SSE_float4(a0 & b0, a1 & b1, a2 & b2, a3 & b3) 
     */    
    SSE_float4 float4::AND(SSE_float4 a, SSE_float4 b)
    {
        return _mm_and_ps(a, b);
    }    

    /**
     * NOT(SSE_float4 a, SSE_float4 b)
     * @breif "Computes the bitwise AND-NOT of the four single-precision, floating-point values of a and b."
     * @param a SSE_float4  
     * @param b SSE_float4 
     * @return SSE_float4(~a0 & b0, ~a1 & b1, ~a2 & b2, ~a3 & b3) 
     */    
    SSE_float4 float4::NOT(SSE_float4 a, SSE_float4 b)
    {
        return _mm_andnot_ps(a, b);
    }

    /**
     * OR(SSE_float4 a, SSE_float4 b)
     * @breif "Computes the bitwise OR of the four single-precision, floating-point values of a and b."
     * @param a SSE_float4  
     * @param b SSE_float4 
     * @return SSE_float4(a0 | b0, a1 | b1, a2 | b2, a3 | b3) 
     */   
    SSE_float4 float4::OR(SSE_float4 a, SSE_float4 b)
    {
        return _mm_or_ps(a, b);
    }

    /**
     * XOR(SSE_float4 a, SSE_float4 b)
     * @breif "Computes bitwise EXOR (exclusive-or) of the four single-precision, floating-point values of a and b."
     * @param a SSE_float4  
     * @param b SSE_float4 
     * @return SSE_float4(a0 ^ b0, a1 ^ b1, a2 ^ b2, a3 ^ b3) 
     */  
    SSE_float4 float4::XOR(SSE_float4 a, SSE_float4 b)
    {
        return _mm_xor_ps(a, b);
    }   

    /**
     * Float_Shift_Flip(SSE_float4 a)
     * @breif "Flips the SSE_float4 data"
     * @param a SSE_float4  
     * @return SSE_float4(a3, a2, a1, a0) 
     */     
    SSE_float4 float4::Float_Shift_Flip(SSE_float4 a)
    {
        return _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 1, 2, 3));
    }
 
     /**
     * Float_Shift_Left(SSE_float4 a)
     * @breif "Shifts the SSE_float4 data left"
     * @param a SSE_float4  
     * @return SSE_float4(a1, a2, a3, a0) 
     */       
    SSE_float4 float4::Float_Shift_Left(SSE_float4 a)
    {
        return _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1));
    }

     /**
     * Float_Shift_Right(SSE_float4 a)
     * @breif "Shifts the SSE_float4 data right"
     * @param a SSE_float4  
     * @return SSE_float4(a1, a2, a3, a0) 
     */           
    SSE_float4 float4::Float_Shift_Right(SSE_float4 a)
    {
        return _mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 1, 0, 3));
    }

    SSE_float4 float4::Equal(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmpeq_ps(a, b);
    }
    
    SSE_float4 float4::LessThan(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmplt_ps(a, b);
    }
    
    SSE_float4 float4::LessThanOrEqual(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmple_ps(a, b);
    }
    
    SSE_float4 float4::GreaterThan(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmpgt_ps(a, b);
    }
    
    SSE_float4 float4::GreaterThanOrEqual(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmpge_ps(a, b);
    }
    
    SSE_float4 float4::InEqual(SSE_float4 a, SSE_float4 b)
    {
        return _mm_cmpneq_ps(a, b);
    }

}  // NS OpenEngine
}  // NS Math
}  // NS SSE

#endif //_OE_FLOAT4_IMPL_H
