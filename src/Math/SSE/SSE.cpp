#include <Math/SSE/SSE.h>
#include <Math/SSE/float4.h>

// Need EPS
#include <Math/Math.h>
extern float EPS;

namespace OpenEngine {
namespace Math {
namespace SSE {

 //for == pyrecision checks
float4 SSE_EPS(EPS); //extern EPS is defined in Math\\Math.h
    
}  // NS OpenEngine
}  // NS Math
}  // NS SSE
