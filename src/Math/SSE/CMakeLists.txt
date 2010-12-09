ADD_LIBRARY(OpenEngine_Math_SSE
  float4.h
  float4_impl.h
  Gcc_Specific_Macros.h
  SSE_Status.h
  SSE.cpp
  SSE.h
  FloatMatrix.h
  Matrix44.h
  FloatVector.h
  Vector4.h
  Vector3.h
  SSEMath.h
)

TARGET_LINK_LIBRARIES(OpenEngine_Math
)
