#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x1vMax_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
  tmp[0] = 0.5*(A[9]*B[9]+A[8]*B[8]+A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.004761904761904762*(92.22255689363637*(A[4]*B[8]+B[4]*A[8])+105.0*(A[5]*B[7]+B[5]*A[7])+93.91485505499116*(A[3]*B[6]+B[3]*A[6]+A[1]*B[4]+B[1]*A[4])+105.0*(A[2]*B[3]+B[2]*A[3]+A[0]*B[1]+B[0]*A[1])); 
  tmp[2] = 0.004761904761904762*(92.22255689363637*(A[5]*B[9]+B[5]*A[9])+93.91485505499116*(A[3]*B[7]+B[3]*A[7])+105.0*(A[4]*B[6]+B[4]*A[6])+93.91485505499116*(A[2]*B[5]+B[2]*A[5])+105.0*(A[1]*B[3]+B[1]*A[3]+A[0]*B[2]+B[0]*A[2])); 
  tmp[3] = 0.004761904761904762*(92.2225568936364*(A[7]*B[9]+B[7]*A[9]+A[6]*B[8]+B[6]*A[8])+(84.0*A[6]+93.91485505499116*A[2])*B[7]+(84.0*B[6]+93.91485505499116*B[2])*A[7]+93.91485505499116*(A[1]*B[6]+B[1]*A[6]+A[3]*B[5]+B[3]*A[5]+A[3]*B[4]+B[3]*A[4])+105.0*(A[0]*B[3]+B[0]*A[3]+A[1]*B[2]+B[1]*A[2])); 
  tmp[4] = 0.004761904761904762*((62.60990336999411*A[8]+92.22255689363637*A[1])*B[8]+92.22255689363637*B[1]*A[8]+93.91485505499116*A[7]*B[7]+(67.0820393249937*A[6]+105.0*A[2])*B[6]+105.0*B[2]*A[6]+(67.0820393249937*A[4]+105.0*A[0])*B[4]+105.0*B[0]*A[4]+93.91485505499116*(A[3]*B[3]+A[1]*B[1])); 
  tmp[5] = 0.004761904761904762*((62.60990336999411*A[9]+92.22255689363637*A[2])*B[9]+92.22255689363637*B[2]*A[9]+(67.0820393249937*A[7]+105.0*A[1])*B[7]+105.0*B[1]*A[7]+93.91485505499116*A[6]*B[6]+(67.0820393249937*A[5]+105.0*A[0])*B[5]+105.0*B[0]*A[5]+93.91485505499116*(A[3]*B[3]+A[2]*B[2])); 
  tmp[6] = 0.004761904761904762*(92.2225568936364*(A[3]*B[8]+B[3]*A[8])+84.0*(A[3]*B[7]+B[3]*A[7])+(93.91485505499116*A[5]+67.0820393249937*A[4]+105.0*A[0])*B[6]+(93.91485505499116*B[5]+67.0820393249937*B[4]+105.0*B[0])*A[6]+105.0*(A[2]*B[4]+B[2]*A[4])+93.91485505499116*(A[1]*B[3]+B[1]*A[3])); 
  tmp[7] = 0.004761904761904762*(92.2225568936364*(A[3]*B[9]+B[3]*A[9])+(67.0820393249937*A[5]+93.91485505499116*A[4]+105.0*A[0])*B[7]+(67.0820393249937*B[5]+93.91485505499116*B[4]+105.0*B[0])*A[7]+84.0*(A[3]*B[6]+B[3]*A[6])+105.0*(A[1]*B[5]+B[1]*A[5])+93.91485505499116*(A[2]*B[3]+B[2]*A[3])); 
  tmp[8] = 0.004761904761904762*((62.60990336999411*A[4]+105.0*A[0])*B[8]+(62.60990336999411*B[4]+105.0*B[0])*A[8]+92.2225568936364*(A[3]*B[6]+B[3]*A[6])+92.22255689363637*(A[1]*B[4]+B[1]*A[4])); 
  tmp[9] = 0.004761904761904762*((62.60990336999411*A[5]+105.0*A[0])*B[9]+(62.60990336999411*B[5]+105.0*B[0])*A[9]+92.2225568936364*(A[3]*B[7]+B[3]*A[7])+92.22255689363637*(A[2]*B[5]+B[2]*A[5])); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x1vMax_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
  tmp[0] = 0.7071067811865476*A[3]*B[8]+0.7071067811865476*A[2]*B[4]+0.7071067811865476*A[1]*B[1]+0.7071067811865476*A[0]*B[0]; 
  tmp[1] = 0.6210590034081187*A[2]*B[8]+0.6210590034081187*A[3]*B[4]+0.632455532033676*A[1]*B[4]+0.632455532033676*B[1]*A[2]+0.7071067811865476*A[0]*B[1]+0.7071067811865476*B[0]*A[1]; 
  tmp[2] = 0.7071067811865477*A[2]*B[6]+0.7071067811865476*A[1]*B[3]+0.7071067811865476*A[0]*B[2]; 
  tmp[3] = 0.621059003408119*A[3]*B[6]+0.632455532033676*A[1]*B[6]+0.632455532033676*A[2]*B[3]+0.7071067811865476*A[0]*B[3]+0.7071067811865476*A[1]*B[2]; 
  tmp[4] = 0.4216370213557841*A[3]*B[8]+0.6210590034081188*A[1]*B[8]+0.4517539514526258*A[2]*B[4]+0.7071067811865476*A[0]*B[4]+0.6210590034081188*B[1]*A[3]+0.7071067811865476*B[0]*A[2]+0.632455532033676*A[1]*B[1]; 
  tmp[5] = 0.7071067811865477*A[1]*B[7]+0.7071067811865476*A[0]*B[5]; 
  tmp[6] = 0.4517539514526258*A[2]*B[6]+0.7071067811865476*A[0]*B[6]+0.621059003408119*A[3]*B[3]+0.632455532033676*A[1]*B[3]+0.7071067811865476*A[2]*B[2]; 
  tmp[7] = 0.632455532033676*A[2]*B[7]+0.7071067811865476*A[0]*B[7]+0.7071067811865477*A[1]*B[5]; 
  tmp[8] = 0.4216370213557841*A[2]*B[8]+0.7071067811865476*A[0]*B[8]+0.4216370213557841*A[3]*B[4]+0.6210590034081188*A[1]*B[4]+0.7071067811865476*B[0]*A[3]+0.6210590034081188*B[1]*A[2]; 
  tmp[9] = 0.7071067811865475*A[0]*B[9]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x1vMax_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.870828693386971*A[3])+1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.870828693386971*A[3]+1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
 
  double As[4]; 
  double Bs[10]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = 0.0; 
    Bs[4] = 0.0; 
    Bs[5] = B[5]; 
    Bs[6] = 0.0; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = B[9]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = B[7]; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(10,10); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,4) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,5) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(1,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(1,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(1,4) = 0.7071067811865475*As[1]; 
  data->AEM_D(1,5) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(1,6) = 0.7071067811865475*As[2]; 
  data->AEM_D(1,7) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  data->AEM_D(2,6) = 0.7071067811865475*As[2]; 
  data->AEM_D(2,7) = 0.6210590034081187*As[3]+0.632455532033676*As[1]; 
  data->AEM_D(3,2) = 0.7071067811865475*As[3]; 
  data->AEM_D(3,3) = 0.6210590034081186*As[2]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,10,1) = data->u_D; 
 
} 
 
