#include <iostream>
#include <fstream>
#include<math.h>
#include <Rcpp.h>
using namespace Rcpp; 

#define dim 25
#define dimIv 55
#define dimOut 105

//Note: All pieces of code beginning with a @ will be replaced by the required code by R before compiling
//For instance @AddDim will be replaced by the dimension of the model

// utility function that converts a Rcpp::List to a double**
// WARNING: do not forget to free the double** after use!
template<typename T>
void RcppListToPptr(Rcpp::List L, T**& pptr) {
	for (unsigned int it=0; it<L.size(); it++) {
		std::vector<double> tempVec = L[it];
		pptr[it] = (double*) malloc(sizeof(*pptr[it]) * tempVec.size());
		for (unsigned int it2=0; it2<tempVec.size(); it2++) {
			pptr[it][it2] = tempVec[it2];
		}
	}
}

void Func(double t, double* y, double* parms, double* ydot, double* x, double** dataExogVar, double** exogSamplingTime, int nExogVar, int* comptExogVar) {

ydot[0] = parms[0] * y[0];
ydot[5] = y[6] * y[5];
ydot[6] = parms[48] * y[6];
ydot[8] = parms[7] * y[8];
ydot[9] = parms[12] * y[9];
ydot[10] = parms[49] * y[10];
ydot[12] = parms[46] * y[11] - (parms[46] * parms[43]/parms[44] + parms[47]) * y[12] + parms[47] * y[13] * parms[44]/parms[45];
ydot[13] = parms[47] * y[12] - parms[47] * y[13] * parms[44]/parms[45];
ydot[15] = parms[42] * (y[14] - y[15])/parms[37];
x[4] = parms[2] * (1.0 - y[3]/parms[9]);
x[10] = y[4]/parms[5];
x[20] = std::min(pow((y[9]/y[8]), (1.0/(parms[35] - 1.0))), 1.0);
x[22] = 0.0;
x[27] = -parms[68] * (parms[65] * y[17] * (1.0 - y[17]/parms[63]));
x[28] = -parms[69] * (parms[66] * y[18] * (1.0 - y[18]/parms[64]));
x[33] = (parms[54]/log(2.0)) * log(y[11]/parms[43]);
x[34] = 0.5 + t/(parms[13] * 2.0);
x[36] = 5.0/(1.0/parms[36] + parms[41] * (parms[39] - parms[40]));
x[37] = std::min(y[9], y[8]);
x[40] = 100.0 * y[19]/388.0;
x[45] = 0.0;
x[51] = y[20] + y[21] + y[22];
x[52] = log(2.0)/parms[51];
x[53] = parms[75] * (y[20] * parms[68] + y[21] * parms[69] + y[22] * parms[70]);
x[54] = 0.0;
ydot[3] = y[3] * x[4];
ydot[23] = -parms[52] * y[23] + x[53];
x[5] = x[10]/y[0];
x[19] = (y[5] * y[8]/(1000.0 * parms[35])) * pow(x[20], parms[35]);
x[23] = parms[34] * x[22];
x[26] = y[23] * parms[52] + (1.0 - parms[75]) * (y[20] * parms[68] + y[21] * parms[69] + y[22] * parms[70]);
x[35] = x[33] + x[34];
x[38] = parms[50] * x[10] * y[5] * (1.0 - x[20]);
x[46] = (1.0 - x[45]) * y[10]/parms[70];
x[47] = y[10] * x[45];
ydot[14] = x[35]/x[36] - parms[54]/(x[36] * parms[39]) * y[14] - (parms[42]/x[36]) * (y[14] - y[15]);
ydot[16] = -y[16] * x[46]/y[19];
x[6] = x[5]/y[3];
x[9] = parms[4] + x[23];
x[24] = 1.0 - (1.0 - x[22])/(1.0 - x[23]);
x[29] = -parms[70] * (parms[67] * y[19] * (1.0 - y[19]/y[16]) - x[46]);
x[39] = std::min(pow(((parms[78] + y[9])/parms[79] * (1.0 - parms[76])/parms[76]), (parms[76] - 1.0)) * x[38]/parms[77], y[24]);
ydot[24] = -x[39];
x[3] = parms[20] + parms[21] * x[6];
x[11] = x[10] * (1.0 - x[24]) * (1.0 - x[19]);
x[21] = x[37] * parms[14] * x[39];
x[30] = x[27] + x[28] + x[29];
x[41] = pow((x[38]/parms[77]), (1.0/(1.0 - parms[76]))) * pow(x[39], (-parms[76]/(1.0 - parms[76]))) + std::max(x[39] - y[24], 0.0);
ydot[1] = x[3] * y[1];
x[7] = y[1] * x[5]/(y[7] * x[11]);
x[8] = y[1] * x[5]/(y[7] * x[11]);
x[12] = y[2]/(y[7] * x[11]);
x[14] = y[7] * x[11] - y[1] * x[5] - parms[6] * y[2] - y[7] * x[9] * y[4] - y[7] * x[21];
x[42] = x[41] * parms[84];
x[43] = x[41] * (1.0 - parms[85] - parms[84]);
x[44] = x[41] * parms[85];
x[13] = parms[22] * (parms[23] * (x[8] + parms[24]) - 1.0);
x[15] = std::max(0.0, x[14]/(y[7] * x[11]));
x[16] = x[14]/(y[7] * y[4]);
x[25] = x[39] + x[42] * parms[68] + x[43] * parms[69] + x[44] * parms[70];
x[48] = y[20] + x[42];
x[49] = y[21] + x[43];
x[50] = y[22] + x[44];
ydot[7] = x[13] * y[7];
ydot[17] = parms[65] * y[17] * (1.0 - y[17]/parms[63]) - x[48];
ydot[18] = parms[66] * y[18] * (1.0 - y[18]/parms[64]) - x[49];
ydot[19] = parms[67] * y[19] * (1.0 - y[19]/y[16]) - x[50] - x[46];
x[0] = std::max(parms[17], std::min(parms[18], parms[15] + parms[16] * x[15]));
x[2] = std::max(parms[27], std::min(parms[28], parms[25] + parms[26] * x[16]));
x[31] = x[25] + x[26] + x[30];
ydot[11] = x[31]/3.666 - y[11] * parms[46] + parms[46] * y[12] * parms[43]/parms[44];
x[1] = parms[19] * x[0] * x[11];
x[17] = x[2] * y[7] * y[4];
x[32] = x[31] - x[30];
ydot[4] = x[1] - x[9] * y[4] - x[54];
ydot[20] = parms[71] * x[51] * (((x[1] - x[9] * y[4] - x[54]) * y[0] - y[4] * parms[0] * y[0])/(y[4] * y[0]) + 1.0) - y[20];
ydot[21] = parms[72] * x[51] * (((x[1] - x[9] * y[4] - x[54]) * y[0] - y[4] * parms[0] * y[0])/(y[4] * y[0]) + 1.0) - y[21];
ydot[22] = parms[73] * x[51] * (((x[1] - x[9] * y[4] - x[54]) * y[0] - y[4] * parms[0] * y[0])/(y[4] * y[0]) + 1.0) - y[22];
x[18] = x[14] - x[17];
ydot[2] = y[7] * x[1] - x[18] - (y[7] * x[9] * y[4]);
}
	
Rcpp::NumericMatrix RK4(int nt, 
                      double byT,
                      std::vector<double> Ry0,
                      std::vector<double> Rparms, 
                      double** dataExogVar,
                      double** exogSamplingTime, 
                      int nExogVar) {
	int it, it1;
	double *y = &Ry0[0];
	double *parms = &Rparms[0];
	double y1[dim], y2[dim], y3[dim], ydot0[dim], ydot1[dim], ydot2[dim], ydot3[dim], ydots[dim], x0[dimIv], x1[dimIv], x2[dimIv], x3[dimIv];
	Rcpp::NumericMatrix out(nt, dimOut);

	for (it=0; it<dim;it++) { //init out vector
		out(0, it)=y[it];
	}
	int comptExogVar[nExogVar];
	for (it=0; it<nExogVar; it++) comptExogVar[it]=1;

	// get intermediateVar and compute distance at t=0 //
	Func(0, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar); 
						for (it1=0; it1<dim; it1++) {
							out(0, dim+it1) = ydot0[it1];
						}
						for (it1=0; it1<dimIv; it1++) {
							out(0, 2*dim+it1) = x0[it1];
						}
						 
	
	for (it=0; it<nExogVar; it++) comptExogVar[it]=1;
	
	for (it=0; it<(nt-1); it++) {

			
			

			Func(it*byT, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);

			for (it1=0; it1<dim; it1++)
				y1[it1] = y[it1] + ydot0[it1]*0.5*byT;
			Func((it + 0.5)*byT, y1, parms, ydot1, x1, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);

			for (it1=0; it1<dim; it1++)
				y2[it1] = y[it1] + ydot1[it1]*0.5*byT;
			Func((it + 0.5)*byT, y2, parms, ydot2, x2, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);

			for (it1=0; it1<dim; it1++)
				y3[it1] = y[it1] + ydot2[it1]*byT;
			Func((it+1)*byT, y3, parms, ydot3, x3, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);

			for (it1=0; it1<dim; it1++) {
				ydots[it1] = (ydot0[it1] + 2.0*ydot1[it1] + 2.0*ydot2[it1] + ydot3[it1])/6.0;
				y[it1] = y[it1] + byT*ydots[it1];
				out(it+1, it1) = y[it1];
			}
			
			for(it1=0;it1<dim;it1++){
							out(it+1, dim+it1) = ydots[it1];
						}
						for(it1=0;it1<dimIv;it1++){
							out(it+1, 2*dim+it1) = (x0[it1] + 2.0*x1[it1] + 2.0*x2[it1] + x3[it1])/6.0;
						}
				
	}
	return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RK4(int nt, 
                        double byT,
                        std::vector<double> Ry0,
                        std::vector<double> Rparms, 
                        Rcpp::List RdataExogVar,
                        Rcpp::List RexogSamplingTime) {
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);
	int nExogVar = RdataExogVar.size();
	Rcpp::NumericMatrix out = RK4(nt, byT, Ry0, Rparms, dataExogVar, exogSamplingTime, nExogVar);
	for (unsigned int it=0; it<RdataExogVar.size(); it++) {
		free(dataExogVar[it]);
		free(exogSamplingTime[it]);
	}
	free(dataExogVar);
	free(exogSamplingTime);
	
	return out;
}