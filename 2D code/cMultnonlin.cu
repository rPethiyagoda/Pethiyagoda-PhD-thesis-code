#define PI 3.141592653589793238462643
#define blocDim 256
#define powOfTwo 4
#define timerCount 10
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

__global__ void guts(const long int N,const double * x,double * b)
{
	long int threadPos;

	//__syncthreads();

	threadPos = blockIdx.y*gridDim.x+blockIdx.x;
	//threadPos = threadIdx.x;
	if(threadPos<(N-1)){
		if(threadPos>0){
			b[threadPos] = 1.0/2.0*(x[threadPos-1]+x[threadPos]);
		}else{
			b[threadPos] = 1.0/2.0*x[threadPos];
		}
	//	threadPos += blockDim.x;
	}
	
	
}
