#define PI 3.141592653589793238462643
#define blocDim 512
#define powOfTwo 4
#define timerCount 10


// Main GPU kernel
__global__ void guts(const long int N,
	const double * theta,const double * thetaHalf,
	const double * tauHalf,const double * tauHalfDiff,
	const double * phi,const double * phiHalf,
	const double * weiPhi,const double b,
	const double alpha,const double beta,
	const double Fr,double * Func)
{
	__shared__ double pH,tH,modi,thisBlock[blocDim];
	__shared__ long int k;
	double I1,FuncTemp = 0;
	long int i,threadPos;

	// First thread initialises variables
	if(threadIdx.x==0){
		k = blockIdx.y*gridDim.x+blockIdx.x;
		modi = beta*PI/alpha;
		if(k<(N-1)){
			pH = phiHalf[k];
			tH = thetaHalf[k];
		}
	}

	// After initialising variables
	__syncthreads();

	if(k<(N-1)){
		// Have each thread sum over some values of the double integrals
		threadPos = threadIdx.x;
		i = threadPos;
		// Each loop is one collocation point
		while(threadPos<N){

			// Calculate contributions to the integral
			I1 = weiPhi[i]*(theta[i]-tH)/(1-exp(modi*(pH-phi[i])));//

			// Accumulate integral contributions
			FuncTemp += I1;
		
			// Update collocation point
			threadPos += blockDim.x;
			i = threadPos;
		}


		// Add total contribution from thread to a storage vector
		thisBlock[threadIdx.x] = FuncTemp;

		// All threads finished evaluating the BIE
		__syncthreads();
   
		// Sum up all thread contributions
		for(i=blocDim/2;i>0;i=i/2){
			if(threadIdx.x<i){
				thisBlock[threadIdx.x] += thisBlock[threadIdx.x+i];
			}
			__syncthreads();
		}

		// Store complete BIE in correct loction of output vector
		if(threadIdx.x==(blockDim.x-1)){
			Func[k+(N-1)] = tauHalf[k]
				-(1 +(2*exp(modi*pH)+b+1)/(2*(b-1))*log((exp(modi*pH)+1)/(exp(modi*pH)+b)))
				-beta/alpha*thisBlock[0]
				-tH/PI*log(abs((exp(modi*phi[N-1])-exp(modi*pH))/(exp(modi*phi[0])-exp(modi*pH))));
		}


		// Fist block calculates Bernoulli's equation
		if(k==0){

			// Have each thread compute Bernoulli's equation for a mesh half point
			i = threadIdx.x;
			while(i<(N-1)){

				// Bernoulli's Equation
				Func[i]=Fr*Fr*exp(3*tauHalf[i])*tauHalfDiff[i]+sin(thetaHalf[i]);
		
				i += blockDim.x;
			}

		}
	}
	
}
