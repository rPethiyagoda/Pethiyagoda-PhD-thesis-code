function [Func,Flag] = F2D(variables,Fr,N,deltaPhi,b,alpha,beta,intMethod)
% F2D The function that needs to be minimised. Takes
% the following inputs
% variables - vector of tau and theta
% Fr - the Froude number
% N - the number of nodes in the phi direction
% deltaPhi - the distance between nodes in the phi direction
% b - b from the hump definition
% alpha - nondimensional upstream depth
% beta - mappting type, -1: maps upstream to origin, 1: maps downstream to origin
% intMethod - a string storing the method of integration

% Define global variables for metrics and the GPU kernel
global NonLinFuncRunCount NonLinFuncRunTime kern gpu;

NonLinFuncRunCount = NonLinFuncRunCount + 1;
NonLinFuncRunTimeTemp = tic;

% Set success
Flag = 0;
tau = [0;variables(1:(N-1))];
theta = [0;variables(N:(2*(N-1)))];
phi = (-(N-1)/2:(N-1)/2)'*deltaPhi;

% Calcualte half mesh points
phiHalf = (phi(1:N-1)+phi(2:N))./2;
tauHalf = (tau(1:N-1)+tau(2:N))./2;
tauHalfDiff = (tau(2:N)-tau(1:(N-1)))/deltaPhi;
thetaHalf = (theta(1:N-1)+theta(2:N))./2;

% Use auxiliary function to compute correct integration weightings
weightP=intWeight(phi,intMethod);

if gpu
    Func = zeros(2*(N-1),1);
    % Compute function on the GPU
    [dFunc] = feval(kern,N*ones(1,1,'int64'),...
        theta,thetaHalf,tauHalf,tauHalfDiff,...
        phi,phiHalf,...
        weightP,b,alpha,beta,Fr,Func);
    % Copy results back to host memory
    Func = gather(dFunc);

else
    modi = beta*pi/alpha;
    %Bernoulli's equation
    Func1 = Fr^2*exp(3*tauHalf).*tauHalfDiff+sin(thetaHalf);
    
    % Boundary integral equation
    Func2 = zeros(N-1,1);
    parfor i=1:(N-1)
        integrand = (theta-thetaHalf(i))./(1-exp(modi*(phiHalf(i)-phi)));%
        integral = sum(weightP.*integrand);
        Func2(i) = tauHalf(i) - 1 -(2*exp(modi*phiHalf(i))+b+1)/(2*(b-1))*log((exp(modi*phiHalf(i))+1)/(exp(modi*phiHalf(i))+b))...
            -beta/alpha*integral-1/pi*thetaHalf(i)*log(abs((exp(modi*phi(end))-exp(modi*phiHalf(i)))/(exp(modi*phi(1))-exp(modi*phiHalf(i)))));
    end

    Func = [Func1;Func2];

end


NonLinFuncRunTime = NonLinFuncRunTime + toc(NonLinFuncRunTimeTemp);
end