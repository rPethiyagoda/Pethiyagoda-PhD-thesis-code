function [x,y,xB,yB]=getXY(yDash,N,deltaPhi,varargin)

tau = [0;yDash(1:(N-1))];
theta = [0;yDash(N:(2*N-2))];
phi = (-(N-1)/2:(N-1)/2)'*deltaPhi;
tauHalf = (tau(1:N-1)+tau(2:N))./2;
thetaHalf = (theta(1:N-1)+theta(2:N))./2;

y = zeros(N,1);
x0 = phi(1);
x = x0*ones(N,1);

for i=2:N
    y(i) = y(i-1)+deltaPhi*(exp(-tauHalf(i-1))*sin(thetaHalf(i-1)));
    x(i) = x(i-1)+deltaPhi*(exp(-tauHalf(i-1))*cos(thetaHalf(i-1)));
end

xB = [];
yB = [];

if length(varargin)==3
    b = varargin{1};
    alpha = varargin{2};
    beta = varargin{3};
    modi = beta*pi/alpha;
    weightP=intWeight(phi,'trap');
    int1 = @(phiPrime) intlenFuncX(phiPrime,b,modi,weightP,theta,phi);
    int2 = @(phiPrime) intlenFuncY(phiPrime,b,modi,weightP,theta,phi);
    
    humpPoints = 100;
    
    xB = zeros(humpPoints+2,1);
    yB = -alpha*ones(humpPoints+2,1);
    
    xB(1) = phi(1);
    
    len = integral(int1,phi(1),alpha/(beta*pi)*log(b));
    xB(2) = xB(1) + len;
    
    bphis = linspace(alpha/(beta*pi)*log(b),0,humpPoints);
    
    for i=1:humpPoints-1
        xB(i+2) = xB(i+1) + integral(int1,bphis(i),bphis(i+1));
        yB(i+2) = yB(i+1) + integral(int2,bphis(i),bphis(i+1));
    end
    
    xB(humpPoints+2) = max(x);
    yB(humpPoints+2) = yB(humpPoints+1);
end

xshift = x((N-1)/2+1);
xB = xB-xshift;
x = x-xshift;
end

function int = intlenFuncX(phiPrime,b,modi,weightP,theta,phi)
n = length(phiPrime);
int = zeros(size(phiPrime));
for i=1:n
    int(i) = exp(-tauFunc(phiPrime(i),b,modi,weightP,theta,phi)).*cos(thetaFunc(phiPrime(i),b,modi));
end
end

function int = intlenFuncY(phiPrime,b,modi,weightP,theta,phi)
n = length(phiPrime);
int = zeros(size(phiPrime));
for i=1:n
    int(i) = exp(-tauFunc(phiPrime(i),b,modi,weightP,theta,phi)).*sin(thetaFunc(phiPrime(i),b,modi));
end
end

function tau = tauFunc(phiPrime,b,modi,weightP,theta,phi)
tau = -1 - (2*-exp(modi*phiPrime)+b+1)/(2*(b-1))*log(abs((-exp(modi*phiPrime)+1)/(-exp(modi*phiPrime)+b)))...
    +modi/pi*sum(weightP.*theta./(1+exp(modi*(phiPrime-phi))));
end

function t = thetaFunc(phiPrime,b,modi)
t = pi/(b-1)*exp(modi*phiPrime)-pi*(b+1)/(2*(b-1));
t(phiPrime<log(b)/modi) = 0;
t(phiPrime>0) = 0;
end
