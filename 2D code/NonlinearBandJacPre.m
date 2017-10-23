function [Flag] = NonlinearBandJacPre1(variables,Fx,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band)

global LinPreRunTime PreSetupCount data bMult cMult gpu;

fprintf('making preconditioner');

LinAlgPreRunTimeTemp = tic;

Flag = 0;

%% Initialise variables/ generate the mesh
% | A B |
% | C D |
ps = (N-1);

tau = [0;variables(1:(N-1))];
theta = [0;variables(N:(2*(N-1)))];
% b = variables(2*N-1);
% a = sqrt(b);
% b = a^2;

tauHalfDiff = (tau(2:N)-tau(1:(N-1)))/deltaPhi;
tauHalf = (tau(1:N-1)+tau(2:N))./2;
thetaHalf = (theta(1:N-1)+theta(2:N))./2;

phi = (-(N-1)/2:(N-1)/2)'*deltaPhi;


phiHalf = (phi(1:N-1)+phi(2:N))./2;

weiPhi=intWeight(phi,intMethod)';

fprintf('.');
%% Sub matrix A
part1 = 3/2*exp(3*tauHalf).*tauHalfDiff;
part2 = 1/(deltaPhi)*exp(3*tauHalf);
temp=Fr^2*[part1-part2,part1+part2];
bandTemp = temp;
bandTemp((band+1):(band-1):end,1)=0;
A = spdiags(bandTemp,[1,0],N-1,N-1)';

fprintf('.');
%% Sub matrix B/C base
if gpu
    gy = ceil((N-1)/60000);
    gx = ceil((N-1)/gy);

    kernbMult = parallel.gpu.CUDAKernel('bMultnonlin.ptx','bMultnonlin.cu','guts');
    kernbMult.ThreadBlockSize = [1,1,1];               
    kernbMult.GridSize = [gx,gy,1];
    kernbMult.SharedMemorySize = 10000;

    bMult = @(r) gather(feval(kernbMult,N*ones(1,1,'int64'),thetaHalf,r,zeros(ps,1)));

    kerncMult = parallel.gpu.CUDAKernel('cMultnonlin.ptx','cMultnonlin.cu','guts');
    kerncMult.ThreadBlockSize = [1,1,1];               
    kerncMult.GridSize = [gx,gy,1];
    kerncMult.SharedMemorySize = 10000;

    cMult = @(r) gather(feval(kerncMult,N*ones(1,1,'int64'),r,zeros(ps,1)));

else
    B = spdiags(repmat(cos(thetaHalf)/2,1,2),[1,0],N-1,N-1)';
    bMult = @(r) B*r;
    C = spdiags(1/2*ones(N-1,2),[1,0],N-1,N-1)';
    cMult = @(r) C*r;
end

fprintf('.');
%% Sub matrix D
kl = band;
ku = kl;
centre = kl+ku+1;

for k = 1:(N-1)
    data.D(:,k) = 0;
    temp = weiPhi./(1-exp(beta*pi/alpha*(phiHalf(k)-phi')));
    tempCol = weiPhi(k+1)./(1-exp(beta*pi/alpha*(phiHalf-phi(k+1))));
    
    bandStart = max((k-band),1);
    bandEnd = min((k+band),(N-1));
    sb = centre + bandStart - k;
    eb = centre + bandEnd - k;
    
    data.D(sb:eb,k) = -beta/alpha*tempCol(bandStart:bandEnd);
    
    
    extra = 1/2*(beta/alpha*sum(temp)...
        -1/pi*log(abs((exp(beta*pi/alpha*phi(end))-exp(beta*pi/alpha*phiHalf(k)))...
        /(exp(beta*pi/alpha*phi(1))-exp(beta*pi/alpha*phiHalf(k))))));
    
    data.D(centre,k) = data.D(centre,k) + extra;
    if k>1
        data.D(centre+1,k-1) = data.D(centre+1,k-1) + extra;
    end
    
end


fprintf('.');

%% Form the Schur compliment

for k=1:(N-1)
    temp = zeros(N-1,1);
    temp(k) = 1;
    tempPart = A\bMult(temp);
    Dtemp = cMult(tempPart);

    
    bandEnd = min((k+band),(N-1));
    eb = centre + bandEnd - k;
    
    data.D(centre:eb,k) = data.D(centre:eb,k)-Dtemp(k:bandEnd);
end

fprintf('.');

%% Factor matrices
ipvtD = zeros(ps,1,'int64');
factorBand(data.D,ipvtD,kl,ku);
fprintf('.');

%% Assign the proconditioner to a global variable
data.A = A;
data.ipvtD = ipvtD;
data.kl = kl;
data.ps = ps;


PreSetupCount = PreSetupCount+1;
LinPreRunTime = LinPreRunTime + toc(LinAlgPreRunTimeTemp);
fprintf('.');
disp('finished');
end

function B = convert2band(A)
% Lower triangular
[m,n] = size(A);

B = zeros(m,n);
for i=1:n
    B(1:(end-i+1),i)=A(i:end,i);
end

end

function K = K3(x,xHalf)

xDiff = x-xHalf;

K=log(abs(xDiff));


end


function P = multiBlock(A,k);
P = A;
for i = 2:k
    P = blkdiag(P,A);
end

end

% Determines one half of the closed form integral
function val=closedIntt(s,t,A,B,C)

% if abs(t)<10^-15
%     val = 0;
% else
    val= t/sqrt(A).*log(2*A*s+B*t+2*sqrt(A*(A*s.^2+B*s.*t+C*t.^2)));
% end
    val(abs(t)<10^-15) = 0;
end

% Determines another half of the closed form integral
function val=closedInts(s,t,A,B,C)

% if abs(s)<10^-15
%     val = 0;
% else
    val= s/sqrt(C).*log(2*C*t+B*s+2*sqrt(C*(A*s.^2+B*s.*t+C*t.^2)));
% end
    val(abs(s)<10^-15) = 0;
end

% Determines the closed form integral
function singInt=definitInt(sLower,sUpper,tLower,tUpper,A,B,C)

singInt =   (closedIntt(sUpper,tUpper,A,B,C)+closedInts(sUpper,tUpper,A,B,C))-...
            (closedIntt(sLower,tUpper,A,B,C)+closedInts(sLower,tUpper,A,B,C))-...
            (closedIntt(sUpper,tLower,A,B,C)+closedInts(sUpper,tLower,A,B,C))+...
            (closedIntt(sLower,tLower,A,B,C)+closedInts(sLower,tLower,A,B,C));

end


