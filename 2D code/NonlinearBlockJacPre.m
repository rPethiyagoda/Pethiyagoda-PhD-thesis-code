function Flag = NonlinearBlockJacPre1(variables,Fx,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band)
global LinPreRunTime data bMult cMult gpu mkl;
fprintf('making preconditioner');

LinAlgPreRunTimeTemp = tic;

Flag = 0;

%% Initialise variables/ generate the mesh
% Block matric setup
% | A B |
% | C D |
ps = (N-1);

tau = [0;variables(1:(N-1))];
theta = [0;variables(N:(2*(N-1)))];


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

for k = 1:(N-1)
    temp = weiPhi./(1-exp(beta*pi/alpha*(phiHalf(k)-phi')));
    tempCol = weiPhi(k+1)./(1-exp(beta*pi/alpha*(phiHalf-phi(k+1))));
    data.D(1:(N-1),k) = -beta/alpha*tempCol;
    extra = 1/2*(beta/alpha*sum(temp)...
        -1/pi*log(abs((exp(beta*pi/alpha*phi(end))-exp(beta*pi/alpha*phiHalf(k)))...
        /(exp(beta*pi/alpha*phi(1))-exp(beta*pi/alpha*phiHalf(k))))));
    
    data.D(k,k) = data.D(k,k) + extra;
    if k>1
        data.D(k,k-1) = data.D(k,k-1) + extra;
    end
        
    
    bandStart = max((k-band),1);
    bandEnd = min((k+band),ps);
    data.D(1:(bandStart-1),k)=0;
    data.D((bandEnd+1):end,k)=0;
    
end

fprintf('.');

%% Form the Schur compliment
for i=1:ps
    temp = zeros(ps,1);
    temp(i) = 1;
    tempPart = A\bMult(temp);
    Dtemp = cMult(tempPart);
    
    data.D(:,i) = data.D(:,i)-Dtemp;
    
    
    
end
fprintf('.');

%% Factor matrices
ipvtD = zeros(ps,1,'int64');

if mkl
    factorDense(data.D,ipvtD);
else
    [data.L,data.U] = lu(data.D);
end

fprintf('.');

%% Assign the proconditioner to a global variable
data.A = A;
data.ipvtD = ipvtD;
data.ps = ps;

LinPreRunTime = LinPreRunTime + toc(LinAlgPreRunTimeTemp);
fprintf('.');
disp('finished');
end


function K = K3(x,xHalf)

xDiff = x-xHalf;

K=log(abs(xDiff));


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


