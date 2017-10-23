function yDash = computeSurface(N,deltaPhi,bAll,FrAll,alpha,beta,intMethod,band,preconType,yDash0,folder)
%%
% computeSurface
% N - Dimension of the mesh
% deltaPhi - The phi spacing of the mesh
% bAll - b from the hump definition
% FrAll - The Froude number
% alpha - nondimensional upstream depth
% beta - mappting type, -1: maps upstream to origin, 1: maps downstream to origin
% intMethod - the weighting schemes for numerical integration  - values 'trap', 'simpson', 'boole'
% band - the block-bandwidth of the preconditiontioner values [0,M-1]
% preconType - The type of preconditioner used - values 'band' (banded
% storage preconditioner, requires Intel MKL), 'dense' (dense storage) and
% 'full' (full Newton's method)
% yDash0 - The initial guess
% folder - The folder that the surface data is saved to

global PreSolveCount PreSolveTime LinPreRunTime NonLinFuncRunCount NonLinFuncRunTime PreSetupCount data dataFilepath nonLin kern baseMult gpu mkl;
setTo0(true);
nonLin = 0;
if band>(N+1)
    band=N;
end

bl = length(bAll);
Frl = length(FrAll);

if bl==1
    bAll = bAll*ones(1,Frl );
elseif Frl ==1
    FrAll = FrAll*ones(1,bl);
else
    error('Either hAll or FrAll must be of length 1');
end
b = bAll(1);
Fr = FrAll(1);



% Initialise Function---------------------
if gpu
gy = ceil((N-1)/60000);
gx = ceil((N-1)/gy);
kern = parallel.gpu.CUDAKernel('fGutsGPU.ptx','fGutsGPU.cu','guts');
kern.ThreadBlockSize = [512,1,1];               
kern.GridSize = [gx,gy,1];
kern.SharedMemorySize = 10000;
end

%%
data = [];
data.N = N;
data.J = [];
data.ipvt = [];
data.kl = [];
data.L = [];
data.U = [];
data.A = [];
data.ipvtA = [];
data.base = [];
data.B = [];
data.C = [];
data.Bcoef = [];
data.Ccoef = [];
data.D = [];
data.ipvtD = [];
data.E = [];
data.F = [];
data.G = [];
data.H = [];
data.I = [];
data.ps = [];
data.eEq = [];
data.Linv = [];



%number of equations
neq = 2*N-2;

% Maximum Krylov subspace size
maxl = 100;

% Maximum norm of the search direction
mxnewt = 50;


% Scaling vectors
fscale = [10^3*ones(N-1,1);ones(N-1,1)];
yscale = ones(neq,1);

strategy = 'LineSearch';

%%

% Set global values to zero
setTo0(true);

% Count the number of bootstrapping steps
[~,numEps] = size(bAll);

% Generate the preconditioner
preSetupTimet = tic;
if strcmp(preconType,'dense')
    ps = (N-1);
    data.D = zeros(ps,ps);
    psetfnNonlin = @ (y,yscale,Fy,fscale) NonlinearBlockJacPre(y,Fy,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band);
    psetfnNonlin(yDash0,0,0,0);
elseif strcmp(preconType,'band')
    ps = (N-1);
    kl = band;
    ku = kl;
    data.D = zeros(2*kl+ku+1,ps);
    psetfnNonlin = @ (y,yscale,Fy,fscale) NonlinearBandJacPre(y,Fy,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band);
    psetfnNonlin(yDash0,0,0,0);
end
preSetupTime = toc(preSetupTimet);


% itterate over all bootstrapping steps
yDash=yDash0;
for i=1:numEps
    
    yDash = zeros(2*(N-1),1);
    
    % Set global values (not includin preconditioner set up time) to zero
    setTo0(true);
    JFNKLinBandTime = tic;
    b = bAll(:,i);
    Fr = FrAll(:,i);

    disp([num2str(N),' Delta x=',...
        num2str(deltaPhi(1)),...
        ', b=',num2str(b),', F=',num2str(Fr),', PreType=',preconType,' imxy=',intMethod]);

    
    f = @(r) F2D(r,Fr,N,deltaPhi,b,alpha,beta,intMethod);
    
    % Set up KINSol options depeding on the preconditioner storage type
    if strcmp(preconType,'dense')   % Dense storage
        psetfnNonlin = @ (y,yscale,Fy,fscale) NonlinearBlockJacPre(y,Fy,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band);
        optionsFull = KINSetOptions('MaxNewtonStep',mxnewt, ...
                                'LinearSolver', 'GMRES', ...
                                'KrylovMaxDim', maxl,...
                                'MaxNumSetups',20,...
                                'PrecSetupFn',psetfnNonlin,...
                                'InitialSetup',false,...
                                'PrecSolveFn',@preSolveBlock);
    elseif strcmp(preconType,'band')    %Banded storage
        psetfnNonlin = @ (y,yscale,Fy,fscale) NonlinearBandJacPre(y,Fy,N,deltaPhi,Fr,alpha,beta,intMethod,yscale,fscale,band);
        optionsFull = KINSetOptions('MaxNewtonStep',mxnewt, ...
                                'LinearSolver', 'GMRES', ...
                                'KrylovMaxDim', maxl,...
                                'PrecSetupFn',psetfnNonlin,...
                                'InitialSetup',false,...
                                'PrecSolveFn',@preSolveBnd);
    else                                % Finite difference Jacobian
        dnsjac = @ (y,Fy) finite_difference_jacobian(f,y,Fy);
        optionsFull = KINSetOptions('MaxNewtonStep',mxnewt, ...
                        'JacobianFn',dnsjac,...
                        'LinearSolver', 'Dense');
        
        
    end

    KINInit(f, neq, optionsFull);
    [status, yDash] = KINSol(yDash, strategy, yscale, fscale);
        

    
    JFNKStats = KINGetStats();
    totalJFNKLinBandTime = toc(JFNKLinBandTime);
    KINFree;

    clear optionsJNFKBandLinAlg;


    nz = ((2*band+1)*(N-1)-band*(band+1));
    percent = nz/(N-1)^2*100;
    if i==1
        totalTime = totalJFNKLinBandTime+preSetupTime;
        LinPreRunTime = LinPreRunTime+preSetupTime;
    else
        totalTime = totalJFNKLinBandTime;
    end
    nonLinStep = JFNKStats.nni;
    try
    avgSspace = JFNKStats.LSInfo.nli/nonLinStep;
    disp([num2str(band),' & ',num2str(percent),' & ',num2str(JFNKStats.LSInfo.npe),' & ',num2str(NonLinFuncRunCount),...
        ' & ',num2str(nonLinStep),' & ',num2str(avgSspace),...
        ' & ',num2str(LinPreRunTime),' & ',num2str(PreSolveTime),...
        ' & ',num2str(NonLinFuncRunTime),' & ',num2str(totalTime),'\\']);
    catch
        disp([num2str(band),' & ',num2str(percent),' & ',num2str(NonLinFuncRunCount),...
        ' & ',num2str(nonLinStep),' & ','N/A',...
        ' & ',num2str(LinPreRunTime),' & ',num2str(PreSolveTime),...
        ' & ',num2str(NonLinFuncRunTime),' & ',num2str(totalTime),'\\']);
    end
    
    timestamp = datestr(now,'yyyy-mm-dd (HH.MM.SS)');

    filepath = [folder,'N=',num2str(N),...
        ' b ',num2str(b),' F ',num2str(Fr),...
        ' dphi ',num2str(deltaPhi(1)),...
        ' band ',num2str(band),' im-',intMethod,' - ',timestamp,'.mat'];
    
    save(filepath,'N','deltaPhi','b','Fr','alpha','beta','intMethod','band','preconType','yDash0',...
        'JFNKStats','yDash','LinPreRunTime','NonLinFuncRunCount','NonLinFuncRunTime','PreSolveCount','PreSolveTime',...
        'fscale','maxl','mxnewt','neq','nonLinStep','nz','percent','status','strategy',...
        'timestamp','totalJFNKLinBandTime','totalTime','yscale');
end
end