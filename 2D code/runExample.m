global mkl gpu;

startup_STB();
% Parameters------------------------------
Fr = 0.6; % The Froude number
bAll = 1.5:0.5:6; % Bootstrap up to b=6
alpha = pi; % nondimensional upstream depth
beta = -1;  % mapping type, 1: maps upstream to origin, -1: maps downstream to origin
intMethod = 'trap'; % Method of integration

mkl = false;
gpu = false;

%% Create the mesh-----------------------------
if gpu && mkl
    N = 48001;
    deltaPhi = 40/(N-1);

    band = 48000;
else
    N = 3001;
    deltaPhi = 20/(N-1);

    band = 3000;
end

yDash0 = zeros(2*(N-1),1);

yDash = computeSurface(N,deltaPhi,bAll,Fr,alpha,beta,intMethod,band,'dense',yDash0,'./Stats/');

%% Plot Bottom-------------------------------

[x,y,xB,yB] = getXY(yDash,N,deltaPhi,bAll(end),alpha,beta);
figure;plot(xB,yB,'k');