function [z,Flag]=preSolveBlock(y,yscale,Fy,fscale,r)
% Apply the preconditioner matrix to a vector r

global PreSolveCount PreSolveTime data baseMult mkl;

PreSolveCount = PreSolveCount + 1;
PreSolveTimeTemp = tic;

Flag =0;
l = length(r);
r1 = r(1:l/2);
r2 = r((l/2+1):l);


if mkl
    y1 = solveBand(data.A,data.ipvtA,r1,1,1,1);
    x2 = r2-data.Ccoef*baseMult(y1);
    y2 = solveDense(data.D,data.ipvtD,x2);
    tempz1 = data.Bcoef*baseMult(y2);
    z1 = y1-solveBand(data.A,data.ipvtA,tempz1,1,1,1);
    z = [z1;y2];
else
    y1 = data.A\r1;
    x2 = r2-data.Ccoef*baseMult(y1);
    y2 = data.D\x2;
    tempz1 = data.Bcoef*baseMult(y2);
    z1 = y1-data.A\tempz1;
    z = [z1;y2];
end
PreSolveTime = PreSolveTime + toc(PreSolveTimeTemp);