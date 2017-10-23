function [z,Flag]=preSolveBlock(y,yscale,Fy,fscale,r)
global PreSolveCount PreSolveTime data bMult cMult mkl;

PreSolveCount = PreSolveCount + 1;
PreSolveTimeTemp = tic;
ps = data.ps;
Flag =0;
r1 = r(1:ps);
r2 = r((ps+1):(2*ps));


y1 = data.A\r1;
x2 = r2-cMult(y1);
if mkl
    y2 = solveDense(data.D,data.ipvtD,x2);
else
    y2 = data.U\(data.L\x2);
end
z1 = y1-data.A\(bMult(y2));
z = [z1;y2];


PreSolveTime = PreSolveTime + toc(PreSolveTimeTemp);