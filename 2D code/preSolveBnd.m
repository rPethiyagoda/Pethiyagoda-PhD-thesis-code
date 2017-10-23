function [z,Flag]=preSolveBnd(y,yscale,Fy,fscale,r)

global PreSolveCount PreSolveTime data nonLin dataFilepath bMult cMult;

PreSolveCount = PreSolveCount + 1;
PreSolveTimeTemp = tic;

ps = data.ps;
Flag =0;
r1 = r(1:ps);
r2 = r((ps+1):(2*ps));


% even better
y1 = data.A\r1;
x2 = r2-cMult(y1);
y2 = solveBand(data.D,data.ipvtD,x2,data.kl,data.kl,0);
tempz1 = bMult(y2);
z1 = y1-data.A\tempz1;
z = [z1;y2];

PreSolveTime = PreSolveTime + toc(PreSolveTimeTemp);

end