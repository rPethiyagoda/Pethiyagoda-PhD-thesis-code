function plotSurf(yDash,N,deltaPhi,varargin)
% plotSurf Plots the surface that is the solution to the problem.
% yDash - vector of zetaX
% N - the number of nodes in the y and x directions respectively
% deltaPhi - The phi spacing of the mesh
% [b] - b from the hump definition
% [alpha] - nondimensional upstream depth
% [beta] - mappting type, -1: maps upstream to origin, 1: maps downstream to origin

% Making sure the function is the correct length (it may be less given
% computational simplifications).

[x,y,xB,yB] = getXY(yDash,N,deltaPhi,varargin{:});

figure;
plot(x,y);
hold on;
plot(xB,yB,'k');