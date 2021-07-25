function [A,b,output] = revscavPa(A,b,do);
%REVSCAVPA scavenges your element downwards by reversible scavenging
%   This function adds reversible scavenging of your element. It is based
%   on the assumption that your element scavenges like Th from the Hulten model.

fprintf('%s','revscavPa...')

% load the grid and POC data (units of mmole m-3)
load ([do.highestpath '/data/ao.mat'])
load ESRPa

% unpack the scavenging rate constant
R = do.revscavPa.R;
sedremin = do.revscavPa.sedremin;

% get the POC concentrations from each grid cell
ESRPa = ESRPa(ao.iocn);

% the amount of element E which sinks out depends on the equilibrium
% scavenging constant, multiplied by the POC concentration, times the
% sinking rate divided by the height of the grid cell
sinkout = R*(ESRPa./ao.Height);

% find the equation position (positions in the A matrix) of the grid cells
% which lie below each cell
EQNPOSBELOW = cat(3,ao.EQNPOS(:,:,2:24),zeros(91,180,1));

% define the equation positions that particles are sinking from, as well as
% the volumes and heights of those grid cells
frompos = ao.EQNPOS(EQNPOSBELOW~=0);
fromvol = ao.Vol(frompos);
fromheight = ao.Height(frompos);

% define the equation positions that particles are sinking to, as well as
% the volumes and heights of those grid cells
topos = EQNPOSBELOW(EQNPOSBELOW~=0);
tovol = ao.Vol(topos);

% create the sinkout A matrix, and fill in the diagonal with the magnitude
% of the sinking flux out
sinkoutA = speye(ao.nocn,ao.nocn);
sinkoutA(sinkoutA==1)=sinkout;

% calculate the amount of element transferred into each grid cell by
% sinking with K, the sinking rate devided by the grid cell height from
% which sinking occurs, and a correction for volume
sinkin = sinkout(frompos).*(fromvol./tovol);

% create the A matrix for sinking in
sinkinA = sparse(topos,frompos,sinkin,ao.nocn,ao.nocn);

% find the equation positions of cells which lie on the bottom, and the
% amount of reminerlization there is equal to the amount which sinks out of
% that grid cell
btmeqnpos = ao.EQNPOS(ao.ibtm);
sedreminA = sparse(btmeqnpos,btmeqnpos,sinkout(btmeqnpos),ao.nocn,ao.nocn);

% add the A matrix with the sinking matrices
A = A - sinkoutA + sinkinA;

% add sedimentary remineralization if switched on
if sedremin
    A = A + sedreminA;
end

% package outputs
output.R=R;
output.sinkoutA=sinkoutA;
output.sinkinA=sinkinA;
if sedremin
    output.sedreminA=sedreminA;
end
output.citations=cell(1,1);