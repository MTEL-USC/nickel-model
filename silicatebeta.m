function [A,b,output] = silicatebeta(A,b,do)

%SILICATEBETA sets uptake and remineralization of the element in silicate,
%similar to biobeta

fprintf('%s','silicatebeta...')

% unpack the key parameters
silicatebeta = do.silicatebeta.silicatebeta;

% load the grid and the model output for P uptake (productivity)

load ([do.highestpath '/data/ao.mat'])
load ([do.highestpath '/seth_Ni/silica/SI_UP_HLZ.mat'])
load ([do.highestpath '/seth_Ni/silica/SI_REM_HLZ.mat'])

% convert from molar to micromolar
SI_UP_HLZ = SI_UP_HLZ*1e6;
SI_REM_HLZ = SI_REM_HLZ*1e6;

% establish the relative uptake rate constant for element E (where E uptake
% is given by uptake = KPROD*[E])
KPROD = SI_UP_HLZ*silicatebeta;
kprod = KPROD(ao.iocn);

% create the A matrix for loss of E by biological uptake
prodA = sparse(1:ao.nocn,1:ao.nocn,kprod,ao.nocn,ao.nocn);

% calculate the particle concentration profile as the fraction remaining
% compared to zc; this is the particle concentration at the interfaces
% between boxes so this matrix is also 91x180x25
PARTICLES_zc = nansum((SI_UP_HLZ.*ao.VOL),3); % convert from uM y-1 to umoles y-1
PARTICLES_zc(PARTICLES_zc==0) = 0.001; % make sure that you have at least a teeny-tiny amount of particles everywhere in the ocean, so that you don't get nans when you divide by PARTICLES_zc
NET_REM = cumsum((SI_REM_HLZ.*ao.VOL),3); % calculate the total amount of remineralization above the bottom boundary of each box 
NET_REM = cat(3,zeros(91,180,1),NET_REM); % calculate the remineralization at all box boundaries
PARTICLES_zc = repmat(PARTICLES_zc,[1,1,25]);
PARTICLES = PARTICLES_zc-NET_REM;
PARTICLES = PARTICLES./PARTICLES_zc; %convert to the fractional amount of particles compared to the zc depth

% calculate the particle flux divergence (PFD), which is the difference
% between the amount of particles present at the top of the grid cell
% compared to at the bottom of the grid cell; this is of course equal to
% the particle dissolution which occurred within that grid cell; note that
% we are continuing to express everything about particles in terms of their
% fractional abundance compared to the abundance at zc
PFD = -diff(PARTICLES,1,3);

% no particles dissolve in the top two layers, so set the PFD here to zero
PFD(:,:,1:2)=0;

% set the PFD to zero wherever there is not ocean
PFD = PFD.*ao.OCN;

% because productivity occurs in the top two layers, we are going to
% calculate two different A matrices, one which describes the
% remineralization of element E based on biological productivity which
% occurred in the top layer, and then again similarly for the second layer

% create a shoebox containing the productivity rate constants to which
% remineralization is related; for example when considering the
% remineralization which occurs in grid-column (1,1,:), the
% remineralization at every depth within that grid-column will be related
% to the productivity which occurred in the top grid cell (1,1,1), and
% seprarately to the productivity which occurred in grid cell (1,1,2)
KPROD1 = repmat(KPROD(:,:,1),1,1,24); KPROD2 = repmat(KPROD(:,:,2),1,1,24);

% similarly, the equation position (column in the A matrix) for every grid
% cell in the water column (1,1,:) will depend on the equation position of
% the top grid cell (1,1,1) from which particles originate
FROM1 = repmat(ao.EQNPOS(:,:,1),1,1,24); FROM2 = repmat(ao.EQNPOS(:,:,2),1,1,24);
from1 = FROM1(ao.iocn); from2 = FROM2(ao.iocn);

% the equation position (row in the A matrix) in which we slot the
% remineralization depends on the equation position of the cell to which
% particles go
TO1 = ao.EQNPOS; TO2 = ao.EQNPOS;
to1 = TO1(ao.iocn); to2 = TO2(ao.iocn);

% the magnitude of the remineralization is the surface productivity,
% multiplied by the PFD (fraction of that productivity with remineralizes),
% corrected for the height different between the two grid cells
remin1 = KPROD1(ao.iocn).*PFD(ao.iocn).*ao.height(1)./ao.HEIGHT(ao.iocn); remin2 = KPROD2(ao.iocn).*PFD(ao.iocn).*ao.height(2)./ao.HEIGHT(ao.iocn);

% finally, we build the sparse A matrices for remineralization 
remin1A = sparse(to1,from1,remin1,ao.nocn,ao.nocn); remin2A = sparse(to2,from2,remin2,ao.nocn,ao.nocn);

% modify the A matrix for productivity and water-column remineralization
A = A - prodA + remin1A + remin2A;

% package outputs
output.silicatebeta=silicatebeta;
output.prodA=prodA;
output.remin1A=remin1A;
output.remin2A=remin2A;
output.citations={'biological cycling parameters derived from Thomas Weber, Seth John, Alessandro Tagliabue, and Tim DeVries, Biological uptake and reversible scavenging of zinc in the global ocean, Science, 2018, DOI: 10.1126/science.aap8532.'};