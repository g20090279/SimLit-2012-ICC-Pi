function [H, info] = generate_channel (numAntTx, numAntRx)
% Last revised: Mar 12, 2023
% Zekai Liang, liang@ice.rwth-aachen.de

%% channel model information
numRay = 20;  % consists of 20 rays
% angle of departure: uniformly distributed (-pi/3,pi/3), modeling
% 120-degree sector
angleDep = -pi/3 + (pi/3 - (-pi/3)) * rand(numRay,1);
% angle of arrival: uniformly distributed (0,2pi)
angleArr = 2*pi.*rand(numRay,1);
power = 1/numRay;  % power normalization factor, the power of each path

% complex Guassian channel coefficiet
coeff = 1/sqrt(2)*randn(numRay,1) + 1i*1/sqrt(2)*randn(numRay,1);

%% antenna array information
d = 1/2;   % ULA with half wavelength spacing

% array factor, 
arrFactorTx = exp(1i*(0:numAntTx-1)'*2*pi*d*sin(angleDep.'));  % array factor Tx
arrFactorRx = exp(1i*(0:numAntRx-1)'*2*pi*d*sin(angleArr.'));  % array factor Rx


%% generate channel matrix, frequency flat
H = sqrt(power) * arrFactorRx * diag(coeff) * arrFactorTx';
info = struct(                  ...
    'numRay', numRay,             ...
    'angleDep', angleDep,         ...
    'angleArr', angleArr,         ...
    'power', power,               ...
    'coeff', coeff,               ...
    'spaceAnt', d,                ...
    'arrFactorTx', arrFactorTx,   ...
    'arrFactorRx', arrFactorRx    ...
    );
end