% Main file
% Last revised: Mar 13, 2023
% Zekai Liang, liang@ice.rwth-aachen.de
close all;
%% figure 1
numAntTx = 16;
numAntRx = 8;

[chn, chnInfo] = generate_channel(16,8);  % generate channel
F = chn' * chn;  % channel Hermitian times itself
[~, chnPowInd] = sort(abs(chnInfo.coeff),'descend'); % ordering power of rays

% snr in dB, since signal power is 1, the snr is the inverse of noise power
snrVec = (-10:5:30)';
snrVecLen = length(snrVec);

% SVD
[U,Z,V] = svd(chn);

numStream = [1,2,4];

capSVD = zeros(length(snrVec),length(numStream));
capSMF = zeros(length(snrVec),length(numStream));
capProp= zeros(length(snrVec),length(numStream));

for iStr = 1 : length(numStream) % iteration of different number of streams
    % SVD
    powConsPerStr = 1/numStream(iStr);  % total-power constraint per stream
    bfSVD = V(:,1:numStream(iStr));
    bfSVDnorm = bfSVD ./ vecnorm(bfSVD) * sqrt(powConsPerStr);  % normalize

    % Spatial matched filtering
    powConsPerStr = 1/numStream(iStr);  % total-power constraint
    bfSMF = chnInfo.arrFactorTx(:,chnPowInd(1:numStream(iStr)));
    bfSMFnorm = bfSMF ./ vecnorm(bfSMF) * sqrt(powConsPerStr);

    % Proposed
    powConsPerAnt = 1/(numAntTx*numStream(iStr));   % per-antenna power constraint
    bfProp = randn(numAntTx,numStream(iStr)) + 1i*rand(numAntTx,numStream(iStr));
    bfPropNorm = bfProp ./ abs(bfProp) * sqrt(powConsPerAnt);  % same power constraint for every antenna element
    
    % iteration of different SNR values
    for iSnr = 1:length(snrVec)
        % SVD
        capSVD(iSnr,iStr) = log(abs(det( eye(numStream(iStr)) + 10^(snrVec(iSnr)/10) * (bfSVDnorm'*F*bfSVDnorm) )));
        
        % Spatial matched filtering
        capSMF(iSnr,iStr) = log(abs(det( eye(numStream(iStr)) + 10^(snrVec(iSnr)/10) * (bfSMFnorm'*F*bfSMFnorm) )));
        
        % Proposed
        capProp(iSnr,iStr) = log(abs(det( eye(numStream(iStr)) + 10^(snrVec(iSnr)/10) * (bfPropNorm'*F*bfPropNorm) )));
        e = inf;  % difference of capacity
        eps = 10^-3; % preset threshold for stopping iteration
        while e > eps
            if numStream(iStr) == 1  % single-stream, Table 1
                    for iAntTx = 1:numAntTx
                        f = F(iAntTx,:);
                        f(iAntTx) = 0;  % k != i
                        bfTmp = bfPropNorm;
                        bfTmp(iAntTx) = 0;
                        bfPropNorm(iAntTx) = exp(1i*angle(f*bfTmp)) * sqrt(powConsPerAnt);
                    end
                    % new capacity after iteration update
                    capPropNew = log( 1 + 10^(snrVec(iSnr)/10) * abs(bfPropNorm'*F*bfPropNorm) );
            else % multi-stream
                for kStr = 1:numStream(iStr)
                    Vk = bfPropNorm; % submatrix of beamforming matrix V without k-th column
                    Vk(:,kStr) = [];
                    Ck = eye(numStream(iStr)-1) + Vk'*F*Vk;
                    Gk = F - 10^(snrVec(iSnr)/10)*F*Vk*Ck^(-1)*Vk'*F; % for k-th stream
                    for iAntTx = 1:numAntTx
                        gk = Gk(iAntTx,:);
                        gk(iAntTx) = 0;
                        bfPropTmp = bfPropNorm(:,kStr);
                        bfPropTmp(iAntTx) = 0;
                        bfPropNorm(iAntTx,kStr) = exp(1i*angle(gk*bfPropTmp))*sqrt(powConsPerAnt);
                    end
                end
                capPropNew = log(abs(det( eye(numStream(iStr)) + 10^(snrVec(iSnr)/10)*(bfPropNorm'*F*bfPropNorm) )));
            end
            % difference between last iteration and current iteration
            e = (capPropNew-capProp(iSnr,iStr))/capPropNew;
            capProp(iSnr,iStr) = capPropNew;
        end
    end
end

% plot figure 1
figure();
plot(snrVec,capSVD,'b.-',snrVec,capSMF,'ko-',snrVec,capProp,'rv-');
grid on;