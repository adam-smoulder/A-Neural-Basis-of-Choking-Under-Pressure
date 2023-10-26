function [signalVarMean, signalVarUpperErrorBar, signalVarLowerErrorBar,...
    noiseVarMean, noiseVarUpperErrorBar, noiseVarLowerErrorBar,...
    SNRMean, SNRUpperErrorBar, SNRLowerErrorBar] = snrBarPlot(sigVar,noiseVar,snr,Nr,showSig)
% snrBarPlot takes in six matrices of assumed size [N x R x L] where N is
% the number of samples calculated in a bootstrapped distribution, R is 
% the number of reward conditions (or whatever should be compared closely 
% within of these subplot), and L is the number of latents (what should be
% split across subplots). Inputs are in the order shown:
% 
% Inputs: 
% sigVar - signal variance, meaning the variance of the means across
% classes (e.g. variance of mean firing rates for 8 different reach
% directions)
% noiseVar - noise variance, meaning the mean of variances across classes
% (e.g. take the variance for each of the 8 directions individually, then
% take the mean of those)
% snr - signal to noise ratio (signal var divided by noise var - input
% separately because it's assumed to have been calculated across
% bootstraps)
% Nr - minimum number of trials for a given reward size across directions
% (so how many pts go into the bootstraps), should be R x 1
% showSig - true/false on whether or not to show significance stars
%
% Outputs:  
% sigVarMean - nrewards x nlatents mean signal variances
% sigVarUpperErrorBar - the upper error bar on the 95% CI for signal
% variance. This is not the raw location of the top, but how much above the
% mean it should be.
% sigVarLowerErrorBar - the lower error bar on the 95% CI for the signal
% variance. 
% noiseVarMean, noiseVarUpperErrorBar, noiseVarLowerErrorBar, SNRMean,
% SNRUpperErrorBar, SNRLowerErrorBar are similar to these^
%
% Adam Smoulder, 3/11/19 (edited 3/25/19)

alphaDiv2 = 0.025;
BC = 1; % Bonferonni correction

tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution

figure;

[nboots, nrewards, nlatents] = size(sigVar);
X = reshape((1:nrewards)'+(0:(nrewards+1):((nrewards+1)*(nlatents-1))),[nrewards*nlatents,1]);

% evaluate nonoverlapping confidence intervals within each latent
signalVarMean = squeeze(mean(sigVar)); % R x L
sig2p5per = squeeze(prctile(sigVar,alphaDiv2/BC*100));
sig97p5per = squeeze(prctile(sigVar,(1-alphaDiv2/BC)*100));
noiseVarMean = squeeze(mean(noiseVar));
noise2p5per = squeeze(prctile(noiseVar,alphaDiv2/BC*100));
noise97p5per = squeeze(prctile(noiseVar,(1-alphaDiv2/BC)*100));
SNRMean = squeeze(mean(snr));
snr2p5per = squeeze(prctile(snr,alphaDiv2/BC*100));
snr97p5per = squeeze(prctile(snr,(1-alphaDiv2/BC)*100));
allSigRels = [];
allNoiseRels = [];
allSNRRels = [];
allSigVarPs = [];
allNoiseVarPs = [];
allSNRPs = [];
for i = 1:nlatents
    sigrels = [];
    noiserels = [];
    snrrels = [];
    for j = 1:(nrewards-1)
        for k = (j+1):nrewards
            sigVarDiff = squeeze(sigVar(:,j,i)-sigVar(:,k,i));
            noiseVarDiff = squeeze(noiseVar(:,j,i)-noiseVar(:,k,i));
            snrDiff = squeeze(snr(:,j,i)-snr(:,k,i));
            
            sigVarDenom = sqrt(squeeze(var(sigVar(:,j,i))+var(sigVar(:,k,i))));
            noiseVarDenom = sqrt(squeeze(var(noiseVar(:,j,i))+var(noiseVar(:,k,i))));
            snrDenom = sqrt(squeeze(var(snr(:,j,i))+var(snr(:,k,i))));

            sigVarP = tdist2T(mean(sigVarDiff)/sigVarDenom,Nr(j)+Nr(k)-1);
            noiseVarP = tdist2T(mean(noiseVarDiff)/noiseVarDenom,Nr(j)+Nr(k)-1);
            snrP = tdist2T(mean(snrDiff)/snrDenom,Nr(j)+Nr(k)-1);
            
%             if (sum(sigVarDiff > 0)/length(sigVarDiff) > 0.95) || (sum(sigVarDiff < 0)/length(sigVarDiff) > 0.95)
            if sigVarP <= alphaDiv2/BC || sigVarP >= (1-alphaDiv2/BC)
                sigrels = [sigrels, {[j,k]+(nrewards+1)*(i-1)}];
                allSigVarPs(end+1) = sigVarP;
            end
%             if (sum(noiseVarDiff > 0)/length(noiseVarDiff) > 0.95) || (sum(noiseVarDiff < 0)/length(noiseVarDiff) > 0.95)
            if noiseVarP <= alphaDiv2/BC || noiseVarP >= (1-alphaDiv2/BC)
                noiserels = [noiserels, {[j,k]+(nrewards+1)*(i-1)}];
                allNoiseVarPs(end+1) = noiseVarP;
            end
%             if (sum(snrDiff > 0)/length(snrDiff) > 0.95) || (sum(snrDiff < 0)/length(snrDiff) > 0.95)
            if snrP <= alphaDiv2/BC || snrP >= (1-alphaDiv2/BC)
                snrrels = [snrrels, {[j,k]+(nrewards+1)*(i-1)}];
                allSNRPs(end+1) = snrP;
            end
        end; clear k
    end; clear j
    allSigRels = [allSigRels sigrels];
    allNoiseRels = [allNoiseRels noiserels];
    allSNRRels = [allSNRRels snrrels];
end; clear i

signalVarUpperErrorBar = sig97p5per-signalVarMean;
signalVarLowerErrorBar = sig2p5per-signalVarMean;
noiseVarUpperErrorBar = noise97p5per-noiseVarMean;
noiseVarLowerErrorBar = noise2p5per-noiseVarMean;
SNRUpperErrorBar = snr97p5per-SNRMean;
SNRLowerErrorBar = snr2p5per-SNRMean;

subplot(3,1,1)
hold on
bar(X,signalVarMean(:),1,'facecolor',[0.5 0.75 1])
errorbar(X,signalVarMean(:),signalVarLowerErrorBar(:),signalVarUpperErrorBar(:),'k.','linewidth',1)
if ~isempty(allSigRels) && showSig
%     sigstar(allSigRels,0.05.*ones(length(allSigRels),1));
    sigstar(allSigRels, ternop(allSigVarPs>0.5,1-allSigVarPs,allSigVarPs));
end
%xticks([])
xticks(mean(1:nrewards)+(0:(nrewards+1):((nrewards+1)*(nlatents-1))))
xticklabels(num2str((1:nlatents)'));
ylabel('Signal Var.')
set(gca,'fontsize',18)

subplot(3,1,2)
hold on
bar(X,noiseVarMean(:),1,'facecolor',[1,0.55,0.5])
errorbar(X,noiseVarMean(:),noiseVarLowerErrorBar(:),noiseVarUpperErrorBar(:),'k.','linewidth',1)
if ~isempty(allNoiseRels) && showSig
%     sigstar(allNoiseRels,0.05.*ones(length(allNoiseRels),1));
    sigstar(allNoiseRels, ternop(allNoiseVarPs>0.5,1-allNoiseVarPs,allNoiseVarPs));
end
%xticks([])
xticks(mean(1:nrewards)+(0:(nrewards+1):((nrewards+1)*(nlatents-1))))
xticklabels(num2str((1:nlatents)'));
ylabel('Noise Var.')
set(gca,'fontsize',18)

subplot(3,1,3)
hold on
bar(X,SNRMean(:),1,'facecolor',[0.75 0.75 0.75])
errorbar(X,SNRMean(:),SNRLowerErrorBar(:),SNRUpperErrorBar(:),'k.','linewidth',1)
if ~isempty(allSNRRels) && showSig
%     sigstar(allSNRRels,0.05.*ones(length(allSNRRels),1));
    sigstar(allSNRRels, ternop(allSNRPs>0.5,1-allSNRPs,allSNRPs));
end
xticks(mean(1:nrewards)+(0:(nrewards+1):((nrewards+1)*(nlatents-1))))
xticklabels(num2str((1:nlatents)'));
xlabel('Latent #')
ylabel('SNR')
set(gca,'fontsize',18)

disp('WATCH OUT FOR BUG - if values are super super close to equal, they will be "significantly" equal (if p > 0.95)')
disp('This should be very obvious (e.g. means are near equal, obviously super overlapping CIs) - way too long to fix, so edit it out yourself')
end

