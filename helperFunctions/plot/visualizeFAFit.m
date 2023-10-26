function visualizeFAFit(estMdl, lEstO, varargin)
DAY_DIVISIONS = []; % hacked by AS
warnOpts(assignOpts(varargin)); 

figure('Position', [10 10 1600 900]);
nlatents = size(estMdl.C,2);

meanPlotSlots = 1:4:(floor(nlatents/3)*4+1);

% Plot Means
subplot(nlatents,4, meanPlotSlots)
hold on;
plot(estMdl.d, 'r*'); 
title('Estimated Means'); 

% Plot C*C' values

ccPlotSlots = (max(meanPlotSlots)+4):4:(floor(nlatents*2/3)*4+1);

subplot(nlatents,4,ccPlotSlots)
estCC = estMdl.C*estMdl.C'; 
imagesc(estCC); 
axis equal
colorbar;
title('Est CC'''); 


% Examine fits of psi 
psiPlotSlots = (max(ccPlotSlots)+4):4:(floor(nlatents)*4);
subplot(nlatents,4,psiPlotSlots)

hold on;
plot(diag(estMdl.psi), 'r*'); 
title('Estimated Psi'); 

% Plot reconstruction of latents

latentsToPlot = size(lEstO,2);
%latentsToPlot = 10; % just show top 10
for dI = 1:latentsToPlot
    latentPlotSlots = (2+(dI-1)*4):(dI*4);
    a = subplot(latentsToPlot,4,latentPlotSlots); 
    plot(a, lEstO(:, dI), 'b.');
    hold on; 
    if dI==1, title('Reconstructed Latents', 'Interpreter','none');  end
    if dI==latentsToPlot, xlabel('Trial #'); end
    ylabel(['Lat. ' num2str(dI)])
    if ~isempty(DAY_DIVISIONS)
        for j = 1:(length(DAY_DIVISIONS)-1)
            plot(DAY_DIVISIONS(j)*ones(1,100),linspace(min(lEstO(:,dI)),max(lEstO(:,dI)),100),'r--','linewidth',2) % hacked by AS
        end; clear j
    end
    axis([1,size(lEstO,1),-inf,inf])
    %axis([1,40000,-inf,inf])
    hold off
end