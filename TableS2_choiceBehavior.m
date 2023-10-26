%%% This script plots the choice data for each animal. It also displays
%%% numerical information that's helpful
%%%
%%% Adam Smoulder, 12/14/22

% subjectNames = {'Nelson','Ford','Earl','Prez','Rocky'};
subjectNames = {'Earl','Prez','Rocky'};
choiceFilenames = {...
    ...'FigS1ChoiceData_Nelson_20221213_142353.mat';
    ...'FigS1ChoiceData_Ford_20221213_142353.mat';
    'FigS1ChoiceData_Earl_20221213_142353.mat';
    'FigS1ChoiceData_Prez_20221214_084550.mat';
    'FigS1ChoiceData_Rocky_20221213_144336.mat';
    };
nsubjects = length(choiceFilenames);
figfolder = 'D:\AdamMatlab\~chokingUnderPressure\~forManuscript\';
addpath(genpath(figfolder)) % also has helper functions
dateString = grabDateTimeString;

nrewards = 4;

%% Just reduce it to % correct on a grayscale map

curColors = gray;
figure
for f = 1:nsubjects
    load(choiceFilenames{f})
    
    % Combine across directions and just get % correct
    n = triu(n_rewDirS_rewDirS,1)+triu(n_rewDirS_rewDirS',1);
    ncorr = triu(nS1_rewDirS_rewDirS,1)+triu((n_rewDirS_rewDirS-nS1_rewDirS_rewDirS)',1);
    valueGrid = 100*ncorr./(n+realmin);
    
    subplot(1,nsubjects,f)
    f
    imagesc(valueGrid)
    title(['Monkey ' subjectNames{f}(1)])
    set(gca,'clim',[0 100])
    colormap(curColors);
    colorbar
    xticks(1:nrewards); xticklabels({'S','M','L','J'})
    yticks(1:nrewards); yticklabels({'S','M','L','J'})
    set(gca,'fontsize',12,'fontname','arial','tickdir','out','box','off');
    
    for r1 = 1:nrewards
        for r2 = r1:nrewards
            curText = [num2str(ncorr(r1,r2)) '/' num2str(n(r1,r2))];
            text(r2,r1,curText,'fontsize',12,'fontname','arial','fontweight','bold','horizontalalignment','center','color','k')
        end; clear r2
    end; clear r1
    
    disp('Number selected non-bias (red)')
    disp(nS1_rewDirS_rewDirS)
    disp('Number total')
    disp(n_rewDirS_rewDirS)

    
    % Get ntotal, ncorrect, nwrong, and nwrongbias
    ntotal = sum(sum(n_rewDirS_rewDirS-diag(diag(n_rewDirS_rewDirS)))); % just sum of off diagonal elements
    nnonbias = sum(sum(triu(n_rewDirS_rewDirS,1))); % # correct in dir of bias
    nbias = sum(sum(triu(n_rewDirS_rewDirS',1)));
    assert(ntotal==nnonbias+nbias)
    ncorrectnonbias = sum(sum(triu(nS1_rewDirS_rewDirS,1))); % # correct in dir of bias
    ncorrectbias = sum(sum(triu(n_rewDirS_rewDirS',1))) - sum(sum(triu(nS1_rewDirS_rewDirS',1)));
    ncorrect = ncorrectbias+ncorrectnonbias;
    nwrongbias = nnonbias-ncorrectnonbias;
    nwrongnonbias = nbias-ncorrectbias;
    nwrong = ntotal-ncorrect;
    
    disp([num2str(round(100*ncorrect/ntotal,1)) '% correct (' num2str(ncorrect) '/' num2str(ntotal) ')'])
    disp([num2str(round(100*nwrongbias/nwrong,1)) '% errors in bias dir (' num2str(nwrongbias) '/' num2str(nwrong) ')'])

end; clear f
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) 1494 339])

% Save the figure
figname = ['TableS2Fig_choiceBehavior'];
saveFigAndSvg([figfolder 'recentFigs\'],figname);
saveFigAndSvg([figfolder 'allFigs\'],[figname '_' dateString]);

