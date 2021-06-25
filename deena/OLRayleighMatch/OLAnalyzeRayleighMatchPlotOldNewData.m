function OLAnalyzeRayleighMatchPlotOldNewData(subjID,fileNum,firstFName,secondFName)
% Function that makes a plot of cone responses from two OneLight Rayleigh 
% matching sessions which were previously analyzed separately. Takes in 
% the names of two files produced by OLAnalyzeRayleighMatch, and creates
% and saves a plot which distinguishes between matches from the two files.
%
% Syntax:
%   OLAnalyzeRayleighMatchPlotOldNewData(subjID,fileNum,firstFName,secondFName)
%
% Inputs:
%    subjID       - Character vector for subject ID
%    fileNum      - Integer appended to the end of the name of the saved plot
%                   (to prevent overwriting of previous files)
%    firstFName   - Filepath (to OLAnalyzeRayleighMatch output file)
%    secondFName  - Filepath (to OLAnalyzeRayleighMatch output file)
%
% Outputs:
%    None (makes and saves a figure)
%
% Optional key-value pairs:
%    None

% History:
%   06/09/21  dce       Wrote it

% Load data
firstData = load(firstFName);
secondData = load(secondFName);
matchWls = unique([firstData.lightCombos;secondData.lightCombos],'rows');

% Define results directory
resDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),subjID);
if ~exist(resDir,'dir')
    mkdir(resDir);
end
resFile = fullfile(resDir,[subjID '_analyzeMultipleSessions_' num2str(fileNum) '.pdf']);

% Set up figure 
stdConeDiffPlot = figure();
hold on;
xlim([0 0.35]);
ylim([0.015 0.035]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title([subjID ' Cone Response Difference - Standard Cones'],'interpreter','none');

legendHandles = [];
legendEntries = {};
plotColors = 'rkbgcm';
plotColorNames = {'Red','Black','Blue','Green','Cyan','Magenta'};

% Loop through wavelength combos, and compute and plot excitations for each  
for i = 1:size(matchWls,1)
    % Extract relevant spds
    sessionIndsFirst = all(firstData.lightCombos==matchWls(i,:),2);
    sessionIndsSecond = all(secondData.lightCombos==matchWls(i,:),2);
    measPrimarySpdsTrialFirst = firstData.measPrimarySpds(:,sessionIndsFirst);
    measRefSpdsTrialFirst = firstData.measRefSpds(:,sessionIndsFirst);
    measPrimarySpdsTrialSecond = secondData.measPrimarySpds(:,sessionIndsSecond);
    measRefSpdsTrialSecond = secondData.measRefSpds(:,sessionIndsSecond);
    
    % Compute cone excitations
    primaryResStdFirst = firstData.stdObs.T_cones*measPrimarySpdsTrialFirst;
    refResStdFirst = firstData.stdObs.T_cones*measRefSpdsTrialFirst;
    primaryLMinusMStdFirst = primaryResStdFirst(1,:)-primaryResStdFirst(2,:);
    primaryLPlusMStdFirst = primaryResStdFirst(1,:)+primaryResStdFirst(2,:);
    refLMinusMStdFirst = refResStdFirst(1,:)-refResStdFirst(2,:);
    refLPlusMStdFirst = refResStdFirst(1,:)+refResStdFirst(2,:);
    
    primaryResStdSecond = secondData.stdObs.T_cones*measPrimarySpdsTrialSecond;
    refResStdSecond = secondData.stdObs.T_cones*measRefSpdsTrialSecond;
    primaryLMinusMStdSecond = primaryResStdSecond(1,:)-primaryResStdSecond(2,:);
    primaryLPlusMStdSecond = primaryResStdSecond(1,:)+primaryResStdSecond(2,:);
    refLMinusMStdSecond = refResStdSecond(1,:)-refResStdSecond(2,:);
    refLPlusMStdSecond = refResStdSecond(1,:)+refResStdSecond(2,:);
    
    % Plot cone excitations
    figure(stdConeDiffPlot);
    hold on;
    a  = plot(primaryLMinusMStdFirst./primaryLPlusMStdFirst,...
        primaryLPlusMStdFirst,[plotColors(i) '* ']);
    plot(refLMinusMStdFirst./refLPlusMStdFirst,refLPlusMStdFirst,[plotColors(i) 'o ']);
    
    plot(primaryLMinusMStdSecond./primaryLPlusMStdSecond,...
        primaryLPlusMStdSecond,[plotColors(i) '* ']);
    plot(refLMinusMStdSecond./refLPlusMStdSecond,refLPlusMStdSecond,[plotColors(i) 's '],...
        'MarkerFaceColor',plotColorNames{i});
    
    legendHandles = [legendHandles,a];
    legendEntries{end+1} = num2str(matchWls(i,3));
    
    % Add lines between primary mixture/reference light pairs
    for j = 1:length(refLMinusMStdFirst(1,:))
        plot([primaryLMinusMStdFirst(j)/primaryLPlusMStdFirst(j)...
            refLMinusMStdFirst(j)/refLPlusMStdFirst(j)],...
            [primaryLPlusMStdFirst(j) refLPlusMStdFirst(j)],'k-');
    end
    for j = 1:length(refLMinusMStdSecond(1,:))
        plot([primaryLMinusMStdSecond(j)/primaryLPlusMStdSecond(j)...
            refLMinusMStdSecond(j)/refLPlusMStdSecond(j)],...
            [primaryLPlusMStdSecond(j) refLPlusMStdSecond(j)],'k-');
    end
end
% Finish up figure, and save 
legend(legendHandles,legendEntries);
text(0.03,0.033,'square = second data set, circle = first data set');
NicePlot.exportFigToPDF(resFile,stdConeDiffPlot,300);
end