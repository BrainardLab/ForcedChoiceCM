function OLPlotConeEffectsBackup(primaryCones, testCones, subjectID, sessionNum, varargin)
close all; 

% Parse input 
p = inputParser;
p.addParameter('measured', false, @(x) (islogical(x)));
p.addParameter('average', false, @(x) (islogical(x)));
p.addParameter('fType', 'tiff', @(x) (ischar(x)));
p.parse(varargin{:});

% Create folder for saving figures 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'), 'coneResponsePlots', subjectID, num2str(sessionNum));
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

[~, nMatches] = size(primaryCones); 
for i = 1:nMatches

% Make bar graph of cones
figure();
cones = [testCones(1,i), primaryCones(1,i); testCones(2,i),...
    primaryCones(2,i); testCones(3,i), primaryCones(3,i)];
bar(cones);

% Add axis labels 
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity'); 
legend('Test','Primaries');

% Set plot title 
if measured
    category = 'Measured';
else
    category = 'Nominal';
end
if average 
    theTitle = sprintf('Subject %s %s Average Cone Responses', subjectID, category);
else 
    theTitle = sprintf('Subject %s Match %g %s Cone Responses', subjectID, i, category);
end 
title(theTitle);


file = fullfile(outputDir, strrep(theTitle,' ', '_'));
saveas(gcf, file, p.Results.ftype);
end 
end
