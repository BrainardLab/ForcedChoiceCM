function OLPlotConeEffects(primaryCones, testCones, subjectID, sessionNum, varargin)
% Plot cone effects of subjects' OL Rayleigh matches.
%
% Syntax:
%   OLPlotConeEffects(primaryCones, testCones, subjectID, sessionNum, varargin)
%
% Description
%    Takes in calculated cone responses to test and primary lights from
%    subject matches. For each match, produces a bar graph comparing the
%    response of each cone to the test and primary lights. Can also run
%    this with average cone responses to produce a bar graph with error
%    bars.
%
% Inputs:
%    primaryCones  - 3xn matrix with cone responses for primary lights,
%                    where n is the number of matches. L cone responses are
%                    stored in the first row, M in the second, and S in the
%                    third. 
%    testCones     - 3xn matrix with cone responses for test lights, where
%                    n is the number of matches. L cone responses are
%                    stored in the first row, M in the second, and S in the
%                    third. 
%
% Outputs:
%    Produces a plot comparing test and primary cone responses for each 
%    match. When run with averaged data, produces one plot. 
%
% Optional key-value pairs:
%    'measured'  - logical indicating whether the cone effects were
%                  calculated from radiometer measurements (true) or 
%                  nominal spd data (false). Default is false.
%    'average'   - logical indicating whether the provided cone responses
%                  represent individual matches (false) or averages (true).
%                  Default is false. 
%    'err'       - 1x6 vector containing standard error for each cone for 
%                  test and primary lights. Must input this if 'average' is
%                  true. 
%    'fType'    - character vector containing file type extension for 
%                 saving results. Default is 'tiff'.  

% History 
%    1/22/20   dce  -Modified program from OLTestConeEffects   
% Example: OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')


close all; 

% Parse input 
p = inputParser;
p.addParameter('measured', false, @(x) (islogical(x)));
p.addParameter('average', false, @(x) (islogical(x)));
p.addParameter('err', 0, @(x) (isnumeric(x))); 
p.addParameter('fType', 'tiff', @(x) (ischar(x)));
p.parse(varargin{:});

err = p.Results.err; 

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
if p.Results.measured
    category = 'Measured';
else
    category = 'Nominal';
end
if p.Results.average 
    theTitle = sprintf('Subject %s %s Average Cone Responses', subjectID, category);
else 
    theTitle = sprintf('Subject %s Session %g Match %g %s Cone Responses', subjectID, sessionNum, i, category);
end 
title(theTitle);

% If average, add error bars 
if p.Results.average
    errorbar(1:6, cones, err, err); 
end 

file = fullfile(outputDir, strrep(num2str(theTitle,' ', '_'));
saveas(gcf, file, p.Results.fType);
end 
end

