function OLPlotConeEffects(primaryCones, testCones, subjectID, sessionNum, varargin)
% Plot cone effects of subjects' OL Rayleigh matches.
%
% Syntax:
%   OLPlotConeEffects(primaryCones, testCones, subjectID, sessionNum, varargin)
%
% Description
%    Takes in calculated cone responses to test and primary lights from
%    subject matches (these are calculated with the OLGetConeEffects function).
%    For each match, produces a bar graph comparing the response of each
%    cone to the test and primary lights. Can also run this with average
%    cone responses to produce a bar graph with error bars.
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
%                  test and primary lights. Only displayed if 'average' is
%                  true. Default is zeros(1,6).
%    'fType'     - character vector containing file type extension for
%                  saving results. Default is 'tiff'.

% History
%    1/22/20   dce  - Modified program from OLTestConeEffects
%    4/5/20    dce  - Edited for stype
%
% Example: OLPlotConeEffects(p, t, 'Deena', 10)

% Parse input
p = inputParser;
p.addParameter('measured', false, @(x) (islogical(x)));
p.addParameter('average', false, @(x) (islogical(x)));
p.addParameter('err', zeros(1,6), @(x) (isnumeric(x)));
p.addParameter('fType', 'tiff', @(x) (ischar(x)));
p.parse(varargin{:});

err = p.Results.err;

% Create folder for saving figures, if it does not already exist. Average
% plots are saved in a folder named subjectID/average, while standard plots
% are saved in a folder named subjectId/sessionNum.
if p.Results.average
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
        'coneResponsePlots', subjectID, 'average');
else
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
        'coneResponsePlots', subjectID, num2str(sessionNum));
end
% Create the output directory if it does not exist
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

[~, nMatches] = size(primaryCones); % Number of matches

% Loop through each match and make a bar graph of cones.
for i = 1:nMatches
    % Put cones in order for plotting
    cones = [testCones(1,i), primaryCones(1,i); testCones(2,i),...
        primaryCones(2,i); testCones(3,i), primaryCones(3,i)];
    figure();
    bar(cones);
    hold on;
    
    % Add axis labels
    names ={'L'; 'M'; 'S' };
    set(gca,'xticklabel', names)
    ylabel('Relative Response Intensity');
    
    % If making an average plot, add error bars based on passed-in standard
    % error
    if p.Results.average
        % Vector list of the cones, in the order they are plotted
        errCones = [testCones(1,i), primaryCones(1,i), testCones(2,i),...
            primaryCones(2,i), testCones(3,i), primaryCones(3,i)];
        % Vector of positions for error bars
        errBarPos = [0.85 1.15 1.85 2.15 2.85 3.15];
        % Add error bars
        errorbar(errBarPos, errCones, err, err, 'k. ');
    end 
    
    % Set plot title
    if p.Results.measured
        category = 'Measured';
    else
        category = 'Nominal';
    end 
    if p.Results.average
        % A 0 session number indicates that the average was calculated
        % across sessions 
        if sessionNum == 0   
            theTitle = sprintf('Subject %s %s Average Cone Responses',...
                subjectID, category);
            legend('Test','Primaries', 'Standard Error');
        else
            theTitle = sprintf('Subject %s Session %g %s Average Cone Responses',...
                subjectID, sessionNum, category);
            legend('Test','Primaries', 'Standard Error');
        end
    else
        theTitle = sprintf('Subject %s Session %g Match %g %s Cone Responses',...
            subjectID, sessionNum, i, category);
        legend('Test','Primaries');
    end
end
title(theTitle);

% Save file 
fName = fullfile(outputDir, strrep(theTitle,' ', '_'));
saveas(gcf, fName, p.Results.fType);
end

