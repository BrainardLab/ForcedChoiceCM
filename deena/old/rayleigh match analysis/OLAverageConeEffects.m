function OLAverageConeEffects(subjectID, varargin)
% Calculates and plots average cone responses across trials and files
%
% Syntax:
%   OLAverageConeEffects(subjectID, varargin)
%
% Description
%    Takes in a subject ID and reads all of their Rayleigh matches across
%    trials from various OneLight experiments. Uses the OLGetConeEffects 
%    function to find cone responses for each match, then averages these 
%    cone responses and calculates the mean and standard error for each of 
%    three cones for primary and test lights. These are then plotted using 
%    the OLPlotConeEffects function. The user can also specify a session 
%    number, in which case only matches from that session will be averaged 
%    and plotted.
%
% Inputs:
%    subjectID      - character vector of subject ID
%
% Outputs:
%    Produces a plot comparing test and primary cone responses, which is 
%    saved in a subject-specific directory. 
%
% Optional key-value pairs:
%    'sessionNum'  - integer representing the session number. When 0 is
%                    entered, the program averages across all trials, but
%                    when a different number is entered, the program
%                    averages across that specific trial. Default is 0.
%    'measured'    - logical indicating whether to calculate cone effects 
%                    from radiometer measurements (true) or
%                    nominal spd data (false). Default is false.

% History
%    1/22/20   dce  - Modified program from OLTestConeEffects
%    1/29/20   dce  - Finished writing program 
%    4/5/20    dce  - Modified to incorporate dichromat/trichromat files 

% Examples: OLAverageConeEffects('DCE');
%           OLAverageConeEffects('DCE', 'sessionNum', 5); 

% Parse input 
close all;
p = inputParser;
p.addParameter('sessionNum', 0, @(x) (isnumeric(x)));
p.addParameter('measured', false, @(x) (islogical(x)));
p.parse(varargin{:});
sessionNum = p.Results.sessionNum;

% Arrays for storing aggregate data. These will have three columns for L,
% M, and S cones.
primaryCones = double.empty(0,3);
testCones = double.empty(0,3);

% Find relevant files for subjects
baseDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjectID,...
    subjectID);
if (sessionNum ~= 0)  
    % Find a specific session file to average over. This includes dichromat 
    % and trichromat experiment files, which have the format 
    % 'subjID_dichromat_sessionNum' or 'subjID_trichromat_sessionNum' 
    files = dir([baseDir(1:end-1), '*_' num2str(sessionNum), '.mat']);
else
    % Find all of the subject's files. This will automatically include
    % dichromat and trichromat results 
    files = dir([baseDir, '_*.mat']);
end

% Fill array with cone responses from each file
[row, ~] = size(files);
for i = 1:row
    fName = fullfile(files(i).folder, files(i).name);
    % Skip over measurement files, as these are not passed into
    % OLGetConeEffects to calculate cone responses 
    if (all(fName(end-7:end-4) == 'meas')) 
        continue;
    end
    % Calculate cone responses for each file, either measured or nominal
    [pCones, tCones] = OLGetConeEffects(fName, 'measured', p.Results.measured);
    primaryCones = [primaryCones; pCones'];
    testCones = [testCones; tCones'];
end

% Find mean and standard deviation for primary and test cone responses 
[numMatches, ~] = size(primaryCones); % Number of calculated cone responses

primaryAverages = mean(primaryCones, 1);
primarySD = std(primaryCones, 1);

testAverages = mean(testCones, 1);
testSD = std(testCones, 1);

% Calculate standard error 
err = [testSD(1), primarySD(1), testSD(2), primarySD(2), testSD(3),...
    primarySD(3)] / sqrt(numMatches);

% Make a plot of the average cone effects 
OLPlotConeEffects(primaryAverages', testAverages', subjectID, sessionNum,...
    'average', true, 'measured', p.Results.measured, 'err', err);
end