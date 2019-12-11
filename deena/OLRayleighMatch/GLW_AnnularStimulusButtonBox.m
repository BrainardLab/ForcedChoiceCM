% Annular stimulus which can be moved and resized with a button box. The
% directional pad can be used to adjust the position of the annulus left,
% right, up, or down. The Y and A keys are used to increase/decrease the
% inner diameter, while the B and X keys are used to increase/decrease the
% outer diameter. The back button quits, while the start button toggles
% between coarse and fine adjustment. Based on GLWindow example code
%
% 12/3/13  npc Wrote it.
% 6/23/14  ms  Display an annulus.
% 11/21/19 dce Modified to take inp�qqq'''ut from button box instead of keypresses

% Background
backgroundRGB = [1 1 1]; % Half-on
annulusRGB = [0 0 0]; % Full-on

% The 'inner circle', i.e. the hole
innerCircleDiameter = 78; % px
innerCircleRGB = [1 1 1]; % Assume that it's the same as background, but it can be changed.

% The 'outer circle'
outerCircleDiameter = 303;%200; % px
outerCircleRGB = annulusRGB;

% Fixation cross
fixationCrossDiameter = 10;
fixationCrossRGB = [1 0 0];

% Define step sizes for the navigation
fineStepSize = 1;
coarseStepSize = 30;

% Get information about the displays attached to our system.
displayInfo = mglDescribeDisplays;

% We will present everything to the last display. Get its ID.
lastDisplay = length(displayInfo);

% Get the screen size
screenSizeInPixels = displayInfo(lastDisplay).screenSizePixel;

win = [];
try
    % Create a full-screen GLWindow object
    win = GLWindow( 'SceneDimensions', screenSizeInPixels, ...
        'BackgroundColor', backgroundRGB,...
        'windowID',        lastDisplay);
    
    % Open the window
    win.open;
    
    % Add stimulus image to the GLWindow
    centerPosition = [0 13];
    win.addOval(centerPosition, [outerCircleDiameter outerCircleDiameter], outerCircleRGB, 'Name', 'outerCircle');
    win.addOval(centerPosition, [innerCircleDiameter innerCircleDiameter], innerCircleRGB, 'Name', 'innerCircle');
    win.addOval(centerPosition, [fixationCrossDiameter fixationCrossDiameter], fixationCrossRGB, 'Name', 'fixationCross');
    
    % Render the scene
    win.draw;
    
    % Initialize button box
    gamePad = GamePad();
    
    % Display some information
    disp('Press Back to exit');
    disp('Commands:');
    disp(['Press Start to toggle step size: coarse is ' num2str(coarseStepSize)...
        ' px, fine is ' num2str(fineStepSize) ' px']);
    disp('Use directional pad for positional adjustment (L, R, U, D)');
    disp('Use A (-) and Y (+) for diameter adjustment [inner circle]');
    disp('Use X (-) and B (+) for diameter adjustment [outer circle]');
    
    keepLooping = true;
    currCoarseSteps = true; % Start with coarse steps
    while (keepLooping)
        key = gamePad.getKeyEvent();
        if (~isempty(key))
            if currCoarseSteps
                currStepSize = coarseStepSize;
            else
                currStepSize = fineStepSize;
            end
            switch(key.charCode)
                case 'GP:Back'  % Quit
                    keepLooping = false;
                case 'GP:Start' % Change adjustment step size
                    currCoarseSteps = ~currCoarseSteps;
                case 'GP:East'  % Move annulus right
                    centerPosition(1) = centerPosition(1) + currStepSize;
                case 'GP:West'  % Move annulus left
                    centerPosition(1) = centerPosition(1) - currStepSize;
                case 'GP:North' % Move annulus up
                    centerPosition(2) = centerPosition(2) + currStepSize;
                case 'GP:South' % Move annulus down
                    centerPosition(2) = centerPosition(2) - currStepSize;
                case 'GP:X'     % Decrease outer diameter
                    outerCircleDiameter = outerCircleDiameter - currStepSize;
                case 'GP:B'     % Increase outer diameter
                    outerCircleDiameter = outerCircleDiameter + currStepSize;
                case 'GP:A'     % Decrease inner diameter
                    innerCircleDiameter = innerCircleDiameter - currStepSize;
                case 'GP:Y'     % Increase inner diameter
                    innerCircleDiameter = innerCircleDiameter + currStepSize;
            end
            % Check that diameters are not zero
            if innerCircleDiameter == 0
                innerCircleDiameter = 1;
            end
            
            if outerCircleDiameter == 0
                outerCircleDiameter = 1;
            end
            
            % Update the position and dimensions of the circles
            win.setObjectProperty('outerCircle', 'Center', centerPosition);
            win.setObjectProperty('innerCircle', 'Center', centerPosition);
            win.setObjectProperty('fixationCross', 'Center', centerPosition);
            win.setObjectProperty('outerCircle', 'Dimensions', [outerCircleDiameter outerCircleDiameter]);
            win.setObjectProperty('innerCircle', 'Dimensions', [innerCircleDiameter innerCircleDiameter]);
            
            % Print out the position of the annulus
            fprintf('Outer d: %g, Inner d: %g, x: %g, y: %g\n', outerCircleDiameter, innerCircleDiameter, centerPosition(1), centerPosition(2));
            
            % Draw
            win.draw;
        end
    end
    win.draw;
    
    % Optional setting to save annular stimulus window so it can be used
    % later
    saveWin = false; 
    if saveWin
        file = fullfile(getpref('ForcedChoiceCM', 'rayleighDataDir'), 'OLAnnulusSettings.mat');
        save(file, 'win');
    end 
    
catch e
    disp('An exception was raised');
    
    % Send the error back to the Matlab command window.
    rethrow(e);
end




