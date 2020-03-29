% Script for closing an annular stimulus used in OLRayleighMatch and other
% similar programs.
% History 
%            dce xx/xx/20   - Wrote program 
%            dce 3/29/20     - Edited for style

try
    % Close the window if it exists and is a struct. Otherwise, create a
    % new window and close it
    if exist('annulusData', 'var')
        if isstruct(annulusData) && isfield(annulusData, 'win');
            annulusData.win.close();
        else
            win = GLWindow('windowID', length(mglDescribeDisplays));
            win.close();
        end
    else
        win = GLWindow('windowID', length(mglDescribeDisplays));
        win.close();
    end
    
    % Catch errors and rethrow to the Matlab command window
catch e
    disp('An exception was raised');
    rethrow(e);
end