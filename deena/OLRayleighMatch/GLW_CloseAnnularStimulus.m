try
    if exist('annulusData', 'var')
        annulusData.win.close();
    else
        win = GLWindow('windowID', length(mglDescribeDisplays));
        win.close();
    end
catch e
    disp('An exception was raised');
    
    % Send the error back to the Matlab command window.
    rethrow(e);
end