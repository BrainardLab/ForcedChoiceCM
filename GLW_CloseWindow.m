win = [];
try
    win = GLWindow('windowID', length(mglDescribeDisplays));
    win.close;
catch e
    disp('An exception was raised');
    
    % Close the window if it was succesfully created.
    if ~isempty(win)
        win.close;
    end
    
    % Send the error back to the Matlab command window.
    rethrow(e);
end