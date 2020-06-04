% Counts backwards from the end of a vector to find the first element 
% that's different from the end. 
function val2 = findLastUnique(vec)
val1 = vec(end); 
val2 = NaN; 
ind = length(vec) - 1;
stillLooping = true;
while(stillLooping && ind > 0)
    v2 = vec(ind); 
    if v2 ~= val1
        stillLooping = false; 
        val2 = v2; 
    else 
        ind = ind - 1; 
    end 
end 
end 