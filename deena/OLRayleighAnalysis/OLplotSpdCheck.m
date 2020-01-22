% plots an array of target spds

% de modified on 11/26
function OLplotSpdCheck(wls, spds)
[~, col] = size(spds);
hold on; 
for i = 1:col
    plot(wls,spds(:,i));
end
hold off; 
end 