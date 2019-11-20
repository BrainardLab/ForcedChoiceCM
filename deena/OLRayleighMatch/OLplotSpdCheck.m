% plots an array of target spds 
function OLplotSpdCheck(spds, cal)
[~, col] = size(spds);
hold on; 
for i = 1:col
    plot(SToWls(cal.describe.S),spds(:,i));
end
hold off; 
end 