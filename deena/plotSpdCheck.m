% plots an array of target spds 
function plotSpdCheck(spds, cal)
[~, col] = size(spds)
figure(1);
hold on; 
for i = 1:col
    plot(SToWls(cal.describe.S),spds(:,i));
end
hold off; 
end 