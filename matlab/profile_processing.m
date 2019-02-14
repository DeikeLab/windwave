% This is for ploting the interface curve at different dumping time
path = '~/wave/ak/005/';
cc = jet(16);
eta = cell(1,16); % color curves according to time

for i = 1:16
    str = sprintf('graph%d.dat',i);
    filename = [path,str];
    %eta{1,i} = readinprofile (filename);
    A = readinprofile (filename);
    A_sort = sortrows(A);  % sort the order of point with increasing x (1st column)
    plot (A_sort(:,1), A_sort(:,2), 'LineWidth', 2, 'color', cc(i,:));
    hold on
end

legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16');