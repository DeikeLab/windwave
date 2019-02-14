% This is for plotting the velocity profile interpolated onto a n by n grid

path = '~/wave/ak/005/';
n = 129; % the resolution n+1 

N = 16; % the number og dumps
str = sprintf('fieldnl%d',6); % the name of the field file + #
filename = [path,str];
field = readinfield (filename, 2, inf); 
str = sprintf('eta%d',6); % the name of the interface file + #
filename = [path,str];
A = readinprofile (filename);
A_sort = sortrows(A);  % sort the order of point with increasing x (1st column)

for i = 1:N
    x((1+(i-1)*n):(i*n),1) = field((1+n*(i-1)*(n-1)/N):(n*(i-1)*(n-1)/N+n),1);
    y((1+(i-1)*n):(i*n),1) = field((1+n*(i-1)*(n-1)/N):(n*(i-1)*(n-1)/N+n),2);
    f((1+(i-1)*n):(i*n),1) = field((1+n*(i-1)*(n-1)/N):(n*(i-1)*(n-1)/N+n),3);
    u((1+(i-1)*n):(i*n),1) = field((1+n*(i-1)*(n-1)/N):(n*(i-1)*(n-1)/N+n),4);
end 
    
quiver (x, y, (1-f).*u, zeros(size(x)));
hold on 
plot (A_sort(:,1), A_sort(:,2));