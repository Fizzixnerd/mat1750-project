function w=laplacian(u)
global xim xip yip xim yim h;
w=(u(xip,:,:)+u(xim,:,:)+u(:,yip,:)+u(:,yim,:)-4*u)/(h*h);

% im = [N,1:(N-1)] = circular version of i-1
% ip = [2:N,1]     = circular version of i+1
% N  = number of points in each space direction

