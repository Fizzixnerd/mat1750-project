%initialize_cochlea.m

% Rectangular Domain, Periodic, Cochlea is inside within the boundaries
% Can use target points to suppress movement outside the cochlea
% Keeping domain periodic will be beneficial so we don't have to rewrite
% all the code
% 40 x 4 box, cochlea is 2 from left, 2 from right, 1 from top, 1 from bot

L  = 1.0 ;                      % Length 1mm = 0.1 cm
N  = 8   ;                      % Number of points in 1 mm
xN = 40*N    ;               % Number of points x axis
yN = 4*N    ;                % Number of points y axis
h  =L/N   ;                    % partition
ip =[(2:xN),1] ;              % Plus one
im =[xN,(1:(xN-1))];          % Minus one

Nb =  35*N ;                % Number of basilar membrane points
dtheta = h  ;                 % dtheta
kp =[(2:Nb),Nb + 1] ;         % Plus one theta step
km =[0,(1:(Nb-1))] ;           % Minus one theta
K  = 6e5*exp( - 1.4*(1:Nb) ) ;  % Stiffness

Nw1 = 35*N;                 % Top wall
Nw2 = 3*N;                  % Helicotrema wall
Nw3 = 35*N;                 % Bottom wall
Nw  = Nw1 + Nw2 + Nw3;      % Total

rho=1     ;                   % Density
mu=0.02   ;                   % Viscosity

tmax=4     ;                 % End time
dt=0.01    ;                 % Time step
clockmax=ceil(tmax/dt);


% Immersed Boundary Circle
X = zeros(Nb,2);
for k=0:(Nb-1)
  theta=k*dtheta;
  X(k+1,1)= 2*L + theta;
  X(k+1,2)= 2*L;
end

% Fixed Points
Wall = zeros(Nw,2);
for k = 1:Nw1
    theta = k*h;
    Wall(k,1) = 2*L + theta;
    Wall(k,2) = 3*L;
end
for k = 1:Nw2
    theta = pi*k*h/3;
    Wall(Nw1 + k,1) = 37*L + sin(theta);
    Wall(Nw1 + k,2) = 2*L + cos(theta);
end
for k = 1:Nw3
    theta = k*h;
    Wall(Nw1 + Nw2 + k,1) = 37*L - theta;
    Wall(Nw1 + Nw2 + k,2) = 1*L;
end
%initial Velocity
u=zeros(xN,yN,2);
for j1=0:(N-1)
  x=j1*h;
  u(j1+1,:,2)=sin(2*pi*x/L);
end

%Vorticies
%vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
%dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;

%Lines
%values= (-10*dvorticity):dvorticity:(10*dvorticity);
%valminmax=[min(values),max(values)];

%Grid
xgrid=zeros(xN,yN);
ygrid=zeros(xN,yN);
for j=0:(xN-1)
  xgrid(j+1,:)=j*h;
end
for j=0:(yN-1)
    ygrid(:,j+1)=j*h;
end

% Graph
set(gcf,'double','on')
%contour(xgrid,ygrid,vorticity,values)
hold on
plot(X(:,1),X(:,2),'ko')
plot(Wall(:,1),Wall(:,2),'rs')
axis([0,40*L,0,4*L])
caxis(valminmax)
axis equal
axis manual
drawnow
hold off

