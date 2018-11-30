%initialize.m
L=1.0                       % Length 1mm
N=64                        % Number of points in 1 mm
h=L/N                       % partition
ip=[(2:N),1]                % Plus one
im=[N,(1:(N-1))]            % Minus one
Nb=ceil(pi*(L/2)/(h/2))     % Number of boundary points
dtheta=2*pi/Nb              % dtheta
kp=[(2:Nb),1]               % Plus one theta step
km=[Nb,(1:(Nb-1))]          % Minus one theta
K=1                         % Stiffness
rho=1                       % Density
mu=0.01                     % Viscosity
tmax=4                      % End time
dt=0.01                     % Time step
clockmax=ceil(tmax/dt)


% Immersed Boundary Circle
for k=0:(Nb-1)
  theta=k*dtheta;
  X(k+1,1)=(L/2)+(L/4)*cos(theta);
  X(k+1,2)=(L/2)+(L/4)*sin(theta);
end

%initial Velocity (Vertical component has sin wave sheet
u=zeros(N,N,2);
for j1=0:(N-1)
  x=j1*h;
  u(j1+1,:,2)=sin(2*pi*x/L);
end

%Vorticies
vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;

%Lines
values= (-10*dvorticity):dvorticity:(10*dvorticity);
valminmax=[min(values),max(values)];

%Grid
xgrid=zeros(N,N);
ygrid=zeros(N,N);
for j=0:(N-1)
  xgrid(j+1,:)=j*h;
  ygrid(:,j+1)=j*h;
end

% Graph
set(gcf,'double','on')
contour(xgrid,ygrid,vorticity,values)
hold on
plot(X(:,1),X(:,2),'ko')
axis([0,L,0,L])
caxis(valminmax)
axis equal
axis manual
drawnow
hold off

