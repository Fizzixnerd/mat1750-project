%initialize_cochlea.m

% Rectangular Domain, Periodic, Cochlea is inside within the boundaries
% Can use target points to suppress movement outside the cochlea
% Keeping domain periodic will be beneficial so we don't have to rewrite
% all the code
% 40 x 4 box, cochlea is 2 from left, 2 from right, 1 from top, 1 from bot

% Sound wave: Standing wave x = 2 + 2A sin(4 pi y) cos( omega t)
% omega angular frequency, 2pi times the frequency in hertz

%%% THINGS WE NEED TO PLAY AROUND WITH %%%%%%%%%

omega = 100;          % Frequency - Resemble Real life
A = 0.2;              % Amplitude
K_0  = 6e5;           % Stiffness original: 6e5
Krigid = 5e8;         % Stiffness of Walls
tmax=1     ;          % End time
dt=0.00001    ;        % Time step - In seconds
N  = 8   ;            % Number of points in 1 mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L  = 1.0 ;                      % Length 1mm = 0.1 cm
xN = 40*N    ;               % Number of points x axis
yN = 4*N    ;                % Number of points y axis
h  = L/N   ;                    % partition

xip =[(2:xN),1] ;              % x coord Plus one
xim =[xN,(1:(xN-1))];          % x coord Minus one
yip =[(2:yN),1];               % y coord Plus one
yim =[yN,(1:(yN-1))];          % y coord Minus one

Nb =  35*N  ;                  % Number of basilar membrane points
dtheta = h  ;                 % dtheta
kp = (3:Nb+1) ;               % Plus one theta step (ONLY FOR MIDDLE)
km = (1:Nb-1) ;                 % Minus one theta     (ONLY FOR MIDDLE)

% Stiffness
K  = K_0 *exp(-1.4*h*(1:Nb-1)/35);
K = K';

Nw1 = 35*N;                 % Top wall
Nw2 = 3*N;                  % Helicotrema wall
Nw3 = 35*N;                 % Bottom wall
Nw  = Nw1 + Nw2 + Nw3;      % Total will be Nw + 1

rho=1     ;                   % Density
mu=0.02   ;                   % Viscosity

Oval_y = 2 + (0:N)*h;
Oval_y = Oval_y';

Round_y = 1 + (0:N)*h;
Round_y = Round_y';

clockmax=ceil(tmax/dt);


% Immersed Boundary Basilar Membrane
Xm = zeros(Nb+1,2);
for k=0:Nb
  theta=k*dtheta;
  Xm(k+1,1)= 2*L + theta;
  Xm(k+1,2)= 2*L;
end

% Fixed Points
Wall = zeros(Nw,2);
for k = 1:Nw1
    theta = k*h;
    Wall(k,1) = 2*L + theta;
    Wall(k,2) = 3*L - theta/70;
end
for k = 1:Nw2
    theta = pi*k*h/3;
    Wall(Nw1 + k,1) = 37*L + sin(theta);
    Wall(Nw1 + k,2) = 2*L + 0.5*cos(theta);
end
for k = 1:Nw3
    theta = k*h;
    Wall(Nw1 + Nw2 + k,1) = 37*L - theta;
    Wall(Nw1 + Nw2 + k,2) = 1.5*L - theta/70;
end
Wall = [2 , 3 ; Wall];

Xwall = Wall;

% Oval window
Oval = zeros(N+1,2);
Oval(:,1) = 2;
Oval(:,2) = Oval_y;

Xoval = Oval;

% Round window
Round = zeros(N+1,2);
Round(:,1) = 2;
Round(:,2) = Round_y;

Xround = Round;

%%%%%%%%%%%%%%%%%%%%%%%%%%

u=zeros(xN,yN,2);

%Vorticies
vorticity=(u(xip,:,2)-u(xim,:,2)-u(:,yip,1)+u(:,yim,1))/(2*h);
dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;

%Lines
values= (-10*dvorticity):dvorticity:(10*dvorticity);
valminmax=[min(values),max(values)];

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
contour(xgrid,ygrid,vorticity,values)
hold on
plot(Xm(:,1),Xm(:,2),'ko')
plot(Xwall(:,1),Xwall(:,2),'rs')
plot(Xround(:,1),Xround(:,2),'rs')
plot(Xoval(:,1),Xoval(:,2),'bd')
axis([0,40*L,0,4*L])
%caxis(valminmax)
%axis equal
axis manual
drawnow
hold off

