% ib2D.m
% This script is the main program.
% Cochlea

global dt Nb N h rho mu xip xim yip yim a xN yN;
global kp km dtheta K;

initialize_cochlea  %Initialize variables
init_a      % array a

% interp - interpolation with delta_h
% spread - f = integral F * delta -> f = sum F * delta_h
% fluid  - Use a, does fourier transform, abuses the periodic domain
% init_a - a is used in fourier transform: multiplying by sin

% Each time step
for clock=1:clockmax
  % Two step temporal discretization, Mixed Method
  
  XX=X+(dt/2)*interp(u,X,Nb);           % Take half time step? See what membrane does
  WW = Wall + (dt/2)*interp(u,Wall,Nw); % Target Points
  Oval(:,1) = 2 + A*sin(pi*Oval_y)*sin(omega*clock*dt);
  
  ff=spread(Force(XX),XX,Nb);         % New force at this half step
  FWall = -Krigid*(WW - Wall);          % Target Points
  ffwall = spread(FWall,WW,Nw);
  
  ForceOval = [-A*omega^2*sin(pi*Oval_y)*sin(omega*clock*dt),zeros(N+1,1)];
  ffoval = spread(ForceOval,Oval,N);
  
  ff = ff + ffwall + ffoval;
  [u,uu]=fluid(u,ff);              % Discretized Navier Stokes
  
  X(2:Nb,:) = X(2:Nb,:) + dt*interp(uu,XX(2:Nb,:),Nb-2);         % New membrane position
  %Wall = Wall + dt*interp(uu,WW,Nw); %Necessary?
  
  vorticity=(u(xip,:,2)-u(xim,:,2)-u(:,yip,1)+u(:,yim,1))/(2*h);
  
  dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
  values= (-10*dvorticity):dvorticity:(10*dvorticity);
  valminmax=[min(values),max(values)];
  
  %animation:
  
  contour(xgrid,ygrid,vorticity,values)
  plot(X(:,1),X(:,2),'ko')
  hold on
  plot(Wall(:,1),Wall(:,2), 'rs')
  plot(Oval(:,1),Oval(:,2),'bd')
  axis([0,40*L,0,4*L])
  caxis('auto')
  %axis equal
  axis manual
  drawnow
  hold off
  
end

