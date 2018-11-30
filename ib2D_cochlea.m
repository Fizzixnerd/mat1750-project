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
  
  XXm=Xm+(dt/2)*interp(u,Xm);           % Take half time step? See what membrane does
  XXwall = Xwall + (dt/2)*interp(u,Xwall); % Target Points
  XXoval = Xoval + (dt/2)*interp(u,Xoval);
  XXround = Xround + (dt/2)*interp(u,Xround);
  
  Oval(:,1) = 2 + A*sin(pi*Oval_y)*sin(omega*clock*dt);
  
  Fm=Force(XXm);
  Fm(1,:) = Krigid*([2,2]-Xm(1,:));
  Fm(end,:) = Krigid*([37,2]-Xm(end,:));
  ff = spread(Fm,XXm) + spread(Krigid*(Wall - XXwall),XXwall) + spread(Krigid*(Oval - XXoval),XXoval) ...
                                                            + spread(Krigid*(Round - XXround),XXround);         % New force at this half step
  
  [u,uu]=fluid(u,ff);              % Discretized Navier Stokes
  
  Xm = Xm + dt*interp(uu,XXm);         % New membrane position
  Xwall = Xwall + dt*interp(uu,XXwall);
  Xoval = Xoval + dt*interp(uu,XXoval);
  %Wall = Wall + dt*interp(uu,WW,Nw); %Necessary?
  
  vorticity=(u(xip,:,2)-u(xim,:,2)-u(:,yip,1)+u(:,yim,1))/(2*h);
  
  dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
  values= (-10*dvorticity):dvorticity:(10*dvorticity);
  valminmax=[min(values),max(values)];
  
  %animation:
  if mod(clock, 1) == 0
      contour(xgrid,ygrid,vorticity,values)
      plot(Xm(:,1),Xm(:,2),'ko')
      hold on
      plot(Xwall(:,1),Xwall(:,2), 'rs')
      plot(Xround(:,1),Xround(:,2),'rs')
      plot(Xoval(:,1),Xoval(:,2),'bd')
      %plot(Oval(:,1),Oval(:,2),'bd')
      axis([0,40*L,0,4*L])
      caxis('auto')
      %axis equal
      axis manual
      drawnow
      hold off
  end
end

