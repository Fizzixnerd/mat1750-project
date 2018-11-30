function F=Force(X)
global kp km dtheta K Nb;

F = K.*(X(kp,:)+X(km,:)-2*X(2:Nb,:))/(dtheta*dtheta);
F = [0 , 0 ; F ; 0, 0];