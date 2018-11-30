function f=spread(F,X)
global h N xN yN dtheta;

c=dtheta/(h*h);
f=zeros(xN,yN,2);

for k=1:size(X,1)
  s=X(k,:)/h;
  i=floor(s);
  r=s-i;
  i1=mod((i(1)-1):(i(1)+2),xN)+1;
  i2=mod((i(2)-1):(i(2)+2),yN)+1;
  w=phi1(r(1)).*phi2(r(2));
  f(i1,i2,1)=f(i1,i2,1)+(c*F(k,1))*w;
  f(i1,i2,2)=f(i1,i2,2)+(c*F(k,2))*w;
end

