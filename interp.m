function U=interp(u,X)
global h;
global N xN yN;
U=zeros(size(X,1),2);

%U(1,1) = 2;
%U(Nb+1,1) = 37;
%U(1,2) = 2;
%U(Nb+1,2) = 2;

for k=1:size(X,1)
  s=X(k,:)/h;                           % Approx where X is at which (integer) point
  i=floor(s);                           % Gives the integer right below X
  r=s-i;                                % Diff between X and integer below it
  i1=mod((i(1)-1):(i(1)+2),xN)+1;
  i2=mod((i(2)-1):(i(2)+2),yN)+1;
  w=phi1(r(1)).*phi2(r(2));             % delta_h(x - X) / h^2
  U(k,1)=sum(sum(w.*u(i1,i2,1)));       % sum over all u(x)*delta_h, only counting the u's that matter
  U(k,2)=sum(sum(w.*u(i1,i2,2)));
end

