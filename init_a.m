%init_a.m
%This script initializes the array  a  
%that is used in the fluid solver

global a;
a=zeros(xN,yN,2,2);
for m1=0:(xN-1)
  for m2=0:(yN-1)
    a(m1+1,m2+1,1,1)=1;
    a(m1+1,m2+1,2,2)=1;
  end
end

for m1=0:(xN-1)
  for m2=0:(yN-1)
    if ~(((m1==0)|(m1==xN/2)) & ((m2==0)|(m2==yN/2)))         % If not left corner, middle bot, middle left, middle
      t=[2*pi/xN ; 2*pi/yN].*[m1;m2];
      s=sin(t);                                             % sin(2*pi/L m1), sin(2*pi/L m2)
      ss=(s*s')/(s'*s);                                     % sin(2*pi/L m_i) * sin(2*pi/L m_j) / [ sin^2 (2*pi/L m1) + sin^2(2*pi/L m2) ]
%     a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)-(s*s')/(s'*s);
      a(m1+1,m2+1,1,1)=a(m1+1,m2+1,1,1)-ss(1,1);            % 1 - sin^2(2*pi/L m1)
      a(m1+1,m2+1,1,2)=a(m1+1,m2+1,1,2)-ss(1,2);            % - sin(2*pi/L m1) * sin(2*pi/L m2)
      a(m1+1,m2+1,2,1)=a(m1+1,m2+1,2,1)-ss(2,1);            % - sin(2*pi/L m1) * sin(2*pi/L m2)
      a(m1+1,m2+1,2,2)=a(m1+1,m2+1,2,2)-ss(2,2);            % 1 - sin^2(2*pi/L m2)
    end
  end
end

for m1=0:(xN-1)
  for m2=0:(yN-1)
    t=[pi/xN ; pi/yN].*[m1;m2];
    s=sin(t);
    a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)...
                    /(1+(dt/2)*(mu/rho)*(4/(h*h))*(s'*s));
  end
end

