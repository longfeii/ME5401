function [y,overshoot,st]=check(g,t,x0)
%check by step input
step(g)
y=real(step(g,t));
%overshoot check
overshoot=zeros(2,3);
for j=1:2
for i=1:3
    if y(end,i,j)>0
      overshoot(j,i)=(max(y(:,i,j))-y(end,i,j))/y(end,i,j);
    elseif y(end,i,j)<0
      overshoot(j,i)=(min(y(:,i,j))-y(end,i,j))/y(end,i,j);
%     elseif y(end,i)=0
%         overshoot(1,i)=(min(y(:,i))-y(end,i))/y(end,i);
    end
end
end
%settling time check
st=zeros(2,3);
for m=1:2
for i=1:3
    for j=1:100
    if abs(y(j,i,m)-y(end,i,m))/abs(y(end,i,m))<0.02
         st(m,i)=j/10;
         break
    end
    end
end
end
%non-zero initial state
% figure;
% initial(g,x0);
end