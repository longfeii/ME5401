function [y,overshoot,st]=check2(g,t,x0)
%check by step input
step(g)
y=real(step(g,t));
%overshoot check
overshoot=zeros(1,2);

for i=1:2
    if y(end,i)>0
      overshoot(1,i)=(max(y(:,i))-y(end,i))/y(end,i);
    elseif y(end,i)<0
      overshoot(1,i)=(min(y(:,i))-y(end,i))/y(end,i);
%     elseif y(end,i)=0
%         overshoot(1,i)=(min(y(:,i))-y(end,i))/y(end,i);
    end
end

%settling time check
st=zeros(1,2);
for i=1:2
    for j=1:100
    if abs(y(j,i)-y(end,i))/abs(y(end,i))<0.02
         st(1,i)=j/10;
         break
    end
    end
end
end
%non-zero initial state
% figure;
% initial(g,x0);