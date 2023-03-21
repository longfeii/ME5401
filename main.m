clc;
clear;
a=3;b=2;c=2;d=4;
g=9.8;
Mf=2.14+c/20;
Mr=5.91-b/10;
Mc=1.74;
LFf=0.05;
Lr=0.128;
Lc=0.259;
Jx=0.5+(c-d)/100;
alpha=15.5-a/3+b/2;
gamma=11.5+(a-c)/(b+d+3);
Hf=0.18;
Hr=0.161;
Hc=0.098;
LF=0.133;
LR=0.308+(a-d)/100;
mux=3.33-b/20+a*c/60;
beta=27.5-d/2;
delta=60+(a-b)*c/10;
den=Mf*Hf^2+Mr*Hr^2+Mc*Hc^2+Jx;
a51=Mc*g/den;
a52=(Mf*Hf+Mr*Hr+Mc*Hc)*g/den;
a53=(Mr*Lr*LF+Mc*Lc*LF+Mf*LFf*LR)*g/((LR+LF*den));
a54=-Mc*Hc*alpha/den;
a55=-mux/den;
a56=Mf*Hf*LFf*gamma/den;
b51=Mc*Hc*beta/den;
b52=-Mf*Hf*LFf*delta/den;

A=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1
   0 6.5 -10 -alpha 0 0;
   a51 a52 a53 a54 a55 a56;
   5 -3.6 0 0 0 -gamma];
B=[0 0;
   0 0;
   0 0;
   beta 11.2;
   b51 b52;
   40 delta];
C=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0];
D=zeros(3,2);
t=linspace(0,10,100);
% u=[max(0,min(t-1,1));max(0,min(t-1,1))];
g=ss(A,B,C,zeros(3,2));
x0=[0.2 -0.1 0.15 -1 0.8 0]';
% step(g)%two input, three output six states
% figure;
% initial(g,x0)
[NUM1,DEN1]=ss2tf(A,B,C,D,1); %for the first input with three outputs
roots1=roots(DEN1);%dominant root is about-0.1 unstable point 2.8
[NUM2,DEN2]=ss2tf(A,B,C,D,2);%for the second input with three outputs
%% pole placement for feedback controller
%check the controbility
Q=ctrb(A,B);
rQ=rank(Q);
assert(rQ==6);
syms wn xi 
%satisfy the specific performance
f1=exp(-pi*xi/(1-xi^2)^0.5)-0.1<0;
f2=4/(xi*wn)-5<0;
[wn,xi]=solve([f1 f2]);%xi>0.59 let it be 0.707 while wn*xi>0.8
xi=0.75;
wn=1.2/xi;
%desired poles
deroot=roots([1 2*xi*wn wn^2]);
%new model's six poles
P=double([deroot(1),deroot(2),real(deroot(1))*2,real(deroot(1))*2,real(deroot(1))*3,real(deroot(1))*3,]');%Locate extra poles other than the dominant ones to be 2-5 times faster than the dominant ones:
%palce fucntion
kp=place(A,B,P);
% new state model of place
gp=ss(A-B*kp,B,C,zeros(3,2));
[numkp,denkp]=ss2tf(A-B*kp,B,C,D,1);
roots2=roots(denkp);%dominant root is -1.2
% [y,overshoot,st]=check(gp,t,x0);

% %determine the pole position
% for i=2:5
%     
%             P2=double([deroot(1),deroot(2),real(deroot(1))*i,real(deroot(1))*(i),real(deroot(1))*(i+1),real(deroot(1))*(i+1)]');
%             kp2=place(A,B,P2);
%             gp2=ss(A-B*kp2,B,C,D);
%             [y,overshoot,st]=check(gp2,t,x0);
%             ps(i,1:3)=st(1,:);ps(i,4:6)=st(2,:);
%             po(i,1:3)=overshoot(1,:);po(i,4:5)=overshoot(2,:);
% end

%own code
%unity rank method
% urs=zeros(100);
% uro=zeros(100);
% 
% for i=1:50
%     w=1;
%     for j=0:50
%     q=[i,j]';
%     Qunirank=ctrb(A ,B*q);
%     rQunirank=rank(Qunirank);
%     assert(rQunirank==6);
%     syms kur [1 6]
%     syms s
%     Aur=A-B*q*kur;
%     coef1=coeffs(det(s*eye(6)-Aur),s);
%     coef2=coeffs((s-P(1))*(s-P(2))*(s-P(3))*(s-P(4))*(s-P(5))*(s-P(6)),s);
%     kurs=solve(coef2==coef1);
%     kurd=double(struct2array(kurs));
%     Kur=q*kurd;
%     % new state model of unity rank method
%     gur=ss(A-B*Kur,B,C,zeros(3,2));
%     % [numkur,denkur]=ss2tf(A-B*Kur,B,C,D,1);
%     % rootsur=roots(denkur);
%     [y,overshoot,st]=check(gur,t,x0);
%     urs(i,w:(w+2))=st(1,:);urs(i,(w+3):(w+5))=st(2,:);
%     uro(i,w:(w+2))=overshoot(1,:);uro(i,(w+3):(w+5))=overshoot(2,:);
%     w=w+6;
%     end
% end

%full rank method
lli=Q(:,1:6);
Cfrr=rank(lli);%rank is 6 so the first 6th columns are independent
Cfr=[Q(:,1),Q(:,3),Q(:,5),Q(:,2),Q(:,4),Q(:,6)];
inCfr=inv(Cfr);
T=[inCfr(3,:);inCfr(3,:)*A;inCfr(3,:)*A^2;inCfr(6,:);inCfr(6,:)*A;inCfr(6,:)*A^2];
Abar=round(T*A/T,3);
Bbar=round(T*B,3);
syms kbar [2,6] 
syms x
% % one way of expression of Ad
Ad11=sym2poly((x-P(1))*(x-P(2))*(x-P(5)));
Ad1=[-Ad11(4),-Ad11(3),-Ad11(2),0,0,0] ;
Ad22=sym2poly((x-P(3))*(x-P(4))*(x-P(6)));
Ad2=[0,0,0,-Ad22(4),-Ad22(3),-Ad22(2)];
canonical_form=Abar-Bbar*kbar;
f1fr=canonical_form(3,:)-Ad1;
f2fr=canonical_form(6,:)-Ad2;
[kbar]=solve([f1fr,f2fr]);
kbard1=double(struct2array(kbar));
kbard=[kbard1(1:6);kbard1(7:12)];
kfr=kbard*T;
% new state model of full rank method
gfr=ss(A-B*kfr,B,eye(6),zeros(6,2));
[numkfr,denkfr]=ss2tf(A-B*kfr,B,C,D,1);
rootsfr=roots(denkfr)
% [y,overshoot,st]=check(gfr,t,x0);

%another way of expression of Ad
% Ada11=sym2poly((x-P(1))*(x-P(2))*(x-P(5))*(x-P(3))*(x-P(4))*(x-P(6)));
% Ada1=[-Ada11(7),-Ada11(6),-Ada11(5),-Ada11(4),-Ada11(3),-Ada11(2)] ;
% canonical_form=Abar-Bbar*kbar;
% f1fra=canonical_form(6,:)-Ada1;
% f2fra=canonical_form(3,:)-[0,0,0,1,0,0];
% [kbara]=solve([f1fra,f2fra]);
% kbarda1=double(struct2array(kbara));
% kbarda=[kbarda1(1:6);kbarda1(7:12)];
% kfra=kbarda*T
% % new state model of full rank method
% gfra=ss(A-B*kfra,B,C,zeros(3,2));
% [numkfr,denkfr]=ss2tf(A-B*kfra,B,C,D,1);
% rootsfr=roots(denkfr)
% [y,overshoot,st]=check(gfra,t,x0);

%non-zero initial state
% figure;
% initial(gfr,x0);
%% LQR method
% by the lqr 
overshoot=zeros(50,6);
st=zeros(50,6);
for x=1:50
Q=[x 0 0 0 0 0;
   0 x 0 0 0 0;
   0 0 x 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1];
R=[1 0;
   0 1];
N=zeros();
Hlqr=sqrt(Q);
Woh=obsv(A,Hlqr);
rWoh=rank(Woh);
assert(rWoh==6);
[klqr,Plqr,e] = lqr(g,Q,R,N);
glqr=ss(A-B*klqr,B,C,zeros(3,2));
% %check by step input
y=real(step(glqr,t));
%overshoot check

for i=1:3
    if y(end,i,1)>0
      overshoot(x,i)=(max(y(:,i,1))-y(end,i,1))/y(end,i,1);
    elseif y(end,i,1)<0
      overshoot(x,i)=(min(y(:,i,1))-y(end,i,1))/y(end,i,1);
    end
end
%settling time check

for i=1:3
    for j=1:100
    if abs(y(j,i,1)-y(end,i,1))/abs(y(end,i,1))<0.02
         st(x,i)=j/10;
         break
    end
    end
end
for i=1:3
    if y(end,i,2)>0
      overshoot(x,i+3)=(max(y(:,i,2))-y(end,i,2))/y(end,i,2);
    elseif y(end,i,2)<0
      overshoot(x,i+3)=(min(y(:,i,2))-y(end,i,2))/y(end,i,2);
    end
end
%settling time check

for i=1:3
    for j=1:100
    if abs(y(j,i,2)-y(end,i,2))/abs(y(end,i,2))<0.02
         st(x,i+3)=j/10;
         break
    end
    end
end

end

%after the determination of Q and R
Q=[3 0 0 0 0 0;
   0 3 0 0 0 0;
   0 0 3 0 0 0;
   0 0 0 4 0 0;
   0 0 0 0 0.5 0;
   0 0 0 0 0 1];
R=[1 0;
    0 1];
%by icare
[Pd,kicare]=icare(A,B,Q,R,[],[],[]);
%by own code
%solve ARE directly
% flqr=A'*Pdo+Pdo*A+Q-Pdo*B*inv(R)*B'*Pdo;
% Pdo=

%solve ARE by eigenvalue-eigenvector based algorithm
Gamma=[A,-B*inv(R)*B';
       -Q,-A'];
[V,ev]=eig(Gamma);
ev=sum(ev);
v=V(:,find(real(ev)<0));
P2=v(7:12,:)/v(1:6,:);
klqr2=real(inv(R)*B'*P2)
glqr2=ss(A-B*klqr2,B,C,zeros(3,2));
[numkfr,denkfr]=ss2tf(A-B*klqr2,B,C,D,1);
rootsfr=roots(denkfr)
[y,overshoot,st]=check(glqr2,t,x0);
%non-zero initial state
% figure;
% initial(glqr2,x0);
%% 3 state observer
%full order observers
figure;
% for i=3:5
Qo=ctrb(A',C');
rQo=rank(Qo);
assert(rQo==6);
%six poles for obsever 3times faster than controller poles
Pob=3.*P;
syms ko [3 6]
syms x
%full rank method
llio=Qo(:,1:6);
Cfrro=rank(llio);%rank is 6 so the first 6th columns are independent
Cfro=[Qo(:,1),Qo(:,4),Qo(:,2),Qo(:,5),Qo(:,3),Qo(:,6)];
inCfro=inv(Cfro);
To=[inCfro(2,:);inCfro(2,:)*A';inCfro(4,:);inCfro(4,:)*A';inCfro(6,:);inCfro(6,:)*A'];
Abaro=round(To*A'/To,3);
Bbaro=round(To*C',3);
% one way of expression of Ad
Ad11=sym2poly((x-Pob(1))*(x-Pob(2)));
Ad1=[-Ad11(3),-Ad11(2),0,0,0,0] ;
Ad22=sym2poly((x-P(3))*(x-P(4)));
Ad2=[0,0,-Ad22(3),-Ad22(2),0,0];
Ad33=sym2poly((x-P(5))*(x-P(6)));
Ad3=[0,0,0,0,-Ad33(3),-Ad33(2)];
canonical_form=Abaro-Bbaro*ko;
f1fr=canonical_form(2,:)-Ad1;
f2fr=canonical_form(4,:)-Ad2;
f3fr=canonical_form(6,:)-Ad3;
[ko]=solve([f1fr,f2fr,f3fr]);
kod1=double(struct2array(ko));
kod=[kod1(1:6);kod1(7:12);kod1(13:18)];
koo=kod*To;
L=koo'
% new state model of full rank method
gfro=ss([A-B*klqr2,B*klqr2;zeros(6),A-L*C],[B;zeros(6,2)],eye(12),zeros(12,2));
% [y,overshoot,st]=check(gfro,t,x0);

initial(gfro,[x0;x0]);
hold on
% legend(num2str(i))
% initial(gfr,x0)
% end
%reduced order observer
%dimension of y is 3 while A=6x6 so construct a three order observer
syms To [3 6]
Go=eye(3);
D1=[-3.6,0,0;
    0,-3.6,0;
    0,0,-4.8];
To=solve(To*A-D1*To-Go*C==0);
Td=double(struct2array(To));
Tr=[Td(1:6);Td(7:12);Td(13:18)];
Eo=Tr*B;
%check the rank
Aro=[C;Tr];
ror=rank(Aro);
Aroi=inv(Aro);
%x1=y1 x2=y2 x3=y3

%% decouplinig controller
C2=[1 0 0 0 0 0;
   0 0 1 0 0 0];

Bx=zeros(2);
SI=zeros(2,1);
for i=1:2
    si=1;
    while 1
         if C2(i,:)*(A^(si-1))*B==0
              si=si+1;
         else 
              Bx(i,:)=C2(i,:)*(A^(si-1))*B;
              SI(i)=si;
              break
         end
    end
end
% two si are all 2
F=inv(Bx);
%set the poles for d(t) and pi(t) to be -1.2 due to the performance
% (A-P(1)*eye(6))*(A-P(2)*eye(6))*(A-P(3)*eye(6))*(A-P(4)*eye(6))*(A-P(5)*eye(6))*(A-P(6)*eye(6))
Cxx=zeros(2,6);
for i=1:2
    Cxx(i,:)=C2(i,:)*(A^2+2.4*A+1.44*eye(6));
end
kdecoup=F*Cxx
gdecoup=ss(A-B*kdecoup,B*F,C2,zeros(2));
[numk1,denk1]=ss2tf(A,B,C2,zeros(2),1);% the first col for before transfer function
[numk2,denk2]=ss2tf(A,B,C2,zeros(2),2);
[numkd1,denkd1]=ss2tf(A-B*kdecoup,B*F,C2,zeros(2),1);%only related to u(1)
[numkd2,denkd2]=ss2tf(A-B*kdecoup,B*F,C2,zeros(2),2);
Gde=tf(gdecoup)
 
%check by step input
% figure;
% step(gdecoup);
% %non-zero initial state
% figure;
% initial(gdecoup,x0);
syms s x
%output control (only know the output)

Ns=[poly2sym(numk1(1,:)),poly2sym(numk2(1,:)); 
   poly2sym(numk1(2,:)),poly2sym(numk2(2,:))];% their denominators are sameGtf=tf(ss(A,B,C2,zeros(2,2)))
kd=adjoint(Ns);
xll=(det(Ns)*eye(2)/poly2sym(denk1(:)));
[n1,d1]=numden(xll(1,1));
n3=coeffs(n1);
% [z1,p1,k1]=tf2zp(coeffs(n1),coeffs(d1))
[n2,d2]=numden(xll(2,2));
% [z2,p2,k2]=tf2zp(coeffs(n2),coeffs(d2))
root11=roots(sym2poly(n1));
root22=roots(sym2poly(d1));
rootns=roots(sym2poly(det(Ns)));
deno=(x+1.2)^3;
ks=[1/deno,0;0,1/deno];
h1=[1/(x+1.2),0;0,1/(x+1.2)];
H1=xll*ks;
Hd=inv(eye(2)+H1)*H1;
% Ns=[25*s^4+414*s^3+494*s^2-7324*s+10047,11*s^4+182*s^3-238*s^2-5434*s+14267;
%      40*s^4+806*s^3+2168*s^2-11097*s-10767,60*s^4+1213*s^3+3166*s^2-17089*s-11458]
% kd=adjoint(Ns);
% ds=s^6+31.8*s^5+285.4*s^4+307.6*s^3-3462.3*s^2-1175.2*s-2030;
% GKd=det(Ns)*eye(2)/ds
% [n1,d1]=numden(GKd(1,1))
% [n2,d2]=numden(GKd(2,2))
% root11=roots(sym2poly(n1))
% root22=roots(sym2poly(d1))% unstable poles
%% special set point
%can onlly measure the three outputs
ysp=-0.1*C*inv(A)*B*[-0.5+(a-b)/20;0.1+(b-c/(a+d+10))];
%Agument the system
Qa5=[A,B;C(1:2,:),zeros(2,2)];
ras5=rank(Qa5);%the rank is equal to 8
Aa5=[A,zeros(6,2);-C(1:2,:),zeros(2)];
Ba5=[B;zeros(2,2)];
Ca5=[C(1:2,:),zeros(2)];
Q5=eye(8);
Q5=diag([1 1 1 1 1 1 3 3]);
R=[1 0;
   0 1];
%solve ARE by eigenvalue-eigenvector based algorithm
Gamma=[Aa5,-Ba5*inv(R)*Ba5';
       -Q5,-Aa5'];
[V,ev]=eig(Gamma);
ev=sum(ev);
v=V(:,find(real(ev)<0));
P2=v(9:16,:)/v(1:8,:);
ka5=real(inv(R)*Ba5'*P2)
Aa52=[A-B*ka5(:,1:6),-B*ka5(:,7:8);-C(1:2,:),zeros(2)];
t1 = 0:0.05:20;
u1 = zeros(1,length(t1));
u1(t1>=10) = 1;
u2=ones(1,length(t1));
u=[u1;u2];
w=[-1,1]';
Ba52=[B*w,zeros(6,1);0,ysp(1);0,ysp(2)];
Ba52c=[zeros(6,1);ysp(1);ysp(2)];
ga5=ss(Aa52,Ba52,Ca5,zeros(2,2));
% x0_s=[x0;ysp(1)-x0(1);ysp(2)-x0(2)];
lsim(ga5,u,t1,zeros(8,1))
ga5c=ss(Aa52,Ba52c,Ca5,zeros(2,1));
t=linspace(0,10,100);
[y,overshoot,st]=check2(ga5c,t,x0) ;


%for the third output
Qa53=[A,B;C(3,:),zeros(1,2)];
ras53=rank(Qa53)%the rank is equal to 7
Aa53=[A,zeros(6,1);-C(3,:),zeros(1)];
Ba53=[B;zeros(1,2)];
Ca53=[C(3,:),zeros(1)];
Q53=diag([1 1 1 1 1 1 3]);
%solve ARE by eigenvalue-eigenvector based algorithm
Gamma=[Aa53,-Ba53*inv(R)*Ba53';
       -Q53,-Aa53'];
[V,ev]=eig(Gamma);
ev=sum(ev);
v=V(:,find(real(ev)<0));
P2=v(8:14,:)/v(1:7,:);
ka53=real(inv(R)*Ba53'*P2)
Aa532=[A-B*ka53(:,1:6),-B*ka53(:,7);-C(3,:),zeros(1)];
t = 0:0.05:20;
u1 = zeros(1,length(t));
u1(t>=10) = 1;
u2=ones(1,length(t));
u=[u1;u2];
Ba532=[B*w,zeros(6,1);0,ysp(3)];
ga53=ss(Aa532,Ba532,Ca53,zeros(1,2));
figure;
lsim(ga53,u,t,zeros(7,1))
%% mutivariable integral control
%Augment the system
Qa6=[A,B;C,zeros(3,2)];
ras6=rank(Qa6')%the rank is equal to 8
% syms ka1 [2,6] 
% syms ka2 [2,2]
%for the first two outputs
Aa6=[A,zeros(6,2);-C(1:2,:),zeros(2,2)];
Ba6=[B;zeros(2,2)];
Ca6=[C(1:2,:),zeros(2,2)];
%solve ARE by eigenvalue-eigenvector based algorithm
Gamma=[Aa6,-Ba6*inv(R)*Ba6';
       -Q5,-Aa6'];
[V,ev]=eig(Gamma);
ev=sum(ev);
v=V(:,find(real(ev)<0));
P2=v(9:16,:)/v(1:8,:);
ka6=real(inv(R)*Ba6'*P2)
Aa62=[A-B*ka6(:,1:6),-B*ka6(:,7:8);-C(1:2,:),zeros(2,2)];
w=[-1,1]';
ga6=ss(Aa62,[B*w;1;-1],Ca6,zeros(2,1));%upper three rows are disturbence lower two rows are set points
step(ga6)
%for the third output
Qa63=[A,B;C(3,:),zeros(1,2)];
ras63=rank(Qa63);%the rank is equal to 7
Aa63=[A,zeros(6,1);-C(3,:),zeros(1)];
Ba63=[B;zeros(1,2)];
Ca63=[C(3,:),zeros(1)];
%solve ARE by eigenvalue-eigenvector based algorithm
Gamma=[Aa63,-Ba63*inv(R)*Ba63';
       -Q53,-Aa63'];
[V,ev]=eig(Gamma);
ev=sum(ev);
v=V(:,find(real(ev)<0));
P2=v(8:14,:)/v(1:7,:);
ka63=real(inv(R)*Ba63'*P2)
Aa632=[A-B*ka63(:,1:6),-B*ka63(:,7);-C(3,:),zeros(1)];
ga63=ss(Aa632,[B*w;1.5],Ca63,zeros(1));%upper three rows are disturbence lower two rows are set points
figure;
step(ga63)
%% PID control
% Gamma=[A-B*kdecoup,-B*F*inv(R)*(B*F)';
%        -Q,-(A-B*kdecoup)'];
% [V,ev]=eig(Gamma);
% ev=sum(ev);
% v=V(:,find(real(ev)<0));
% P2=v(7:12,:)/v(1:6,:);
% klqr2d=real(inv(R)*(B*F)'*P2)
% glqr2d=ss(A-B*kdecoup-B*F*klqr2d,B*F,C2,zeros(2,2));
% pole(glqr2d)
% [y,overshoot,st]=check2(glqr2d,t,x0);
C2=[1 0 0 0 0 0;
   0 0 1 0 0 0];
Bx=zeros(2);
SI=zeros(2,1);
for i=1:2
    si=1;
    while 1
         if C2(i,:)*((A-B*klqr2)^(si-1))*B==0
              si=si+1;
         else 
              Bx(i,:)=C2(i,:)*((A-B*klqr2)^(si-1))*B;
              SI(i)=si;
              break
         end
    end
end
% two si are all 2
F=inv(Bx);
%set the poles for d(t) and pi(t) to be -1.2 due to the performance
Cxx=zeros(2,6);
for i=1:2
    Cxx(i,:)=C2(i,:)*((A-B*klqr2+1.2*eye(6))^2);
end
kdecoup=F*Cxx;
gdecoup=ss(A-B*klqr2-B*kdecoup,B*F,C2,zeros(2));
pole(gdecoup)
tf(gdecoup)
pole(glqr2)