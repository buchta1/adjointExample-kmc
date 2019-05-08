%%%%%%%%%%%%%%spring-mass-damper
clear all;
close all;
folderForPics='data'
nt=2000; %since we're derive-then-discretize v. discretize-then-derive you need large time resolution
T=25;
dt=T/(nt-1);

%data storage
y = zeros(nt,1);
yobs = zeros(nt,1);
y=zeros(nt,1);
d2ydt2=zeros(nt,1);
dydt=zeros(nt,1);
djdm=zeros(nt,1);
djdc=zeros(nt,1);
djdk=zeros(nt,1);
m=zeros(nt,1);
k=zeros(nt,1);
c=zeros(nt,1);
ki=zeros(nt,1);
ci=zeros(nt,1);
mi=zeros(nt,1);
v = zeros(nt,1);
vobs = zeros(nt,1); %this is just a placeholder for generating dummy spring data
t = zeros(nt,1);
insCost = zeros(nt,1);
data=zeros(nt,4);

%observation data generating m-k-c values
%m(:)=1.;
%k(:)=1.;
%c(:)=0.;
for i=1:nt
m(i)=1.*(1.-0.5*double(i)/double(nt)); 
k(i)=1.; 
c(i)=0; 
end;

%intentionally wrong m-k-c for demonstration
for i=1:nt
mi(i)=2; 
ki(i)=1.; 
ci(i)=1.; 
end;

%time vector for plotting purposes later
for i=1:nt
    t(i)=(i-1)*dt;
end;

yobs(1)=3.3;
y(1)=yobs(1);
vobs(1)=1.;
v(1)=vobs(1);

[success,vobs,yobs] = getForward(vobs,yobs,k,m,c,nt,dt);
[success,v,y] = getForward(v,y,ki,mi,ci,nt,dt);
[Cost]=getCost(y,yobs,nt,0.,T);

subplot(2,1,1);
plot(t,yobs,t,y,'o')
ax = gca;
ax.FontSize = 20; 
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('$y$','FontSize',20,'Interpreter','latex')
legend('$y_\mathrm{obs}$','$y_o$','Interpreter','latex')

%now run the adjoint solve
yadj=zeros(nt,1); %this includes the homogeneous IC
vadj=zeros(nt,1);
[success,vadj,yadj] = getAdjoint(vadj,yadj,ki,mi,ci,y,yobs,nt,dt);

subplot(2,1,2);
plot(t,yadj,t,vadj)
ax = gca;
ax.FontSize = 20; 
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('$y^\dagger$,$v^\dagger$','FontSize',20,'Interpreter','latex')
legend('$y^\dagger$','$v^\dagger$','Interpreter','latex')
saveas(gcf,[pwd '/data/initialComparisonAndAdjoint.png'])

data(:,1)=t(:);
data(:,2)=y(:);
data(:,3)=yobs(:);
data(:,4)=yadj(:);
%csvwrite('/Users/davidbuchta/Downloads/springAdjoint',data);

%compute the sensitivity
[d2ydt2] = getSecondDerivative(y,nt,dt);
[dydt] = getFirstDerivative(y,nt,dt);
for i=1:nt
    djdm(i)=-yadj(i)*d2ydt2(i);
    djdc(i)=-yadj(i)*v(i); 
    djdk(i)=-yadj(i)*y(i);
end;


%integrate the local sensitivity
[sensK]=integrate(0.,T,nt,djdk)
[sensC]=integrate(0.,T,nt,djdc)
[sensM]=integrate(0.,T,nt,djdm)
[sensitivity]=integrate(0.,T,nt,djdm.*djdm+djdc.*djdc+djdk.*djdk)

%brute force check the sensitivity using gradient 
[oldCost]=getCost(y,yobs,nt,0.,T)
Nerr=100;
err=zeros(Nerr,2);
for i=1:Nerr
alpha=1e-3/(double(i)^5);
[success,v,y] = getForward(v,y,ki(:)-alpha*djdk(:),mi(:)-alpha*djdm(:),ci(:)-alpha*djdc(:),nt,dt);
[Cost]=getCost(y,yobs,nt,0.,T);
dJdk=(Cost-oldCost)/alpha;
err(i,1)=alpha;
err(i,2)=abs(dJdk+sensitivity)/sensitivity;
end;


figure(2)
loglog(err(:,1),err(:,2))%,err(:,1),err(:,3),err(:,1),err(:,4))
ax = gca;
ax.FontSize = 20; 
xlabel('$\alpha$','FontSize',20,'Interpreter','latex')
ylabel('$\mathcal{E}/||\mathcal{G}||$','FontSize',20,'Interpreter','latex')
legend('error')
saveas(gcf,[pwd '/data/adjointError.png'])


%Use optimize function
knew=zeros(nt,1);
cnew=zeros(nt,1);
mnew=zeros(nt,1);
oldy=zeros(nt,1);
oldv=zeros(nt,1);
newy=yobs(1);
newv=vobs(1);
oldy=yobs(1);
oldv=vobs(1);
[Ji,knew,mnew,cnew] = optimize(v,y,yobs,ki,mi,ci,nt,dt,T,1000);
[success,oldv,oldy] = getForward(oldv,oldy,ki,mi,ci,nt,dt);
[success,newv,newy] = getForward(newv,newy,knew,mnew,cnew,nt,dt);
figure(3)
p=plot(t,yobs,'o',t,oldy,'-',t,newy)
p(2).LineWidth = 1;
p(3).LineWidth = 3;
ax = gca;
ax.FontSize = 20; 
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('$y$','FontSize',20,'Interpreter','latex')
legend('observation','initial guess','optimized')
saveas(gcf,[pwd '/data/firstVersusOptimized.png'])


figure(4)
subplot(3,1,1);
plot(t,mnew,t,m,'o')
ax = gca;
ax.FontSize = 20; 
ylim([0 3.])
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('mass','FontSize',20,'Interpreter','latex')
legend('$m_\mathrm{opt}$','unknown $m$','Interpreter','latex')

subplot(3,1,2);
plot(t,knew,t,k,'o')
ax = gca;
ax.FontSize = 20; 
ylim([0 3.])
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('spring constant','FontSize',20,'Interpreter','latex')
legend('$k_\mathrm{opt}$','unknown $k$','Interpreter','latex')

subplot(3,1,3);
plot(t,cnew,t,c,'o')
ylim([-0.5 2.])
ax = gca;
ax.FontSize = 20; 
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('damper','FontSize',20,'Interpreter','latex')
legend('$c_\mathrm{opt}$','unknown $c$','Interpreter','latex')
set(gcf,'Position',[10 10 750 900])
saveas(gcf,[pwd '/data/k-m-cVersusUnknown.png'])

%naive steepest descent 
%todo: improve the optimizer (Brent's line minimization; conjugate gradients, etc) 
function [Ji,ki,mi,ci] = optimize(v,y,yobs,ko,mo,co,nt,dt,T,nIters)

tmpk=zeros(nt,1);
tmpc=zeros(nt,1);
tmpm=zeros(nt,1);
yadj=zeros(nt,1); 
vadj=zeros(nt,1);
djdm=zeros(nt,1);
djdc=zeros(nt,1);
djdk=zeros(nt,1);
d2ydt2=zeros(nt,1);

tmpk=ko;
tmpc=co;
tmpm=mo;

CostDuringOpt=zeros(nIters,3);
yDuringOpt=zeros(nt,7);
for t=1:nt
yDuringOpt(t,1)=dt*double(t-1); %set the time vector
end


for n=1:nIters
[success,v,y] = getForward(v,y,tmpk,tmpm,tmpc,nt,dt);
[Cost]=getCost(y,yobs,nt,0.,T); 
initialCost=Cost;
if n==1
    firstCost=Cost;
end

%save intermediate optimizations for plotting
if n==1 
 yDuringOpt(:,2)=y(:);, 
end
if n==10 
 yDuringOpt(:,3)=y(:);
end
if n==30 
 yDuringOpt(:,4)=y(:);
end
if n==100 
 yDuringOpt(:,5)=y(:);
end
if n==300 
 yDuringOpt(:,6)=y(:);
end
if n==1000 
 yDuringOpt(:,7)=y(:);
end

[success,vadj,yadj] = getAdjoint(vadj,yadj,tmpk,tmpm,tmpc,y,yobs,nt,dt);

%compute the sensitivity
[d2ydt2] = getSecondDerivative(y,nt,dt);
for i=1:nt
    djdm(i)=-yadj(i)*d2ydt2(i);
    djdc(i)=-yadj(i)*v(i); 
    djdk(i)=-yadj(i)*y(i);
end;
[sensitivity]=integrate(0.,T,nt,djdm.*djdm+djdc.*djdc+djdk.*djdk);

%naive steepest descent (keep the 
for i=1:50
 alpha=5.*sensitivity^0.5/(double(i)^4); %useful starting point for this data minimization 
 [success,v,y] = getForward(v,y,tmpk(:)-alpha*djdk(:),tmpm(:)-alpha*djdm(:),tmpc(:)-alpha*djdc(:),nt,dt);
 [Cost]=getCost(y,yobs,nt,0.,T);
 if Cost < initialCost %save k-m-c if cost reduced
  ki(:)=tmpk(:)-alpha*djdk(:);
  ci(:)=tmpc(:)-alpha*djdc(:);
  mi(:)=tmpm(:)-alpha*djdm(:); 
  initialCost=Cost;
  Ji=Cost;
 end
end
 tmpk(:)=ki(:);
 tmpc(:)=ci(:);
 tmpm(:)=mi(:);
 fprintf('(Iter,Cost,Relative)=(%2.0f,%8.5f,%8.5f)\n',double(n),Ji,Ji/firstCost*100.)
 CostDuringOpt(n,1)=n;
 CostDuringOpt(n,2)=Ji/firstCost;
 CostDuringOpt(n,3)=sensitivity;
  
end
figure(7)
loglog(CostDuringOpt(:,1),CostDuringOpt(:,2),'o',CostDuringOpt(:,1),CostDuringOpt(:,3)/CostDuringOpt(1,3),'o')
legend('$\mathcal{J}_i/\mathcal{J}_o$','$||\mathcal{G}_i||/||\mathcal{G}_o||$','Interpreter','latex')
ax = gca;
ax.FontSize = 20; 
xlabel('$i$ iteration','FontSize',20,'Interpreter','latex')
ylabel('Normalized Cost & Gradient','FontSize',20,'Interpreter','latex')
saveas(gcf,[pwd '/data/CostDuringOpt.png'])

figure(8)
plot(yDuringOpt(:,1),yobs(:),'o',...
    yDuringOpt(:,1),yDuringOpt(:,2),...
    yDuringOpt(:,1),yDuringOpt(:,3),...
    yDuringOpt(:,1),yDuringOpt(:,4),...
    yDuringOpt(:,1),yDuringOpt(:,5),...
    yDuringOpt(:,1),yDuringOpt(:,6),...
    yDuringOpt(:,1),yDuringOpt(:,7))
ax = gca;
ax.FontSize = 20; 
xlabel('$t$','FontSize',20,'Interpreter','latex')
ylabel('$y$','FontSize',20,'Interpreter','latex')
legend('$y_\mathrm{obs}$',...
       '$y_{i=0}$',...
       '$y_{i=10}$',...
       '$y_{i=30}$',...
       '$y_{i=100}$',...
       '$y_{i=300}$',...
       '$y_{i=1000}$',...
       'Interpreter','latex')
saveas(gcf,[pwd '/data/YDuringOpt.png'])

end

function [Cost] = getCost(y,yobs,nt,tmin,tmax)
insCost = zeros(nt,1);
for i=1:nt
    insCost(i)=(y(i)-yobs(i))*(y(i)-yobs(i));
end;
[Cost]=integrate(tmin,tmax,nt,insCost);
end

function [integral] = integrate(tmin,tmax,nt,scalarVec)
dt=(tmax-tmin)/(nt-1);
integral=0.;
integral=integral+dt*0.5*(scalarVec(1)+scalarVec(nt));
for i=2:nt-1
integral=integral+dt*0.5*(2.*scalarVec(i));
end;
end

function [success,vadj,yadj] = getAdjoint(vadj,yadj,k,m,c,y,yobs,nt,dt)
%set the adjoint initial condition at t=T (for our purposes always
%homogeneous)
vadj(nt)=0.;
yadj(nt)=0.;

%initialize storage for the adjoint non-constant coefficients
d2mdt2=zeros(nt,1);
dmdt=zeros(nt,1);
dcdt=zeros(nt,1);
[d2mdt2]=getSecondDerivative(m,nt,dt);
[dmdt]=getFirstDerivative(m,nt,dt);
[dcdt]=getFirstDerivative(c,nt,dt);

for i=nt-1:-1:1 %marching backward in time since adjoint IVP starts at t=T
    
%using second-order Runge-Kutta approximation
k1vadj=-dt*(... 
       -vadj(i+1)/m(i+1)*(2.*dmdt(i+1)-c(i+1)) ... 
       -yadj(i+1)/m(i+1)*(d2mdt2(i+1)-dcdt(i+1)+k(i+1)) ...
       +2.*(y(i+1)-yobs(i+1))/m(i+1));
k1yadj=-dt*vadj(i+1);

%use mid-point approximation for source terms and non-const. coefficients
%we're using constant dt so this is appropriate
ymid=(y(i+1)+y(i))*0.5;
yobsmid=(yobs(i+1)+yobs(i))*0.5;
kmid=(k(i+1)+k(i))*0.5;
cmid=(c(i+1)+c(i))*0.5;
mmid=(m(i+1)+m(i))*0.5;
dmdtmid=(dmdt(i+1)+dmdt(i))*0.5;
d2mdt2mid=(d2mdt2(i+1)+d2mdt2(i))*0.5;
dcdtmid=(dcdt(i+1)+dcdt(i))*0.5;

%finish the rk step
k2vadj=-dt*(... 
       -(vadj(i+1)+k1vadj*0.5)/mmid*(2.*dmdtmid-cmid) ... 
       -(yadj(i+1)+k1yadj*0.5)/mmid*(d2mdt2mid-dcdtmid+kmid) ...
       +2.*(ymid-yobsmid)/mmid);
k2yadj=-dt*(vadj(i+1)+k1vadj*0.5);

%finish the rk2
%recall we're integrating from t=T..0 hence the (i+1) on the right side
vadj(i)=vadj(i+1)+k2vadj;
yadj(i)=yadj(i+1)+k2yadj;

end;
success=1; %Todo at error catching in case doesn't finish
end

%assumes the IC is already set outside of this function
function [success,v,y] = getForward(v,y,k,m,c,nt,dt)
for i=2:nt
k1v=dt*(-k(i-1)*y(i-1)/m(i-1)-c(i-1)*v(i-1)/m(i-1));
k1y=dt*v(i-1);

kmid=(k(i)+k(i-1))*0.5;
cmid=(c(i)+c(i-1))*0.5;
mmid=(m(i)+m(i-1))*0.5;

k2v=dt*(-kmid*(y(i-1)+k1y*0.5)/mmid-cmid*(v(i-1)+k1v*0.5)/mmid);
k2y=dt*(v(i-1)+k1v*0.5);
v(i)=v(i-1)+k2v;
y(i)=y(i-1)+k2y;

end;
success=1.; %Todo at error catching in case doesn't finish
end

function[dfdt] = getFirstDerivative(f,nt,dt)
dfdt(1)=(-f(1)+f(2))/dt;
dfdt(nt)=(f(nt)-f(nt-1))/dt;
for i=2:nt-1
 dfdt(i)=0.5*(-f(i-1)+f(i+1))/dt;
end;
end

function[d2fdt2] = getSecondDerivative(f,nt,dt)
d2fdt2(1)=(f(1)-2.*f(2)+f(3))/dt/dt;
d2fdt2(nt)=(f(nt)-2.*f(nt-1)+f(nt-2))/dt/dt;
for i=2:nt-1
 d2fdt2(i)=(f(i-1)-2.*f(i)+f(i+1))/dt/dt;
end;
end



