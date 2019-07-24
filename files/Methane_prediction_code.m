% Matlab code for model fitting, validation and phase plane analysis in the manuscript: Second-generation stoichiometric mathematical model to
%predict methane emissions from oil sands tailings
%% fitting model to Decane, Octane, Heptane, Hexane  biodegradation data from laboratory cullture
function [vals,params,ci]=fit_to_n_alkanes_data
%dbstop if caught error
data = dataset('xlsfile', 'datatariq', 'Sheet',17); %the data can be obatined from Siddique et al. 2006
time=data.Day(1:7);
T = [time; time; time; time];%time in days
C61=data.C6mmole;%hexane data
C71=data.C7mmole;%heptane data
C81=data.C8mmole;%ocatne data
C101=data.C10mmole;%decane data
C6=C61(1:7);
C7=C71(1:7);
C8=C81(1:7);
C10=C101(1:7);
Y=[C6; C7; C8; C10];
dsid = [ones(length(time),1); 2*ones(length(time),1); 3*ones(length(time),1); 4*ones(length(time),1)];
X = [T dsid];
%params0=[327.615787410708,0.262218847536511,430.327575377768,269.14490458692,90.5884032453207,12.4460621918727,26,60,70,5,0.000573104641468882];%
params0=rand(1,11);
model=@(params,X)subfun(params, X);% @(params) use to indicate taht the function depends only on params variable and others are fixed
lb=[0,0,0,0,0,0,2,50,30,2, 0];
ub=inf*ones(1,length(params0));
[params,resnorm]= lsqcurvefit(model,params0,X,Y,lb,ub);
model1=@(vals,X)subfun1(vals, X);
opts = statset('nlinfit');
pts.RobustWgtFun = 'bisquare';
[vals,r,J,cov,mse] = nlinfit(X,Y,model1,params,opts);
yzero=[vals(11),C6(1),C7(1),C8(1),C10(1)];
y = ode15s(@(t,y)system6(t,y,vals), [0, X(end,1)], yzero);
yffit=deval(y,X(X(:,2)==1));
Yfit=[yffit(2,:)';yffit(3,:)';yffit(4,:)';yffit(5,:)'];
h1=figure(1);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,C6,'rd','LineWidth',5);
xlabel('time (days)')
ylabel('Hexane (mmole)')
set(gca,'TickDir','Out');
%axes('linewidth',5,'box','on','xtick',[],'ytick',[],'hittest','off');
axis([0 322 0 C6(1)])
h2=figure(2);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,C7,'rd','LineWidth',5);
xlabel('time (days)')
ylabel('Heptane (mmole)')
set(gca,'TickDir','Out');
axis([0 322 0 C7(1)])
h3=figure(3);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,C8,'rd','LineWidth',5);
xlabel('time (days)')
ylabel('Octane (mmole)')
set(gca,'TickDir','Out');
axis([0 322 0 C8(1)])
h4=figure(4);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',5);
plot(time,C10,'rd','LineWidth',5);
xlabel('time (days)')
ylabel('Decane (mmole)')
set(gca,'TickDir','Out');
axis([0 322 0 C10(1)])
figure(h1)
line(time,yffit(2,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h2)
line(time,yffit(3,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h3)
line(time,yffit(4,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h4)
line(time,yffit(5,:)','color','r','LineWidth',5)
% yy = ode15s(@system5, [0, X(end)], yzero);
cost_func = 'NMSE';

fit = goodnessOfFit(Yfit,Y,cost_func)
fit1=goodnessOfFit(yffit(2,:)',C6,cost_func)
fit2=goodnessOfFit(yffit(3,:)',C7,cost_func)
fit3=goodnessOfFit(yffit(4,:)',C8,cost_func)
fit4=goodnessOfFit(yffit(5,:)',C10,cost_func)
    function err = subfun(params,X)
        yzero2=[params(11),C6(1),C7(1),C8(1),C10(1)];
        yyy = ode15s(@(ttt,yyy)system5(ttt,yyy,params), [0,X(end,1)], yzero2);
        yint = deval(yyy,X(X(:,2)==1));% X(:,2)==1 refers to the positions where the same column is 1
        yfit=[yint(2,:), yint(3,:), yint(4,:), yint(5,:)]';
        %AA = (Y-yfit).^2;
        err=yfit;
    end
    function yfit1 = subfun1(vals,X)
        yzero1=[vals(11),C6(1),C7(1),C8(1),C10(1)];
        yy = ode15s(@(tt,yy)system6(tt,yy,vals), [0,X(end,1)], yzero1);
        yint1 = deval(yy,X(X(:,2)==1));%
        yfit1=[yint1(2,:), yint1(3,:), yint1(4,:), yint1(5,:)]';
    end

end

function dydt=system5(t,y,params)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=params(1);
Kf=params(2);
Kg_C6=params(3);
Kg_C7=params(4);
Kg_C8=params(5);
Kg_C10=params(6);
tau1=params(7);
tau2=params(8);
tau3=params(9);
tau4=params(10);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau1,  g_C6=0;  else
    g_C6=y(2)/(Kg_C6+y(2));end

if t<tau2,  g_C7=0;   else
    g_C7=y(3)/(Kg_C7+y(3));    end

if t<tau3,    g_C8=0;   else
    g_C8 =y(4)/(Kg_C8+y(4));    end

if t<tau4,   g_C10=0;   else
    g_C10=y(5)/(Kg_C10+y(5));    end


% ODE
dydt=[mu*B*(min(f, g_C6)+min(f, g_C7)+min(f, g_C8)+min(f,g_C10))-d*B;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_C7)+Cin;
    -1/r*mu*B*min(f,g_C8)+Cin;
    -1/r*mu*B*min(f,g_C10)+Cin];
end
function dydt1=system6(t,y,vals)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=vals(1);
Kf=vals(2);
Kg_C6=vals(3);
Kg_C7=vals(4);
Kg_C8=vals(5);
Kg_C10=vals(6);
tau1=vals(7);
tau2=vals(8);
tau3=vals(9);
tau4=vals(10);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau1,  g_C6=0;  else
    g_C6=y(2)/(Kg_C6+y(2));end

if t<tau2,  g_C7=0;   else
    g_C7=y(3)/(Kg_C7+y(3));    end

if t<tau3,    g_C8=0;   else
    g_C8 =y(4)/(Kg_C8+y(4));    end

if t<tau4,   g_C10=0;   else
    g_C10=y(5)/(Kg_C10+y(5));    end



% ODE
dydt1=[mu*B*(min(f, g_C6)+min(f, g_C7)+min(f, g_C8)+min(f,g_C10))-d*B;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_C7)+Cin;
    -1/r*mu*B*min(f,g_C8)+Cin;
    -1/r*mu*B*min(f,g_C10)+Cin];

end

%%
%%
function [vals,ci,params]=fit_to_pentane_data
%dbstop if caught error
% data
data = dataset('xlsfile', 'datatariq', 'Sheet',16);%data can be obtained from Mohamad Shahimin et al. 2016
time=data.Day;
T = [time];
C5=data.C5mmole;
Y=C5;
X = T;
params0=rand(1,2);
%params0=[154.949582600892,200];
model=@(params,X)subfun(params, X);% @(params) use to indicate taht the function depends only on params variable and others are fixed
%lb=[100,20];
lb=zeros(1,length(params0));
ub=inf*ones(1,length(params0));
%ub=[300,200];
[params,RSS]= lsqcurvefit(model,params0,X,Y,lb,ub);
model1=@(vals,X)subfun1(vals, X);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[vals,r,J,cov,mse] = nlinfit(X,Y,model1,params,opts);
ci = nlparci(vals,r,'covar',cov);
n=length(time);
k=length(vals);
%AIC_c=2*k+n*log(resnorm)
AIC_c=2*k+n*log(mse)
yzero=[0.0138,C5(1)];
y = ode15s(@(t,y)system6(t,y,vals), [0, X(end)], yzero);
yffit=deval(y,T);
Yfit=[yffit(2,:)'];
h1=figure(1);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,C5,'rd','linewidth',5);
xlabel('time (days)')
ylabel('Pentane (mmole)')
set(gca,'TickDir','Out');
axis([0 532 0 C5(1)])

cost_func = 'NMSE';
fit = goodnessOfFit(Yfit,Y,cost_func)

    function err = subfun(params,X)
        yzero1=[0.0138,C5(1)];
        yy = ode15s(@(tt,yy)system5(tt,yy,params), [0,X(end)], yzero1);
        yint = deval(yy,X);
        yfit=[yint(2,:)]';
        err=yfit;
    end
    function yfit1 = subfun1(vals,X)
        yzero2=[0.0138,C5(1)];
        yyy = ode15s(@(ttt,yyy)system6(ttt,yyy,vals), [0,X(end)], yzero2);
        yint1 = deval(yyy,X);
        yfit1=[yint1(2,:)]';
    end

end

function dydt=system5(t,y,params)
% initializations
Cin=0; %0 for lab data
mu=2;%1.3*7; %estimate (growth)
%mu=params(11);
d=0;%0.3*7; %estimate (death)
theta=0.2; %estimate (theta)
%theta=params(12);
r=0.4; %estimate (yield coefficient) - 0.1 (Hao), 0.4 (Julia)
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_C5=params(1);
Lag_C5=params(2);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<Lag_C5,  g_C5=0;  else
    g_C5=y(2)/(Kg_C5+y(2));end


% ODE
dydt(1)=mu*B*(min(f, g_C5))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_C5)+Cin;
dydt=[dydt(1),dydt(2)]';
end
function dydt=system6(t,y,vals)
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_C5=vals(1);
Lag_C5=vals(2);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<Lag_C5,  g_C5=0;  else
    g_C5=y(2)/(Kg_C5+y(2));end

% ODE
dydt(1)=mu*B*(min(f, g_C5))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_C5)+Cin;
dydt=[dydt(1),dydt(2)]';
end
%%
%% Code for 3-MC6, 2-MC7, 4-MC7, 2-MC8. Data is available upon request.
function [vals, ci, params]=fit_to_Isoalkanes_data
%dbstop if caught error
%% data
data = dataset('xlsfile', 'datatariq', 'Sheet',13);
time=data.Day;
T = [time; time; time;time];
MC36=data.MC36mmole;
MC27=data.MC27mmole;
MC47=data.MC47mmole;
MC28=data.MC28mmole;
Y=[MC36; MC27 ; MC47;MC28];
dsid = [ones(length(time),1); 2*ones(length(time),1); 3*ones(length(time),1);4*ones(length(time),1)];
X = [T dsid];
params0=rand(1,8);
%params0=[147.763781233381*MC36(1),823.797200113276*MC27(1),179.629199806177*MC47(1),413.434328041012*MC28(1),25,25,25,25];
model=@(params,X)subfun(params, X);% @(params) use to indicate taht the function depends only on params variable and others are fixed
lb=zeros(1,length(params0));
%lb=[0,0,0,0,0,0,0,0];
ub=1000000*ones(1,length(params0));
%ub=[1000000,1000000,1000000,10000000,430,847,665,670];
[params,RSS] = lsqcurvefit(model,params0,X,Y,lb,ub);
model1=@(vals,X)subfun1(vals, X);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[vals,r,J,cov,mse] = nlinfit(X,Y,model1,params,opts);
yzero=[0.000573104641468882,MC36(1),MC27(1),MC47(1),MC28(1)];
y = ode15s(@(t,y)system6(t,y,vals), [0, X(end,1)], yzero);
yffit=deval(y,X(X(:,2)==1));
Yfit=[yffit(2,:)';yffit(3,:)';yffit(4,:)';yffit(5,:)'];
h1=figure(1);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,MC36,'rd','linewidth',5);
xlabel('time (days)')
ylabel('3-MC_6(mmole)')
set(gca,'TickDir','Out');
axis([0 1700 0 MC36(1)])
h2=figure(2);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,MC27,'rd','linewidth' ,5);
xlabel('time (days)')
ylabel('2-MC_7 (mmole)')
set(gca,'TickDir','Out');
axis([0 1700 0 MC27(1)])
h3=figure(3);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,MC47,'rd','linewidth',5);
xlabel('time (days)')
ylabel('4-MC_7 (mmole)')
set(gca,'TickDir','Out');
axis([0 1700 0 MC47(1)])
h4=figure(4);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,MC28,'rd','linewidth',5);
xlabel('time (days)')
ylabel('2-MC_8 (mmole)')
set(gca,'TickDir','Out');
axis([0 1700 0 MC28(1)])
%axis([0 1601 0.11 0.15])
figure(h1)
line(time,yffit(2,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h2)
line(time,yffit(3,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h3)
line(time,yffit(4,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h4)
line(time,yffit(5,:)','color','r','LineWidth',5)
%legend('measured','fitted')
cost_func = 'NMSE';
fit = goodnessOfFit(Yfit,Y,cost_func)
fit1=goodnessOfFit(yffit(2,:)',MC36,cost_func)
fit2=goodnessOfFit(yffit(3,:)',MC27,cost_func)
fit3=goodnessOfFit(yffit(4,:)',MC47,cost_func)
fit4=goodnessOfFit(yffit(5,:)',MC28,cost_func)

    function err = subfun(params,X)
        yzero1=[0.000573104641468882,MC36(1),MC27(1),MC47(1),MC28(1)];
        yy = ode15s(@(tt,yy)system5(tt,yy,params), [0,X(end,1)], yzero1);
        yint = deval(yy,X(X(:,2)==1));% X(:,2)==1 refers to the positions where the same column is 1
        yfit=[yint(2,:), yint(3,:), yint(4,:),yint(5,:)]';
        %AA = (Y-yfit).^2;
        err=yfit;
    end
    function yfit1 = subfun1(vals,X)
        yzero2=[0.000573104641468882,MC36(1),MC27(1),MC47(1),MC28(1)];
        yyy = ode15s(@(ttt,yyy)system6(ttt,yyy,vals), [0,X(end,1)], yzero2);
        yint1 = deval(yyy,X(X(:,2)==1));% X(:,2)==1 refers to the positions where the same column is 1
        yfit1=[yint1(2,:), yint1(3,:), yint1(4,:),yint1(5,:)]';
        %AA = (Y-yfit).^2;
        %err=yfit;
    end

end
function dydt=system5(t,y,params)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_MC36=params(1);
Kg_MC27=params(2);
Kg_MC47=params(3);
Kg_MC28=params(4);
tau_MC36=params(5);
tau_MC27=params(6);
tau_MC47=params(7);
tau_MC28=params(8);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau_MC36, g_MC36=0;  else
    g_MC36=y(2)/(Kg_MC36+y(2));end


if t<tau_MC27,    g_MC27=0;   else
    g_MC27=y(3)/(Kg_MC27+y(3));    end


if t<tau_MC47,   g_MC47=0;   else
    g_MC47 =y(4)/(Kg_MC47+y(4));    end
if t<tau_MC28,    g_MC28=0;   else
    g_MC28=y(5)/(Kg_MC28+y(5));    end



% ODE
dydt(1)=mu*B*(min(f, g_MC36)+min(f, g_MC27)+min(f, g_MC47)+min(f, g_MC28))-d*B;
dydt(2)=-1/r*mu*B*min(f, g_MC36)+Cin;
dydt(3)=-1/r*mu*B*min(f, g_MC27)+Cin;
dydt(4)=-1/r*mu*B*min(f, g_MC47)+Cin;
dydt(5)=-1/r*mu*B*min(f, g_MC28)+Cin;
dydt=[dydt(1),dydt(2),dydt(3),dydt(4),dydt(5)]';
end
function dydt=system6(t,y,vals)
% initializations
Cin=0;
mu=2;

d=0;
theta=0.2;

r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_MC36=vals(1);
Kg_MC27=vals(2);
Kg_MC47=vals(3);
Kg_MC28=vals(4);
tau_MC36=vals(5);
tau_MC27=vals(6);
tau_MC47=vals(7);
tau_MC28=vals(8);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau_MC36, g_MC36=0;  else
    g_MC36=y(2)/(Kg_MC36+y(2));end


if t<tau_MC27,    g_MC27=0;   else
    g_MC27=y(3)/(Kg_MC27+y(3));    end


if t<tau_MC47,   g_MC47=0;   else
    g_MC47 =y(4)/(Kg_MC47+y(4));    end

if t<tau_MC28,    g_MC28=0;   else
    g_MC28=y(5)/(Kg_MC28+y(5));    end




% ODE
dydt(1)=mu*B*(min(f, g_MC36)+min(f, g_MC27)+min(f, g_MC47)+min(f, g_MC28))-d*B;
dydt(2)=-1/r*mu*B*min(f, g_MC36)+Cin;
dydt(3)=-1/r*mu*B*min(f, g_MC27)+Cin;
dydt(4)=-1/r*mu*B*min(f, g_MC47)+Cin;
dydt(5)=-1/r*mu*B*min(f, g_MC28)+Cin;
dydt=[dydt(1),dydt(2),dydt(3),dydt(4),dydt(5)]';
end
%%
%% code for 2-MC5
function [vals,ci,params]=fit_to_2MC5_data
%dbstop if caught error
close all
clear all
clc
format long
T=[0;600;1300;1600];
MC25=[0.357391506149919;0.207287073566953;0.017869575307496;0];%data from Siddique et al. 2015

Y=MC25;
X = T;

params0=rand(1,2);
%params0=[165.940127757038,0.00856575823994166];
%params0=[166,];
model=@(params,X)subfun(params, X);
%lb=[0,0];
lb=zeros(1,length(params0));
ub=inf*ones(1,length(params0));
%ub=[10000,100];
[params,RSS]= lsqcurvefit(model,params0,X,Y,lb,ub);
model1=@(vals,X)subfun1(vals, X);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[vals,r,J,cov,mse] = nlinfit(X,Y,model1,params,opts);

yzero=[vals(2),MC25(1)];
y = ode15s(@(t,y)system6(t,y,vals), [0, X(end)], yzero);
yffit=deval(y,T);
Yfit=[yffit(2,:)'];
h1=figure(1);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(X,MC25,'rd','linewidth',5);
xlabel('time (days)')
ylabel('2-MC_5 (mmole)')
set(gca,'TickDir','Out');
axis([0 1600 0 MC25(1)])
figure(h1)
line(T,yffit(2,:)','color','r','LineWidth',5)
%legend('measured','fitted')
cost_func = 'NMSE';
fit = goodnessOfFit(Yfit,Y,cost_func)

    function err = subfun(params,X)
        yzero1=[params(2),MC25(1)];
        yy = ode15s(@(tt,yy)system5(tt,yy,params), [0,X(end)], yzero1);
        yint = deval(yy,X);
        yfit=[yint(2,:)]';
        err=yfit;
    end
    function yfit1 = subfun1(vals,X)
        yzero2=[vals(2),MC25(1)];
        yyy = ode15s(@(ttt,yyy)system6(ttt,yyy,vals), [0,X(end)], yzero2);
        yint1 = deval(yyy,X);
        yfit1=[yint1(2,:)]';
        
    end

end

function dydt=system5(t,y,params)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_MC25=params(1);
Lag_MC25=23;
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<Lag_MC25,  g_MC25=0;  else
    g_MC25=y(2)/(Kg_MC25+y(2));end



% ODE
dydt(1)=mu*B*(min(f, g_MC25))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_MC25)+Cin;
dydt=[dydt(1),dydt(2)]';
end
function dydt=system6(t,y,vals)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_MC25=vals(1);
Lag_MC25=23;
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<Lag_MC25,  g_MC25=0;  else
    g_MC25=y(2)/(Kg_MC25+y(2));end



% ODE
dydt(1)=mu*B*(min(f, g_MC25))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_MC25)+Cin;
dydt=[dydt(1),dydt(2)]';
end
%%
%%fit to BTEX hydrocarbons
function [vals,params,ci]=fit_to_BTEX_data
%dbstop if caught error
close all
clear all
clc
format long
% data
data = dataset('xlsfile', 'datatariq', 'Sheet',11);%data is from siddique et al. 2007
time=data.Day;
T = [time; time; time];
Toluene=data.Toluenemmole;
mXylenes=data.mXylenemmole;
oXylene=data.oXylenemmole;
Y=[Toluene;mXylenes ; oXylene];
dsid = [ones(length(time),1); 2*ones(length(time),1); 3*ones(length(time),1)];
X = [T dsid];
params0=rand(1,3);
%params0=[4.47784190084274,85.0815935026002,17.5097400254338,30,70,60.0000000000000];
model=@(params,X)subfun(params, X);% @(params) use to indicate taht the function depends only on params variable and others are fixed
lb=zeros(1,length(params0));
%lb=[0,0,0,20,60,17];
%ub=[1000000,1000000,1000000,50,150,60];
ub=inf*ones(1,length(params0));
params = lsqcurvefit(model,params0,X,Y,lb,ub);
model1=@(vals,X)subfun1(vals, X);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[vals,r,J,cov,mse] = nlinfit(X,Y,model1,params,opts);
yzero=[0.000573104641468882,Toluene(1),mXylenes(1),oXylene(1)];
y = ode15s(@(t,y)system6(t,y,vals), [0, X(end,1)], yzero);
yffit=deval(y,X(X(:,2)==1));
Yfit=[yffit(2,:)';yffit(3,:)';yffit(4,:)'];
h1=figure(1);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,Toluene,'rd','linewidth',5);
xlabel('time (days)')
ylabel('Toluene (mmole)')
set(gca,'TickDir','Out');
axis([0 252 0 Toluene(1)])
h2=figure(2);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,mXylenes,'rd','linewidth' ,5);
xlabel('time (days)')
ylabel('m-,p-Xylenes (mmole)')
set(gca,'TickDir','Out');
axis([0 252 0 mXylenes(1)])
h3=figure(3);
set(0,'DefaultAxesFontSize',32,'DefaultTextFontSize',32,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(time,oXylene,'rd','linewidth',5);
xlabel('time (days)')
ylabel('o-Xylene (mmole)')
set(gca,'TickDir','Out');
axis([0 252 0 oXylene(1)])
%axis([0 1601 0.11 0.15])
figure(h1)
line(time,yffit(2,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h2)
line(time,yffit(3,:)','color','r','LineWidth',5)
%legend('measured','fitted')
figure(h3)
line(time,yffit(4,:)','color','r','LineWidth',5)
%legend('measured','fitted')
cost_func = 'NMSE';
fit = goodnessOfFit(Yfit,Y,cost_func)
fit1=goodnessOfFit(yffit(2,:)',Toluene,cost_func)
fit2=goodnessOfFit(yffit(3,:)',mXylenes,cost_func)
fit3=goodnessOfFit(yffit(4,:)',oXylene,cost_func)

    function err = subfun(params,X)
        yzero1=[0.000573104641468882,Toluene(1),mXylenes(1),oXylene(1)];
        yyy = ode15s(@(ttt,yyy)system5(ttt,yyy,params), [0,X(end,1)], yzero1);
        yint = deval(yyy,X(X(:,2)==1));% X(:,2)==1 refers to the positions where the same column is 1
        yfit=[yint(2,:), yint(3,:), yint(4,:)]';
        %AA = (Y-yfit).^2;
        err=yfit;
    end
    function yfit1 = subfun1(vals,X)
        yzero2=[0.000573104641468882,Toluene(1),mXylenes(1),oXylene(1)];
        yy = ode15s(@(tt,yy)system6(tt,yy,vals), [0,X(end,1)], yzero2);
        yint1 = deval(yy,X(X(:,2)==1));%
        yfit1=[yint1(2,:), yint1(3,:), yint1(4,:)]';
    end

end
function dydt=system5(t,y,params)
% initializations
Cin=0;
mu=2;

d=0;
theta=0.2;

r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_Toluene=params(1);
Kg_mXylenes=params(2);
Kg_oXylene=params(3);
tau_Toluene=params(4);
tau_mXylenes=params(5);
tau_oXylene=params(6);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau_Toluene, g_Toluene=0;  else
    g_Toluene=y(2)/(Kg_Toluene+y(2));end

if t<tau_mXylenes,    g_mXylenes=0;   else
    g_mXylenes =y(3)/(Kg_mXylenes+y(3));    end

if t<tau_oXylene,   g_oXylene=0;   else
    g_oXylene =y(4)/(Kg_oXylene+y(4));    end


% ODE
dydt(1)=mu*B*(min(f, g_Toluene)+min(f, g_mXylenes)+min(f, g_oXylene))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_Toluene)+Cin;
dydt(3)=-1/r*mu*B*min(f,g_mXylenes)+Cin;
dydt(4)=-1/r*mu*B*min(f,g_oXylene)+Cin;
dydt=[dydt(1),dydt(2),dydt(3),dydt(4)]';
end
function dydt=system6(t,y,vals)
% initializations
Cin=0;
mu=2;

d=0;
theta=0.2;

r=0.4;

Tn=327.615787410708;
Kf=0.262218847536511;
Kg_Toluene=vals(1);
Kg_mXylenes=vals(2);
Kg_oXylene=vals(3);
tau_Toluene=vals(4);
tau_mXylenes=vals(5);
tau_oXylene=vals(6);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);

% Lag times
if t<tau_Toluene, g_Toluene=0;  else
    g_Toluene=y(2)/(Kg_Toluene+y(2));end

if t<tau_mXylenes,    g_mXylenes=0;   else
    g_mXylenes =y(3)/(Kg_mXylenes+y(3));    end

if t<tau_oXylene,   g_oXylene=0;   else
    g_oXylene =y(4)/(Kg_oXylene+y(4));    end


% ODE
dydt(1)=mu*B*(min(f, g_Toluene)+min(f, g_mXylenes)+min(f, g_oXylene))-d*B;
dydt(2)=-1/r*mu*B*min(f,g_Toluene)+Cin;
dydt(3)=-1/r*mu*B*min(f,g_mXylenes)+Cin;
dydt(4)=-1/r*mu*B*min(f,g_oXylene)+Cin;
dydt=[dydt(1),dydt(2),dydt(3),dydt(4)]';
end
%%
%% methane model validation with methane from parafinic diluent in CNUL MFT (see figure 1B or the manuscript). Data obatiend from Mohamad Shhimin and Siddique 2017a
function methane=validation_comparison_with_CNUL_Methane
%dbstop if caught error

% data
data = dataset('xlsfile', 'datatariq', 'Sheet',22);
X=data.Day;
%CNRLtime=data.Day1;
Y=data.Albianmethane;

C50=0.4269;
C60=0.1126;
MC25=0.2448;
yzero=[0.00329283300847435,C50,C60,MC25,0,0,0];
y1zero=[C50,C60,MC25,0,0,0];
y2zero=[C50,C60,MC25,0,0,0];

yy = ode15s(@system5, [0, X(end)], yzero);
yfit=deval(yy,X);
yy1 = ode15s(@system51, [0, X(end)], y1zero);
yfit1=deval(yy1,X);
yy2 = ode15s(@system52, [0, X(end)], y2zero);
yfit2=deval(yy2,X);
methan=0.65*(yfit(5,:)*4+yfit(6,:)*4.75+yfit(7,:)*4.75);
methan1=0.65*(yfit1(4,:)*4+yfit1(5,:)*4.75+yfit1(6,:)*4.75);
methan2=0.65*(yfit2(4,:)*4+yfit2(5,:)*4.75+yfit2(6,:)*4.75);
cost_func = 'NMSE';

fit = goodnessOfFit(methan',Y,cost_func)
fit1 = goodnessOfFit(methan1',Y,cost_func)
fit2 = goodnessOfFit(methan2',Y,cost_func)

h1=figure(1);
set(0,'DefaultAxesFontSize',28,'DefaultTextFontSize',28,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(X,Y,'rd','linewidth',5);
xlabel('time (days)')
ylabel('CH_4 (mmole)')
set(gca,'TickDir','Out');
axis([0 X(end) 0 max(Y)])
hold on;

plot(X,methan,'k',X,methan1,'--b',X,methan2,':k','LineWidth',5)

end

function dydt=system5(t,y)
% initializations
Cin=0;
mu=2;
d=0;
theta=0.2;
r=0.31;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_C6=365.0518;%430.327575377768;%137.749494059463;
Kg_MC25=166;%165.9401;%166;%130.2141-201.6662
Kg_C5=154.949582600892;%[96.4323835305526];%154.949582600892;
Lag_C6=26;%26;%115
Lag_C5=200;
Lag_MC25=23;
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);
% Lag times
if t<Lag_C5,  g_C5=0;  else
    g_C5=y(2)/(Kg_C5+y(2));end

if t<Lag_C6,  g_C6=0;   else
    g_C6=y(3)/(Kg_C6+y(3));    end
if t<Lag_MC25,    g_MC25=0;   else
    g_MC25=y(4)/(Kg_MC25+y(4));    end


% ODE
dydt=[mu*B*(min(f, g_C5)+min(f, g_C6)+min(f,g_MC25))-d*B;
    -1/r*mu*B*min(f,g_C5)+Cin;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_MC25)+Cin;
    1/r*mu*B*min(f,g_C5);
    1/r*mu*B*min(f,g_C6);
    1/r*mu*B*min(f,g_MC25)];
end
function dydt1=system51(t1,y1)
% initializations

lagC15=294;
lagC16=9*7;
lag2MC51=600;



% Lag times
if t1<lagC15,  Kg_C15=0;   else
    Kg_C15=0.0008576;    end

if t1<lagC16,  Kg_C16=0;  else
    Kg_C16=(0.27/10)/7;end

if t1<lag2MC51,  Kg_2MC51=0;   else
    Kg_2MC51=0.002281;    end



% ODE
dydt1=[-Kg_C15;
    -Kg_C16;
    -Kg_2MC51;
    Kg_C15;
    Kg_C16;
    Kg_2MC51];

end
function dydt2=system52(t2,y2)
% initializations


lagC25=294;
lagC26=15*7;
lag2MC52=600;




% Lag times
if t2<lagC25,  Kg_C25=0;  else
    Kg_C25=0.01117;end
if t2<lagC26,  Kg_C26=0;  else
    Kg_C26=0.09/7;end

if t2<lag2MC52,  Kg_2MC52=0;   else
    Kg_2MC52=0.003501;    end



% ODE
dydt2=[-Kg_C25*y2(1);
    -Kg_C26*y2(2);
    -Kg_2MC52*y2(3);
    Kg_C25*y2(1);
    Kg_C26*y2(2);
    Kg_2MC52*y2(3)];

end
%%
%% validation and comparison using CNRL MFT methane
function methane=validation_comparison_with_CNRL_Methane
%dbstop if caught error
close all
clear all
clc
format long
% data
data = dataset('xlsfile', 'datatariq', 'Sheet',21);
Albiantime=data.Day;
CNRLtime=data.Day1;
Y=data.CNRLmethane(1:end-1);
T=[CNRLtime(1:end-1)];
X = T;

C60=0.105592944998840;

C70=0.151681468915278;
C80=0.091919810907818;

C100=0.014758591608686;

C90=0.026508654295961;

MC250=0.02088651659317;

MC360=0.071849116854605;

MC270=0.036767924363127;

MC470=0.030639936969273;

MC280=0.015593326056448;


yzero=[0.000573104641468882,C60,C70,C80,C90,C100,MC250,MC360,MC270,MC470,MC280,0,0,0,0,0,0,0,0,0,0];
y1zero=[C60,C70,C80,C90,C100,MC250,MC360,MC270,MC470,MC280,0,0,0,0,0,0,0,0,0,0];
y2zero=[C60,C70,C80,C90,C100,MC250,MC360,MC270,MC470,MC280,0,0,0,0,0,0,0,0,0,0];

yy = ode15s(@system5, [0, X(end)], yzero);
yfit=deval(yy,X);
yy1 = ode15s(@system51, [0, X(end)], y1zero);
yfit1=deval(yy1,X);
yy2 = ode15s(@system52, [0, X(end)], y2zero);
yfit2=deval(yy2,X);

methan=.7*(yfit(12,:)*4.75+yfit(13,:)*5.5+yfit(14,:)*6.25+yfit(15,:)*7+yfit(16,:)*7.75+...
    yfit(17,:)*4.75+yfit(18,:)*5.5+yfit(19,:)*6.25+yfit(20,:)*6.25+yfit(21,:)*7);
methan1=.7*(yfit1(11,:)*4.75+yfit1(12,:)*5.5+yfit1(13,:)*6.25+yfit1(14,:)*7+yfit1(15,:)*7.75+...
    yfit1(16,:)*4.75+yfit1(17,:)*5.5+yfit1(18,:)*6.25+yfit1(19,:)*6.25+yfit1(20,:)*7);
methan2=.7*(yfit2(11,:)*4.75+yfit2(12,:)*5.5+yfit2(13,:)*6.25+yfit2(14,:)*7+yfit2(15,:)*7.75+...
    yfit2(16,:)*4.75+yfit2(17,:)*5.5+yfit2(18,:)*6.25+yfit2(19,:)*6.25+yfit2(20,:)*7);

cost_func = 'NMSE';

ffit = goodnessOfFit(methan',Y,cost_func)
fit1 = goodnessOfFit(methan1',Y,cost_func)
fit2 = goodnessOfFit(methan2',Y,cost_func)


h1=figure(1);
set(0,'DefaultAxesFontSize',28,'DefaultTextFontSize',28,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(X,Y,'rd','linewidth',5);
xlabel('time (days)')
ylabel('CH_4 (mmole)')
set(gca,'TickDir','Out');
axis([0 X(end) 0 max(methan)])
hold on;

plot(X,methan,'k',X,methan1,'--b',X,methan2,':k','LineWidth',5)


end

function dydt=system5(t,y)
% initializations
Cin=0;
mu=2.5;

d=0;
theta=0.2;

r=0.25;

Tn=327.615787410708;

Kf=0.262218847536511;

Kg_C6=430.327575377768;

Kg_C7=269.144904586922;

Kg_C8= 90.5884032453207;

Kg_C10= 12.4460621918727;
Kg_MC25=165.9401;
Kg_MC36=1.446139728984238*1e2;
Kg_MC27=3.204133584773790*1e2;
Kg_MC47=1.703470907219932*1e2;
Kg_MC28=3.358911095487180*1e3;
Kg_C9=0.818706893477188;
Lag_C9=70;
Lag_C6=26;
Lag_C7=60;
Lag_C8=70;
Lag_C10=5;
Lag_MC25=23;
Lag_MC36=0.249999913563191*1e2;
Lag_MC27=0.250000143067640*1e2;
Lag_MC47=0.249999968008542*1e2;
Lag_MC28=0.250000105120839*1e2;
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);

% Lag times
if t<Lag_C6,  g_C6=0;  else
    g_C6=y(2)/(Kg_C6+y(2));end

if t<Lag_C7,  g_C7=0;   else
    g_C7=y(3)/(Kg_C7+y(3));    end

if t<Lag_C8,    g_C8=0;   else
    g_C8 =y(4)/(Kg_C8+y(4));    end
if t<Lag_C9,  g_C9=0;  else
    g_C9=y(5)/(Kg_C9+y(5));end

if t<Lag_C10,   g_C10=0;   else
    g_C10=y(6)/(Kg_C10+y(6));    end
if t<Lag_MC25,  g_MC25=0;  else
    g_MC25=y(7)/(Kg_MC25+y(7));end
if t<Lag_MC36, g_MC36=0;  else
    g_MC36=y(8)/(Kg_MC36+y(8));end

if t<Lag_MC27,    g_MC27=0;   else
    g_MC27=y(9)/(Kg_MC27+y(9));    end

if t<Lag_MC47,   g_MC47=0;   else
    g_MC47 =y(10)/(Kg_MC47+y(10));    end

if t<Lag_MC28,    g_MC28=0;   else
    g_MC28=y(11)/(Kg_MC28+y(11));    end




% ODE
dydt=[mu*B*(min(f, g_C6)+min(f, g_C7)+min(f, g_C8)+min(f, g_C9)+min(f,g_C10)+min(f,g_MC25)+min(f, g_MC36)+min(f, g_MC27)+min(f, g_MC47)+min(f, g_MC28))-d*B;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_C7)+Cin;
    -1/r*mu*B*min(f,g_C8)+Cin;
    -1/r*mu*B*min(f,g_C9)+Cin;
    -1/r*mu*B*min(f,g_C10)+Cin;
    -1/r*mu*B*min(f,g_MC25)+Cin;
    -1/r*mu*B*min(f,g_MC36)+Cin;
    -1/r*mu*B*min(f,g_MC27)+Cin;
    -1/r*mu*B*min(f,g_MC47)+Cin;
    -1/r*mu*B*min(f,g_MC28)+Cin;
    1/r*mu*B*min(f,g_C6);
    1/r*mu*B*min(f,g_C7);
    1/r*mu*B*min(f,g_C8);
    1/r*mu*B*min(f,g_C9);
    1/r*mu*B*min(f,g_C10);
    1/r*mu*B*min(f,g_MC25);
    1/r*mu*B*min(f,g_MC36);
    1/r*mu*B*min(f,g_MC27);
    1/r*mu*B*min(f,g_MC47);
    1/r*mu*B*min(f,g_MC28)];
end
function dydt1=system51(t1,y1)
% initializations

Lag_C91=77;
Lag_C61=9*7;
Lag_C71=10*7;
Lag_C81=10*7;
Lag_C101=2*7;
Lag_MC251=600;
Lag_MC361=455;
Lag_MC271=845;
Lag_MC471=665;
Lag_MC281=665;
if t1<Lag_C61,  g_C61=0;  else
    g_C61=(0.27/10)/7;end

if t1<Lag_C71,  g_C71=0;   else
    g_C71=(0.47/10)/7;    end

if t1<Lag_C81,    g_C81=0;   else
    g_C81 =(0.74/10)/7;    end
if t1<Lag_C91,  g_C91=0;  else
    g_C91=2.664e-05;end

if t1<Lag_C101,   g_C101=0;   else
    g_C101=(0.41/10)/7;    end
if t1<Lag_MC251,  g_MC251=0;  else
    g_MC251=0.0002281;end
if t1<Lag_MC361, g_MC361=0;  else
    g_MC361=0.0001816;end

%
if t1<Lag_MC271,    g_MC271=0;   else
    g_MC271=0.00023;    end

%
if t1<Lag_MC471,   g_MC471=0;   else
    g_MC471 =0.0001936;    end

if t1<Lag_MC281,    g_MC281=0;   else
    g_MC281=0.0001772;    end



% ODE
dydt1=[-g_C61;
    -g_C71;
    -g_C81;
    -g_C91;
    -g_C101;
    -g_MC251;
    -g_MC361;
    -g_MC271;
    -g_MC471;
    -g_MC281;
    g_C61;
    g_C71;
    g_C81;
    g_C91;
    g_C101;
    g_MC251;
    g_MC361;
    g_MC271;
    g_MC471;
    g_MC281];

end
function dydt2=system52(t2,y2)
%% initializations
Lag_C92=77;
Lag_C62=15*7;
Lag_C72=17*7;
Lag_C82=12*7;
Lag_C102=5*7;
Lag_MC252=600;
Lag_MC362=455;
Lag_MC272=845;
Lag_MC472=665;
Lag_MC282=665;

if t2<Lag_C62,  g_C62=0;  else
    g_C62=0.09/7;end

if t2<Lag_C72,  g_C72=0;   else
    g_C72=0.19/7;    end

if t2<Lag_C82,    g_C82=0;   else
    g_C82 =0.34/7;    end
if t2<Lag_C92,  g_C92=0;  else
    g_C92=0.01276;end

if t2<Lag_C102,   g_C102=0;   else
    g_C102=0.22/7;    end
if t2<Lag_MC252,  g_MC252=0;  else
    g_MC252=0.003501;end
if t2<Lag_MC362, g_MC362=0;  else
    g_MC362=0.003849;end


if t2<Lag_MC272,    g_MC272=0;   else
    g_MC272=0.005258;    end

if t2<Lag_MC472,   g_MC472=0;   else
    g_MC472 =0.005663;    end

if t2<Lag_MC282,    g_MC282=0;   else
    g_MC282=0.0006584;    end



% ODE
dydt2=[-g_C62*y2(1);
    -g_C72*y2(2);
    -g_C82*y2(3);
    -g_C92*y2(4);
    -g_C102*y2(5);
    -g_MC252*y2(6);
    -g_MC362*y2(7);
    -g_MC272*y2(8);
    -g_MC472*y2(9);
    -g_MC282*y2(10);
    g_C62*y2(1);
    g_C72*y2(2);
    g_C82*y2(3);
    g_C92*y2(4);
    g_C102*y2(5);
    g_MC252*y2(6);
    g_MC362*y2(7);
    g_MC272*y2(8);
    g_MC472*y2(9);
    g_MC282*y2(10)];


end
%%
%% Validation and comparison to other models using the data on Table S1.
function methane=Validation_with_data_TabeS1
%dbstop if caug

% data
data = dataset('xlsfile', 'datatariq', 'Sheet',17);
time=data.daykathy(1:32);
M=data.kathyfinal(1:32);
T = [time];
Y=M;

X = T;

yzero=[0.01,0.005801810164771,0.033928749625786,0.040269631445330, ...
    0.049924028652051,0.032969103240392,0.013187641296157,0,0,0,0,0,0];
y1zero=[0.005801810164771,0.033928749625786,0.040269631445330, ...
    0.049924028652051,0.032969103240392,0.013187641296157,0,0,0,0,0,0];
y2zero=[0.005801810164771,0.033928749625786,0.040269631445330, ...
    0.049924028652051,0.032969103240392,0.013187641296157,0,0,0,0,0,0];

y = ode15s(@system5, [0, X(end)], yzero);

yfit=deval(y,X);
y1 = ode15s(@system51, [0, X(end)], y1zero);

yfit1=deval(y1,X);
y2= ode15s(@system52, [0, X(end)], y2zero);

yfit2=deval(y2,X);
methan=(yfit(8,:)*4.75+yfit(9,:)*5.5+yfit(10,:)*6.25+yfit(11,:)*4.50+yfit(12,:)*5.25+yfit(13,:)*5.25)*0.6;
methan1=(yfit1(7,:)*4.75+yfit1(8,:)*5.5+yfit1(9,:)*6.25+yfit1(10,:)*4.50+yfit1(11,:)*5.25+yfit1(12,:)*5.25)*0.6;
methan2=(yfit2(7,:)*4.75+yfit2(8,:)*5.5+yfit2(9,:)*6.25+yfit2(10,:)*4.50+yfit2(11,:)*5.25+yfit2(12,:)*5.25)*0.6;
cost_func =  'NRMSE';
fitourmodel = goodnessOfFit(Y,methan',cost_func)
fitzeroorder = goodnessOfFit(Y,methan1',cost_func)
fitfirstorder = goodnessOfFit(Y,methan2',cost_func)
h1=figure(1);
set(0,'DefaultAxesFontSize',28,'DefaultTextFontSize',28,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(X,Y,'rd','linewidth',5);
xlabel('time (days)')
ylabel('CH_4 (mmol)')
set(gca,'TickDir','Out');
axis([0 X(end) 0 max(methan)])
hold on;

plot(X,methan,'k',X,methan1,'--b',X,methan2,':k','LineWidth',5)


end

function dydt=system5(t,y)
% initializations
Cin=0;
mu=3;

d=0;
theta=0.25;

r=0.2;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_Toluene=4.802;

Kg_mXylenes=93.2133;

Kg_oXylene=20.8375;
Kg_C6=100.06;

Kg_C7=238.94;

Kg_C8= 73.2553;
tau1=20.0000000000000;
tau2=60;
tau3=70;
tau_Toluene=30;
tau_mXylenes=70;
tau_oXylene=60;
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);


% Lag times
if t<tau1,  g_C6=0;  else
    g_C6=y(2)/(Kg_C6+y(2));end

if t<tau2,  g_C7=0;   else
    g_C7=y(3)/(Kg_C7+y(3));    end

if t<tau3,    g_C8=0;   else
    g_C8 =y(4)/(Kg_C8+y(4));    end
if t<tau_Toluene, g_Toluene=0;  else
    g_Toluene=y(2)/(Kg_Toluene+y(2));end

if t<tau_mXylenes,    g_mXylenes=0;   else
    g_mXylenes =y(3)/(Kg_mXylenes+y(3));    end

if t<tau_oXylene,   g_oXylene=0;   else
    g_oXylene =y(4)/(Kg_oXylene+y(4));    end


dydt=[mu*B*(min(f, g_C6)+min(f, g_C7)+min(f, g_C8)+min(f, g_Toluene)+min(f, g_mXylenes)+min(f, g_oXylene))-d*B;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_C7)+Cin;
    -1/r*mu*B*min(f,g_C8)+Cin;
    -1/r*mu*B*min(f,g_Toluene)+Cin;
    -1/r*mu*B*min(f,g_mXylenes)+Cin;
    -1/r*mu*B*min(f,g_oXylene)+Cin;
    1/r*mu*B*min(f,g_C6);
    1/r*mu*B*min(f,g_C7);
    1/r*mu*B*min(f,g_C8);
    1/r*mu*B*min(f,g_Toluene);
    1/r*mu*B*min(f,g_mXylenes);
    1/r*mu*B*min(f,g_oXylene)];
end

function dydt1=system51(t1,y1)
% initializations
tau111=15*7;%9*7;
tau121=10*7;
tau131=10*7;
tau11=8*7;
tau12=16*7;
tau13=14*7;




% Lag times
if t1<tau111,  Kg_C16=0;  else
    Kg_C16=(0.27/10)/7;end

if t1<tau121,  Kg_C17=0;   else
    Kg_C17=(0.47/10)/7;    end

if t1<tau131,    Kg_C18=0;   else
    Kg_C18= (0.74/10)/7;    end
if t1<tau11,  Kg_Toluene=0;  else
    Kg_Toluene=(0.07/10)/7;end

if t1<tau12,  Kg_mXylenes=0;   else
    Kg_mXylenes=(0.07/10)/7;    end

if t1<tau13,    Kg_oXylene=0;   else
    Kg_oXylene= (0.09/10)/7;    end




% ODE
dydt1=[-Kg_C16;
    -Kg_C17;
    -Kg_C18;
    -Kg_Toluene;
    -Kg_mXylenes;
    -Kg_oXylene;
    Kg_C16;
    Kg_C17;
    Kg_C18;
    Kg_Toluene;
    Kg_mXylenes;
    Kg_oXylene];

end
function dydt2=system52(t2,y2)
% initializations
tau211=15*7;
tau221=17*7;
tau231=12*7;
tau21=7*7;
tau22=17*7;
tau23=14*7;



% Lag times
if t2<tau211,  Kg_C26=0;  else
    Kg_C26=0.09/7;end

if t2<tau221,  Kg_C27=0;   else
    Kg_C27=0.19/7;    end

if t2<tau231,    Kg_C28=0;   else
    Kg_C28= 0.34/7;    end
if t2<tau21,  Kg_Toluene=0;  else
    Kg_Toluene=0.09/7;end

if t2<tau22,  Kg_mXylenes=0;   else
    Kg_mXylenes=0.06/7;    end

if t2<tau23,    Kg_oXylene=0;   else
    Kg_oXylene= 0.15/7;    end




% ODE
dydt2=[-Kg_C26*y2(1);
    -Kg_C27*y2(2);
    -Kg_C28*y2(3);
    -Kg_Toluene*y2(4);
    -Kg_mXylenes*y2(5);
    -Kg_oXylene*y2(6);
    Kg_C26*y2(1);
    Kg_C27*y2(2);
    Kg_C28*y2(3);
    Kg_Toluene*y2(4);
    Kg_mXylenes*y2(5);
    Kg_oXylene*y2(6)];

end

%%
%%Conparison of our model predcitions with the predictiosn of otehr models
%%in the literature
function methane=comparison_with_syncrude_alkane
%dbstop if caught error
% data
data = dataset('xlsfile', 'datatariq', 'Sheet',17);
time=data.Daymethane(1:end-21);
M=data.Alkanemethane(1:end-21);
T = [time];
Y=M;

X = T;
yzero=[0.000509607249615741,0.782092444,0.918325884,0.861381872,0.599130025,0,0,0,0];
y1zero=[0.782092444,0.918325884,0.861381872,0.599130025,0,0,0,0];
y2zero=[0.782092444,0.918325884,0.861381872,0.599130025,0,0,0,0];

y = ode15s(@system5, [0, X(end)], yzero);

yfit=deval(y,X);
y1 = ode15s(@system51, [0, X(end)], y1zero);

yfit1=deval(y1,X);
y2= ode15s(@system52, [0, X(end)], y2zero);

yfit2=deval(y2,X);
methan=(yfit(6,:)*4.75+yfit(7,:)*5.5+yfit(8,:)*6.25+yfit(9,:)*7.75)*0.8;
methan1=(yfit1(5,:)*4.75+yfit1(6,:)*5.5+yfit1(7,:)*6.25+yfit1(8,:)*7.75)*0.8;
methan2=(yfit2(5,:)*4.75+yfit2(6,:)*5.5+yfit2(7,:)*6.25+yfit2(8,:)*7.75)*0.8;
cost_func = 'NMSE';

fitourmodel = goodnessOfFit(methan',Y,cost_func)
fitzeroorder = goodnessOfFit(methan1',Y,cost_func)
fitfirstorder = goodnessOfFit(methan2',Y,cost_func)
h1=figure(1);
set(0,'DefaultAxesFontSize',28,'DefaultTextFontSize',28,'DefaultAxesLineWidth',5,'DefaultLineMarkerSize', 18);
set(gca,'FontSize',32,'LineWidth',32);
plot(X,Y,'rd','linewidth',5);
xlabel('time (days)')
ylabel('CH_4 (mmole)')
set(gca,'TickDir','Out');
hold on
axis([0 175 0 16])

plot(X,methan,'k',X,methan1,'--b',X,methan2,':g','LineWidth',5)


end

function dydt=system5(t,y)
% initializations
Cin=0;
mu=1.95;
d=0;
theta=0.2;
r=0.4;
Tn=327.615787410708;
Kf=0.262218847536511;
Kg_C6=365.0518;
Kg_C7=269.144904586922;
Kg_C8= 107.9215;
Kg_C10= 12.4460621918727;
tau1=23.0000000000000;
tau2=60;
tau3=70;
tau4=5;

B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);

% Lag times
if t<tau1,  g_C6=0;  else
    g_C6=y(2)/(Kg_C6+y(2));end

if t<tau2,  g_C7=0;   else
    g_C7=y(3)/(Kg_C7+y(3));    end

if t<tau3,    g_C8=0;   else
    g_C8 =y(4)/(Kg_C8+y(4));    end

if t<tau4,   g_C10=0;   else
    g_C10=y(5)/(Kg_C10+y(5));    end



% ODE
dydt=[mu*B*(min(f, g_C6)+min(f, g_C7)+min(f, g_C8)+min(f,g_C10))-d*B;
    -1/r*mu*B*min(f,g_C6)+Cin;
    -1/r*mu*B*min(f,g_C7)+Cin;
    -1/r*mu*B*min(f,g_C8)+Cin;
    -1/r*mu*B*min(f,g_C10)+Cin;
    1/r*mu*B*min(f,g_C6);
    1/r*mu*B*min(f,g_C7);
    1/r*mu*B*min(f,g_C8);
    1/r*mu*B*min(f,g_C10)];

end
function dydt1=system51(t1,y1)
% initializations

tau11=15*7;
tau12=10*7;
tau13=10*7;
tau14=2*7;



% Lag times
if t1<tau11,  Kg_C16=0;  else
    Kg_C16=(0.27/10)/7;end

if t1<tau12,  Kg_C17=0;   else
    Kg_C17=(0.47/10)/7;    end

if t1<tau13,    Kg_C18=0;   else
    Kg_C18= (0.74/10)/7;    end

if t1<tau14,   Kg_C110=0;   else
    Kg_C110= (0.41/10)/7;    end

% ODE
dydt1=[-Kg_C16;
    -Kg_C17;
    -Kg_C18;
    -Kg_C110;
    Kg_C16;
    Kg_C17;
    Kg_C18;
    Kg_C110];

end
function dydt2=system52(t2,y2)
% initializations

tau21=15*7;
tau22=17*7;
tau23=12*7;
tau24=5*7;



% Lag times
if t2<tau21,  Kg_C26=0;  else
    Kg_C26=0.09/7;end

if t2<tau22,  Kg_C27=0;   else
    Kg_C27=0.19/7;    end

if t2<tau23,    Kg_C28=0;   else
    Kg_C28= 0.34/7;    end

if t2<tau24,   Kg_C210=0;   else
    Kg_C210= 0.22/7;    end



% ODE
dydt2=[-Kg_C26*y2(1);
    -Kg_C27*y2(2);
    -Kg_C28*y2(3);
    -Kg_C210*y2(4);
    Kg_C26*y2(1);
    Kg_C27*y2(2);
    Kg_C28*y2(3);
    Kg_C210*y2(4)];

end
%%
%%comparison of our model predictiosn with predictions from BTEX methane
%%Tariq et al. 2017
function methane=comparison_BTEX_methane
TM=data.Day;
MM=data.Ch4;
params=1e3*[5.136263100000000,0.000627780000000,0.245673866473879,0.227677675810791,0.059986041033172,0.164435081004224,0.000038252728460];
model=@(vals,TM)subfun(vals, TM);
params0=rand(1,2);
lb=[0,0];
ub=[1,1];
vals = lsqcurvefit(model,params0,TM,MM,lb,ub);
yzero=[params(8),643.472755,918.0739142,983.1050734,852.502112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
[tt, yy] = ode15s(@(tt,yy)system5(tt,yy,vals), [0, TM(end)], yzero);
M=(yy(:,17)*4*0.8)+(yy(:,18)*4.75*0.8)+(yy(:,19)*4*0.8)+(yy(:,20)*4.75*0.8)+(yy(:,21)*4.75*0.8);
h1=figure(1);
set(0,'DefaultAxesFontSize',18,'DefaultTextFontSize',18);
set(gca,'FontSize',18);
plot(TM,MM,'rd','linewidth',5);
errorbar(TM,MM,MMM,'rs','LineWidth',3)
xlabel('time (days)')
ylabel('CH_4 (mmol)')
set(gca,'TickDir','Out');
axis([0 1569 0 2.9])
figure(h1)
line(tt, M,'color','r','LineWidth',3)


    function err = subfun(vals,TM)
        yzero=[params(8),0.29244629,0.112555117,0.146916147,0.244836389,0.132281272,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
        [t,y] = ode15s(@(t,y)system5(t,y,vals), [0,TM(end)], yzero);
        Methan=(y(:,17)*4*0.8)+(y(:,18)*4.75*0.8)+(y(:,19)*4*0.8)+(y(:,20)*4.75*0.8)+(y(:,21)*4.75*0.8);
        yint = interp1(t,Methan,TM);
        err=yint;
    end
end

function dydt=system5(t,y,vals)
params=[4962.3591,0.5608,84.2174,29.075,3392.9329,243.2047,3865.9088,0.00065];
% initializations
Cin=0;
mu=1.3;
d=0;
theta=0.3;
r=0.1;
Tn=params(1);
Kf=params(2);
Kg_npentane=params(3);
Kg_nhexane=params(4);
Kg_methylbutane=params(5);
Kg_methylpentane_2=params(6);
Kg_methylpentane_3=params(7);
B=y(1);
f=(Tn-theta*B)/(Kf+Tn-theta*B);
rho=vals(1);
epsilon=vals(2);




g_npentane=y(2)/(Kg_npentane+y(2));


g_nhexane=y(3)/(Kg_nhexane+y(3));


g_methylbutane =y(4)/(Kg_methylbutane+y(4));


g_methylpentane_2=y(5)/(Kg_methylpentane_2+y(5));


g_methylpentane_3=y(6)/(Kg_methylpentane_3+y(6));


% ODE
dydt=[mu*B*(min(f, g_npentane)+min(f, g_nhexane)+min(f, g_methylbutane)+min(f,g_methylpentane_2)+min(f,g_methylpentane_3))-d*B;
    -1/r*mu*B*min(f,g_npentane)+Cin;
    -1/r*mu*B*min(f,g_nhexane)+Cin;
    -1/r*mu*B*min(f,g_methylbutane)+Cin;
    -1/r*mu*B*min(f,g_methylpentane_2)+Cin;
    -1/r*mu*B*min(f,g_methylpentane_3)+Cin;
    1/r*mu*B*min(f,g_npentane)-rho*y(7);
    1/r*mu*B*min(f,g_nhexane)-rho*y(8);
    1/r*mu*B*min(f,g_methylbutane)-rho*y(9);
    1/r*mu*B*min(f,g_methylpentane_2)-rho*y(10);
    1/r*mu*B*min(f,g_methylpentane_3)-rho*y(11);
    rho*y(7)-epsilon*y(12);
    rho*y(8)-epsilon*y(13);
    rho*y(9)-epsilon*y(14);
    rho*y(10)-epsilon*y(15);
    rho*y(11)-epsilon*y(16);
    epsilon*y(12);
    epsilon*y(13);
    epsilon*y(14);
    epsilon*y(15);
    epsilon*y(16)];
end


%%