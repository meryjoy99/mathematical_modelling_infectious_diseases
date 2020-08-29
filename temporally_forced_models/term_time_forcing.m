function [Output] = Program_5_2(beta0,b1,gamma,mu,S0,I0,MaxTime)
%
% 
%
% Program_5_1( beta0, b1, gamma, mu, S0, I0, MaxTime)
%      This is the MATLAB version of program 5.2 from page 171 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
% It is the simple SIR epidemic with corrected term-time forcing of the transmission rate.
% Note: setting b1 too high can cause numerical difficulties.
%
% This code can also be used to generate bifurcation diagrams, by setting
% b1 equal to a vector of seasonality rates. The bifurcation diagram is
% constructed using extrapolated initial conditions. Try:
% Program_5_2(17/13,[0:0.01:0.25],1/13,1/(50*365),1/17,1e-4,20*365);
%


% Sets up default parameters if necessary.
if nargin == 0
   beta0=17/13;
   b1=0.25;
   gamma=1/13.0;
   mu=1/(50*365.0);
   S0=1/17;
   I0=1e-4;
   MaxTime=10*365;
end

% Checks all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if I0<=0 
    error('Initial level of infecteds (%g) is less than or equal to zero',I0);
end

if beta0<=0 
    error('Transmission rate beta0 (%g) is less than or equal to zero',beta0);
end

if b1<0 
    error('Seasonality beta1 (%g) is less than zero',beta1);
end

if b1>=1
    error('Seasonality beta1 (%g) is greater than or equal to one',beta1);
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end

if mu<0 
    error('Birth / Death rate gamma (%g) is less than zero',mu);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end
    
if S0+I0>1
    warning('Initial level of susceptibles+infecteds (%g+%g=%g) is greater than one',S0,I0,S0+I0);
end

if beta0<gamma+mu
    warning('Basic reproductive ratio (R_0=%g) is less than one',beta0/(gamma+mu));
end

S=S0; I=I0; R=1-S-I;

if length(b1)==1
    
    % Calculate Average Effect of Forcing and Correct for it.
    Ave=0;
    for t=0.5:365
        Ave=Ave+(1+b1*Term(t));
    end
    beta0=beta0/(Ave/365);

    % The main iteration
    options = odeset('RelTol', 1e-5);
    [t, pop]=FORCED_ODE([0 MaxTime],[S I R],options,[beta0 b1 gamma mu]);

    T=t/365; S=pop(:,1); I=pop(:,2); R=pop(:,3);

    % plots the graphs with scaled colours
    subplot(3,1,1)
    h=plot(T,S,'-g');
    xlabel 'Time (years)';
    ylabel 'Susceptibles'

    subplot(3,1,2)
    h=plot(T,I,'-r');
    xlabel 'Time (years)';
    ylabel 'Infectious'

    subplot(3,1,3)
    h=plot(T,R,'-k');
    xlabel 'Time (years)';
    ylabel 'Recovereds'
    
    Output=[t,S,I,R];

else
    if MaxTime<3650
        MaxTime=3650;
    end
    Last=[S0 I0 R];
    for loop=1:length(b1)
        B1=b1(loop);
        Ave=0;
        for t=0.5:365
            Ave=Ave+(1+B1*Term(t));
        end
        Beta0=beta0/(Ave/365);

        options = odeset('RelTol', 1e-5);
        [t, pop]=FORCED_ODE([0 MaxTime],[Last],options,[Beta0 B1 gamma mu]);
    
        Last=[pop(end,1) pop(end,2) pop(end,3)];
        Bifur_I(loop,1:10)=interp1(t,pop(:,2),MaxTime-[0:9]*365,'linear');
        
        semilogy(b1(1:loop),Bifur_I(1:loop,:),'.k');
        xlabel 'Seasonality, b_1'
        ylabel 'Level of Infection' 
        set(gca,'XLim',[min(b1) max(b1)]); drawnow;
    end
    
    semilogy(b1,Bifur_I,'.k');
    xlabel 'Seasonality, b_1'
    ylabel 'Level of Infection'
    
    Output=[reshape(b1,length(b1),1) Bifur_I];
end

% Calculates the differential rates used in the integration.

function [t, pop]=FORCED_ODE(Time, Initial, options, parameter);

t=[]; pop=[];
for Year=0:(max(Time-1)/365)
    [t0, pop0]=ode45(@Diff_5_2,Year*365+[0 6],Initial,options,[(parameter(1)-parameter(2)) parameter(3:4)]);
    [t1, pop1]=ode45(@Diff_5_2,Year*365+[6 100],pop0(end,:),options,[(parameter(1)+parameter(2)) parameter(3:4)]);
    [t2, pop2]=ode45(@Diff_5_2,Year*365+[100 115],pop1(end,:),options,[(parameter(1)-parameter(2)) parameter(3:4)]);
    [t3, pop3]=ode45(@Diff_5_2,Year*365+[115 200],pop2(end,:),options,[(parameter(1)+parameter(2)) parameter(3:4)]);
    [t4, pop4]=ode45(@Diff_5_2,Year*365+[200 251],pop3(end,:),options,[(parameter(1)-parameter(2)) parameter(3:4)]);
    [t5, pop5]=ode45(@Diff_5_2,Year*365+[251 300],pop4(end,:),options,[(parameter(1)+parameter(2)) parameter(3:4)]);
    [t6, pop6]=ode45(@Diff_5_2,Year*365+[300 307],pop5(end,:),options,[(parameter(1)-parameter(2)) parameter(3:4)]);
    [t7, pop7]=ode45(@Diff_5_2,Year*365+[307 356],pop6(end,:),options,[(parameter(1)+parameter(2)) parameter(3:4)]);
    [t8, pop8]=ode45(@Diff_5_2,Year*365+[356 365],pop7(end,:),options,[(parameter(1)-parameter(2)) parameter(3:4)]);
    
    t=[t; t0(1:(end-1)); t1(1:(end-1)); t2(1:(end-1)); t3(1:(end-1)); t4(1:(end-1)); t5(1:(end-1)); t6(1:(end-1)); t7(1:(end-1)); t8(1:(end-1))];
    pop=[pop; pop0(1:(end-1),:); pop1(1:(end-1),:); pop2(1:(end-1),:); pop3(1:(end-1),:); pop4(1:(end-1),:); pop5(1:(end-1),:); pop6(1:(end-1),:); pop7(1:(end-1),:); pop8(1:(end-1),:)];
    Initial=pop8(end,:);
end
t=[t;t8(end)];
pop=[pop;pop8(end,:)];


function dPop=Diff_5_2(t,pop, parameter)

beta=parameter(1); gamma=parameter(2); mu=parameter(3);
S=pop(1); I=pop(2); R=pop(3);

dPop=zeros(3,1);

dPop(1)= mu -beta*S*I - mu*S;
dPop(2)= beta*S*I - gamma*I - mu*I;
dPop(3)= gamma*I - mu*R;


function term=Term(t)

t=mod(t,365);

if (t<6) | (t>100 & t<115) | (t>200 & t<251) | (t>300 & t<307) | (t>356)
    term=-1;
else
    term=1;
end