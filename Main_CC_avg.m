%% Program for human environment interaction during firewood transportation
clear all
tic
global d BETA Np
%% Parameters given
Np=10;% Number of Patches
BETA=5*10^-1;%[5*10^-1,5*10^-1,5*10^-1,5*10^-1];
d_low=0.003;
d_high=0.02;
N=200;
dt=1/52;
time=0:dt:N;
Nt=length(time);
for i=1:Np
    K(i)=5000;
    eps(i)=0.3*dt;
    r(i)=0.06*dt;
end
beta=BETA*dt;
%% Parameters varied
Cl=6.75;%Cl(2)=6.75;Cl(3)=6.75;Cl(4)=6.75;
Ct=5;%Ct(2)=5;Ct(3)=5;Ct(4)=5;
C=Cl+Ct;
% a=0.01;b=1/N;c1=0.01;c=0.01;
n=0.1;
f=0.1;
b=-1/N;
sig=0.1*dt;%*(exp(b*time));
%% Initial conditions
S=zeros(Nt,Np);
I=zeros(Nt,Np);
L=zeros(Nt,Np);
num_SIT=zeros(Nt,Np);

Nimx=1;

for Ni=1:Nimx
SNi(1:Nt,1:Np,Ni)=0;
INi(1:Nt,1:Np,Ni)=0;
LNi(1:Nt,1:Np,Ni)=0;
R=zeros(Nt,Np);
T=zeros(Nt,Np);
P_LT=zeros(Nt,Np);
P_TL=zeros(Nt,Np);
P_S=zeros(Nt,Np);
P_SI=zeros(Nt,Np);
P_IR=zeros(Nt,Np);
P_SIT=zeros(Nt,Np);
num_SI=zeros(Nt,Np);
num_SITNi(1:Nt,1:Np,Ni)=0;
num_S=zeros(Nt,Np);
num_IR=zeros(Nt,Np);
num_LT=zeros(Nt,Np);
num_TL=zeros(Nt,Np);
U_L=zeros(Nt,Np);
U_T=zeros(Nt,Np);
NS=zeros(Nt,1);
NL=zeros(Nt,Np);
NT=zeros(Nt,Np);
for i=1:Np
    if i==1
        SNi(1,i,Ni)=K(i)-15;
    else
        SNi(1,i,Ni)=K(i);
    end
    INi(1,i,Ni)=K(i)-SNi(1,i,Ni);
    NS(i)=1000;% Total number of strategists
    NL(1,i)=100;% Local strategists
    NT(1,i)=NS(i)-NL(1,i);% Transport strategists
    LNi(1,i,Ni)=NL(1,i)/NS(i);%Proportion of local strategists
    T(1,i)=NT(1,i)/NS(i);%1-L(1,i);%Proportion of transport strategists
end
%% Human Environment System Model
for t=1:Nt-1 % Time in weeks
    for i=1:Np % Patch
        P_SI(t,i)=beta*INi(t,i,Ni)/K(i); %Forest model - Probability of S to become I
        P_S(t,i)=r(i)*(1-(SNi(t,i,Ni)+INi(t,i,Ni))/K(i)); %Forest model - Probability of new S
        P_IR(t,i)=eps(i);
        
        U_L(t,i)=-Cl+n*(LNi(t,i,Ni)-0.5);
        U_T(t,i)=-Ct+n*(0.5-LNi(t,i,Ni))-f*INi(t,i,Ni);
        if U_L(t,i)>U_T(t,i)%UL>UT %
            P_LT(t,i)=0;
            P_TL(t,i)=sig*(U_L(t,i)-U_T(t,i)); %UL-UT);%
        else
            P_TL(t,i)=0;
            P_LT(t,i)=sig*(U_T(t,i)-U_L(t,i)); %UT-UL);%
        end
        
        for j=1:Np
            if and(i<=round(Np/2),j<=Np/2)
                P_SIT(t,j)=beta*d_high*T(t,j)*INi(t,j,Ni)/K(j);% Interaction of Patches. T = 1-L (Transport Strategists = 1- Local Strategists)
            else
                P_SIT(t,j)=beta*d_low*T(t,j)*INi(t,j,Ni)/K(j);
            end
        end

        num_SI(t,i)=binornd(SNi(t,i,Ni),P_SI(t,i));
         num_S(t,i)=binornd(SNi(t,i,Ni),P_S(t,i));
        num_IR(t,i)=binornd(INi(t,i,Ni),P_IR(t,i));
        num_LT(t,i)=binornd(NL(t,i),P_LT(t,i));
        num_TL(t,i)=binornd(NT(t,i),P_TL(t,i));
%         for j=1:Np
%             if not(i==j)
%                num_SIT(t,i)=num_SIT(t,i)+binornd(S(t,i),P_SIT(t,j));
%             end
%         end
        SIT=zeros(1,Np);
        for j=1:Np
            if not(i==j)
                SIT(j)=binornd(SNi(t,i,Ni),P_SIT(t,j));
            end
        end
        num_SITNi(t,i,Ni)=sum(SIT(:));
        SNi(t+1,i,Ni)=SNi(t,i,Ni)-num_SI(t,i)+num_S(t,i)-num_SITNi(t,i,Ni); % Susceptables
        INi(t+1,i,Ni)=INi(t,i,Ni)+num_SI(t,i)-num_IR(t,i)+num_SITNi(t,i,Ni); %Infectives
        R(t+1,i)=R(t,i)+num_IR(t,i);

        if INi(t+1,i,Ni)<0
            INi(t+1,i,Ni)=0;
        end
        NL(t+1,i)=NL(t,i)-num_LT(t,i)+num_TL(t,i);
        NT(t+1,i)=NT(t,i)+num_LT(t,i)-num_TL(t,i);%NS(i)-NL(t+1,i)
        LNi(t+1,i,Ni)=NL(t+1,i)/NS(i);
        T(t+1,i)=NT(t+1,i)/NS(i);%1-L(t+1,i);
    end
end
Ni
end

for t=1:Nt-1 % Time in weeks
    for i=1:Np % Patch
        S(t,i)=sum(SNi(t,i,:))/Nimx;
        I(t,i)=sum(INi(t,i,:))/Nimx;
        L(t,i)=sum(LNi(t,i,:))/Nimx;
      if length(nonzeros(num_SITNi(t,i,:)))>0
        num_SIT(t,i)=sum(num_SITNi(t,i,:))/length(nonzeros(num_SITNi(t,i,:)));
      else
        num_SIT(t,i)=0;
      end
    end
end

CrossTimes=zeros(Nt,Np);
for i=1:Np
    k=1;
    for t=1:Nt-1
        if (num_SIT(t,i)>0)
            CrossTimes(k,i)=(t+1)/52;
            k=k+1;
        end
    end
    CrossN(i)=k-1;%sum(length(nonzeros(num_SIT(:,i))));
    if CrossN(i)==0
        CrossMid(i)=0;
        CrossVar(i)=0;
    else
        CrossMid(i)=sum(CrossTimes(1:k-1,i))/CrossN(i);
        CrossVar(i)=sum((CrossTimes(1:k-1,i)-CrossMid(i)).^2)/CrossN(i);
        CrossSD(i)=sqrt(CrossVar(i));
    end
end

%% Figures for Stochastic Model
CODE_Figure(time(1:Nt),S,'Susceptible');%,S2,S3,S4);
CODE_Figure(time(1:Nt),I,'Infested');%,S2,S3,S4);
CODE_Figure(time(1:Nt),L,'Strategists');%,S2,S3,S4);
CODE_Figure(time(1:Nt),num_SIT*Nimx,'Crosspatch infestation');%,S2,S3,S4);
CODE_Figure(time(1:max(CrossN)),CrossTimes(1:max(CrossN),:),'Crosspatch infestation');%,S2,S3,S4);
CODE_BAR(1:Np,CrossN,CrossMid,CrossVar,CrossSD);
toc