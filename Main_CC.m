%% Program for human environment interaction during firewood transportation
clear all
tic
global d BETA Np
%% Parameters given
Np=10;% Number of Patches
BETA=5*10^-1;%[5*10^-1,5*10^-1,5*10^-1,5*10^-1];
for i=1:Np
    if i<=Np/2
        d(i)=0.02;
    else
        d(i)=0.005;
    end
end

N=300;
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
S(:,:)=zeros(Nt,Np);
I=zeros(Nt,Np);
R=zeros(Nt,Np);
L=zeros(Nt,Np);
T=zeros(Nt,Np);
P_LT=zeros(Nt,Np);
P_TL=zeros(Nt,Np);
P_S=zeros(Nt,Np);
P_SI=zeros(Nt,Np);
P_IR=zeros(Nt,Np);
P_SIT=zeros(Nt,Np);
num_SIT=zeros(Nt,Np);
num_SI=zeros(Nt,Np);
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
        S(1,i)=K(i)-15;
    else
        S(1,i)=K(i);
    end
    I(1,i)=K(i)-S(1,i);
    NS(i)=1000;% Total number of strategists
    NL(1,i)=100;% Local strategists
    NT(1,i)=NS(i)-NL(1,i);% Transport strategists
    L(1,i)=NL(1,i)/NS(i);%Proportion of local strategists
    T(1,i)=NT(1,i)/NS(i);%1-L(1,i);%Proportion of transport strategists
end
%% Human Environment System Model
for t=1:Nt-1 % Time in weeks
    for i=1:Np % Patch
        P_SI(t,i)=beta*I(t,i)/K(i); %Forest model - Probability of S to become I
        P_S(t,i)=r(i)*(1-(S(t,i)+I(t,i))/K(i)); %Forest model - Probability of new S
        P_IR(t,i)=eps(i);
        
        U_L(t,i)=-Cl+n*(L(t,i)-0.5);
        U_T(t,i)=-Ct+n*(0.5-L(t,i))-f*I(t,i);
        if U_L(t,i)>U_T(t,i)%UL>UT %
            P_LT(t,i)=0;
            P_TL(t,i)=sig*(U_L(t,i)-U_T(t,i)); %UL-UT);%
        else
            P_TL(t,i)=0;
            P_LT(t,i)=sig*(U_T(t,i)-U_L(t,i)); %UT-UL);%
        end
        
        for j=1:Np
            P_SIT(t,j)=beta*d(j)*T(t,j)*I(t,j)/K(j);% Interaction of Patches. T = 1-L (Transport Strategists = 1- Local Strategists)
        end

        num_SI(t,i)=binornd(S(t,i),P_SI(t,i));
         num_S(t,i)=binornd(S(t,i),P_S(t,i));
        num_IR(t,i)=binornd(I(t,i),P_IR(t,i));
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
                SIT(j)=binornd(S(t,i),P_SIT(t,j));
            end
        end
        num_SIT(t,i)=sum(SIT(:));
        S(t+1,i)=S(t,i)-num_SI(t,i)+num_S(t,i)-num_SIT(t,i); % Susceptables
        I(t+1,i)=I(t,i)+num_SI(t,i)-num_IR(t,i)+num_SIT(t,i); %Infectives
        R(t+1,i)=R(t,i)+num_IR(t,i);

        if I(t+1,i)<0
            I(t+1,i)=0;
        end
        NL(t+1,i)=NL(t,i)-num_LT(t,i)+num_TL(t,i);
        NT(t+1,i)=NT(t,i)+num_LT(t,i)-num_TL(t,i);%NS(i)-NL(t+1,i)
        L(t+1,i)=NL(t+1,i)/NS(i);
        T(t+1,i)=NT(t+1,i)/NS(i);%1-L(t+1,i);
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
CODE_Figure(time(1:Nt),num_SIT,'Crosspatch infestation');%,S2,S3,S4);
CODE_Figure(time(1:max(CrossN))*52,CrossTimes(1:max(CrossN),:),'Crosspatch infestation');%,S2,S3,S4);
CODE_BAR(1:Np,CrossN,CrossMid,CrossVar,CrossSD);
toc