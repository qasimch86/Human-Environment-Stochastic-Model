%% Program for human environment interaction during firewood transportation
clear all
tic
global d BETA
%% Parameters given
Np=10;% Number of Patches
BETA=5*10^-1;%[5*10^-1,5*10^-1,5*10^-1,5*10^-1];
d_high=0.02;
d_low=0.003;
N=100;
dt=1/52;
time=0:dt:N;
Nt=length(time);
for i=1:Np
    K(i)=5000;
    r(i)=0.06*dt;
end
beta=BETA*dt;
%% Parameters varied
Cl=6.75;%Cl(2)=6.75;Cl(3)=6.75;Cl(4)=6.75;
Ct=5;%Ct(2)=5;Ct(3)=5;Ct(4)=5;
C=Cl+Ct;
% a=0.01;b=1/N;c1=0.01;c=0.01;
d_eps=0.05;epsmax=1/3;leps=1/3:d_eps:epsmax;
Nimx=1;
% T_cross1=zeros(length(leps),Nimx);T_cross2=zeros(length(leps),Nimx);
% T_cross3=zeros(length(leps),Nimx);T_cross4=zeros(length(leps),Nimx);
for keps=1:length(leps)
for Ni=1:Nimx
n=0.1;%leps(keps);% Social norms
eps=leps(keps)*dt;
sig=0.1*dt;
f=0.1;
%% Initial conditions
S=zeros(Nt,Np);I=zeros(Nt,Np);R=zeros(Nt,Np);L=zeros(Nt,Np);T=zeros(Nt,Np);P_LT=0;
P_TL=0;P_S=0;P_SI=0;P_IR=0;P_SIT=zeros(1,Np);
num_SIT=zeros(1,Np);num_SI=zeros(Nt,Np);num_OS=zeros(Nt,Np);num_IR=zeros(Nt,Np);num_LT=zeros(Nt,Np);
num_TL=zeros(Nt,Np);U_L=zeros(Nt,Np);U_T=zeros(Nt,Np);NS=zeros(Nt,1);NL=zeros(Nt,Np);NT=zeros(Nt,Np);
tcross=zeros(1,Np);T_cross_avg=zeros(Np,length(leps));
for i=1:Np
    if i==1
        S(1,i)=K(i)-15;
    else
        S(1,i)=K(i);
    end
    I(1,i)=K(i)-S(1,i);
    NS(i)=1000;% Total number of strategists
    NL(1,i)=700;% Local strategists
    NT(1,i)=NS(i)-NL(1,i);% Transport strategists
    L(1,i)=NL(1,i)/NS(i);%Proportion of local strategists
    T(1,i)=NT(1,i)/NS(i);%1-L(1,i);%Proportion of transport strategists
end
tcross=zeros(1,Np)+Nt;
flag(1:Np)=0;
%% Human Environment System Model
for t=1:Nt-1 % Time in weeks
    for i=1:Np % Patch
        P_SI=beta*I(t,i)/K(i); %Forest model - Probability of S to become I
        P_S=r(i)*(1-(S(t,i)+I(t,i))/K(i)); %Forest model - Probability of new S
        P_IR=eps;
        U_L=-Cl+n*(L(t,i)-0.5);
        U_T=-Ct+n*(0.5-L(t,i))-f*I(t,i);
        if U_L>U_T%UL>UT %
            P_LT=0;
            P_TL=sig*(U_L-U_T); %UL-UT);%
        else
            P_TL=0;
            P_LT=sig*(U_T-U_L); %UT-UL);%
        end
        for j=1:Np
            if and(i<=round(Np/2), j<=round(Np/2))
                P_SIT(j)=beta*d_high*T(t,j)*I(t,j)/K(j);% Interaction of Patches. T = 1-L (Transport Strategists = 1- Local Strategists)
            else
                P_SIT(j)=beta*d_low*T(t,j)*I(t,j)/K(j);
            end
        end
        num_SI=binornd(S(t,i),P_SI);
        num_OS=binornd(S(t,i),P_S);
        num_IR=binornd(I(t,i),P_IR);
        num_LT=binornd(NL(t,i),P_LT);
        num_TL=binornd(NT(t,i),P_TL);
        SIT=zeros(1,Np);
        for j=1:Np
            if not(i==j)
                SIT(j)=binornd(S(t,i),P_SIT(j));
            end
        end
        num_SIT=sum(SIT(:));
        if and(num_SIT>0,flag(i)==0)
            tcross(i)=t+1;
            flag(i)=1;
        end
        S(t+1,i)=S(t,i)-num_SI+num_OS-num_SIT; % Susceptables
        I(t+1,i)=I(t,i)+num_SI-num_IR+num_SIT; %Infectives
        R(t+1,i)=R(t,i)+num_IR;

        if I(t+1,i)<0
            I(t+1,i)=0;
        end
        NL(t+1,i)=NL(t,i)-num_LT+num_TL;
        NT(t+1,i)=NT(t,i)+num_LT-num_TL;%NS(i)-NL(t+1,i)
        L(t+1,i)=NL(t+1,i)/NS(i);
        T(t+1,i)=NT(t+1,i)/NS(i);%1-L(t+1,i);
    end
    if sum(flag)==Np
        break;
    end
end

for i=1:Np
    T_cross(i,keps,Ni)=tcross(i);
end
% T_cross1(keps,Ni)=tcross(1);T_cross2(keps,Ni)=tcross(2);
% T_cross3(keps,Ni)=tcross(3);T_cross4(keps,Ni)=tcross(4);

%     if sum(num_SIT(:,k))==0
%         T(k,keps,keps)=t+1;
%         T_cross(k,keps,keps)=t+1;
%     end
end
keps
end
for keps=1:length(leps)
    for i=1:Np
        if isempty(nonzeros(T_cross(i,keps,:)))
            T_cross_avg(i,keps)=Nt;
        else
            T_cross_avg(i,keps)=sum(T_cross(i,keps,:))/length(nonzeros(T_cross(i,keps,:)));% average
        end
    end
end
% for keps=1:length(leps)
%     if isempty(nonzeros(T_cross1(keps,:)))
%         T_cross1_avg(keps)=Nt;
%     else
%         T_cross1_avg(keps)=sum(T_cross1(keps,:))/length(nonzeros(T_cross1(keps,:)));% average
%     end
%     if isempty(nonzeros(T_cross2(keps,:)))
%         T_cross2_avg(keps)=Nt;
%     else
%         T_cross2_avg(keps)=sum(T_cross2(keps,:))/length(nonzeros(T_cross2(keps,:)));% average
%     end
%     if isempty(nonzeros(T_cross3(keps,:)))
%         T_cross3_avg(keps)=Nt;
%     else
%         T_cross3_avg(keps)=sum(T_cross3(keps,:))/length(nonzeros(T_cross3(keps,:)));% average
%     end
%     if isempty(nonzeros(T_cross4(keps,:)))
%         T_cross4_avg(keps)=Nt;
%     else
%         T_cross4_avg(keps)=sum(T_cross4(keps,:))/length(nonzeros(T_cross4(keps,:)));% average
%     end
% end
% figure
% hold on
% for i=1:Np
%     plot(leps,T_cross_avg(i,:)/52,'Color',[1/i i/Np (Np-i)/Np]);
% end
CODE_Figure(leps,T_cross_avg/52,'Time to Crosspatch');%,S2,S3,S4);
toc