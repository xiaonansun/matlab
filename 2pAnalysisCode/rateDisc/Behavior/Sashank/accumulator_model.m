function [timeSteps,a,c]= accumulator_model(stim,params,nTimeSteps)

if~exist('nTimeSteps','var')
    nTimeSteps=1000;
end

sigma_a=params(1); %Accumulator noise
sigma_s=params(2); %Stimulus noise
sigma_i=params(3); %Noise in initial value of a
lambda_a=params(4); %Accumulator leak
lambda_s=params(5); %1/Sensory adaptation tau
phi=params(6);%Constant determining if adaptation is facilitating or depressing
B=params(7); %Level of the bound


tMax=1;
dt=1/nTimeSteps;
timeSteps=dt:dt:tMax;

stimL=round((nTimeSteps)*stim{1});% L&R stimulus timestamps
stimR=round((nTimeSteps)*stim{2});%
events=zeros(1,nTimeSteps);

events(stimL(2:end))=events(stimL(2:end))-1;
events(stimR(2:end))=events(stimR(2:end))+1;


a=zeros(1,nTimeSteps);
c=zeros(1,nTimeSteps);
c(1)=1;
a(1)=normrnd(0,sigma_i); %Noise in initial accumulator value


for t=2:nTimeSteps
    %Weiner noise
    dw=normrnd(0,sigma_a*sqrt(dt));
    
    %Stimulus+adaptiation dynamics

    ds=c(t-1)*(events(t)*(normrnd(1,sigma_s)));
    dc=(lambda_s*(1-c(t-1))+(phi-1)*c(t-1)*norm(events(t)))*dt;
    c(t)=c(t-1)+dc;
    
    %Accumulator dynamics
    da=dw+ds+lambda_a*a(t-1)*dt;
    a(t)=a(t-1)+da;
end

end

