%% Set Parameters

nTimeSteps=1000;

sigma_a=0.5; %Accumulator noise S.D in 1 sec
sigma_s=0.5; %Stimulus noise S.D
sigma_i=0.1; %Noise in initial value of a
lambda_a=-0.5; %Accumulator leak, negative is leaky, positive is unstable
lambda_s=1000; %1/Sensory adaptation tau
phi=-10;%Constant determining if adaptation is facilitating (>1) or depressing (<1).
B=8; %Level of the bound

%Pick stimulus from your favorite example trial. This one is from msm23
stim=bhv.stimEvents{1};

%% Single instantiation of the model
[t,a,c]=accumulator_model(stim,[sigma_a,sigma_s,sigma_i,lambda_a,lambda_s,phi,B],nTimeSteps);


figure
subplot(2,1,1)
plot(t,a)
%plot(t,c)


hold on
plot(t,zeros(1,length(t)),'k--')
%Plot bounds if needed
%plot(t,B*ones(1,length(t)),'r--')
%plot(t,-B*ones(1,length(t)),'r--')

hold off
xlim([0,1]);
%ylim([0,1]);
ylim([-10,10]);
xlabel('Time(Seconds)');
ylabel('Accumulator value a(t)');
title('Single instantiation of noisy accumulator');

clear as cs
% Many instantiations of the model, mean.
for iter=1:5
[ts,as(iter,:),cs(iter,:)]=accumulator_model(bhv.stimEvents{2},[sigma_a,sigma_s,sigma_i,lambda_a,lambda_s,phi,B],nTimeSteps);
end

subplot(2,1,2)

plot(ts,as,'color',[0.8,0.8,0.8])
hold on
f=plot(ts,mean(as),'b');
legend(f,'Mean accumulator value')
plot(t,zeros(1,length(t)),'k--')
%plot(t,B*ones(1,length(t)),'r--')
%plot(t,-B*ones(1,length(t)),'r--')

hold off
xlim([0,1]);
%ylim([0,1]);
ylim([-10,10]);
xlabel('Time(Seconds)');
ylabel('Accumulator value a(t)');
title('Dynamics of noisy accumulator');