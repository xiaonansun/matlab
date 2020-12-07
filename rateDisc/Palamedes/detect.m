%This function simulates choices of an observer with a user-defined
%detection curve.
function choice=detect(mu,stimIntensity)

x=[0:0.01:1];%normalized stim range

%mu=0.55;%threshold intensity
beta=0.05;%slope
lambda=0.01;%lapse rate
gamma=0.5;%guess rate

if rand<(-lambda+gamma+(1-gamma-lambda)*exp(-exp(-(stimIntensity-mu)/beta)))
    choice=1;
else
    choice=0;
end

end
