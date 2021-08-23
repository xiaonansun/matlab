function data = rateDisc_flipHS(data,cMode)
%code to collapse to hemispheres. second input determines whether the two
%sides are averaged or subtracted from one another.

if nargin == 1 %default is average over both sides
    cMode = 'mean';
end

dSize = size(data);
data = reshape(data,dSize(1),dSize(2),[]);
data = cat(4, data(:,1:round(dSize(2)/2),:), fliplr(data(:,round(dSize(2)/2)+1:end,:)));

if strcmpi(cMode,'mean')
    data = mean(data,4);
elseif strcmpi(cMode,'diff')
    data = diff(data,[],4);
elseif strcmpi(cMode,'sum')
    data = sum(data,4);
end

dSize(2) = size(data,2);
data = reshape(data,dSize); %bring back to original format