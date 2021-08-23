function [var_est,r2_est] = variance_components(fPath)
%% Load and reshape data
% disp(strcat(animal,':',rec));
% files_dir=[savedir filesep animal filesep rec filesep];
load([fPath 'AC.mat'],'A','C');

% reshape matrices: neural_flat = [flattened time x num components],
% A_flat = [flattened pixels x num components]
neural_flat=C';
A_flat=reshape(A,size(A,1)*size(A,2),[]);
indstotake=~isnan(A_flat(:,1));
A_flat=A_flat(indstotake,:);

%% LQ decomposition to calculate norms in an efficient way
[~,r]=qr(neural_flat,0);
l=r';

%% Calculate estimate of R^2 and variance of each component
Al=A_flat*l;
% 2 norm of entire data
norm_Al=trace(Al'*Al);
var_est=zeros(size(neural_flat,2),1);
r2_est=zeros(size(neural_flat,2),1);
for i=1:size(neural_flat,2)
    % reconstructed data using just one component
    Ali=A_flat(:,i)*l(i,:);
    resi=Al-Ali; % residual using just one component
    var_est(i)=trace(Ali'*Ali)./norm_Al; %||rec||^2/||data||^2
    r2_est(i)=1-(trace(resi'*resi)./norm_Al); % 1-||res||^2/||data||^2
end

%% Plot
% sort to find highest variance components, and plot top 6
[var_est_sorted,var_est_sort_inds]=sort(var_est,'descend');
figure;
for i=1:6
    subplot(2,3,i);
    imagesc(A(:,:,var_est_sort_inds(i)));
    title(strcat('Variance = ',num2str(100*var_est_sorted(i))));
    axis image off;
end

% sort to find highest R^2 components, and plot top 6
[r2_est_sorted,r2_est_sort_inds]=sort(r2_est,'descend');
figure;
for i=1:6
    subplot(2,3,i);
    imagesc(A(:,:,r2_est_sort_inds(i)));
    title(strcat('R2 = ',num2str(100*r2_est_sorted(i))));
    axis image off;
end
