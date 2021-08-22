function t = parallel_example(iter)

if nargin==0, iter = 8; end

disp('Start sim')

t0 = tic;
parfor idx = 1:iter
     A(idx) = idx;
     pause(2)
end
t = toc(t0);

disp('Sim Completed')
