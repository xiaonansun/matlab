function output_matrix= twoP_indexMatrix(input_vector)
% This function generates a matrix that can be used for comparing cell
% types. The input should be a logical vector representing neuron indices 
% (e.g. of a specific cell type). The output matrix will permutations of
% the input logical vector consisting of three columns: (1) all cells
% represented by a column of ones, (2) the input boolean vector representing the
% indices of the sub population, and (3) the logical negation of the input vector

if ~isvector(input_vector)
disp(['Input must be a logical vector, this input is ' num2str(size(input_vector,1)) ' by ' num2str(size(input_vector,1)) '.']);
return
end

inVec = input_vector;

if isrow(inVec)
    inVec = inVec';
end

if ~islogical(inVec)
    inVec = logical(inVec);
end

output_matrix = [ones(length(inVec),1,'logical') inVec ~inVec];