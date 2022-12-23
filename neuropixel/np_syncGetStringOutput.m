function outString = np_syncGetStringOutput(inMat)

%% Converts input matrix to string vector
D = inMat';
D = D(:)';
C = num2str(D);
C(isspace(C))=[];

outString = C;