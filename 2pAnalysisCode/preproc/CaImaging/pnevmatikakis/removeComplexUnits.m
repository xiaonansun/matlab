function [A, C, S, P, YrA,AY,AA] = removeComplexUnits(A, C, S, f, complexTol, P, YrA,AY,AA)
% Modified substantially by MTK. Now handles S, and has a threshold for the
% size of complex values: it "repairs" small complex values (discarding the
% complex part), and discards large ones.

%% Remove neurons with large complex values

% Find complex units
complexC = (max(abs(imag(C)), [], 2) > complexTol);
nComplex = sum(complexC);

% Check that background isn't complex
if any(~isreal(f))
    error('Inferred activity for background is complex');
end

if nComplex > 0
    A = A(:, ~complexC);
    C = C(~complexC, :);
    S = S(~complexC, :);
    
    if isfield(P, 'gn'); P.gn = P.gn(~complexC); end
    if isfield(P, 'b'); P.b = P.b(~complexC); end
    if isfield(P, 'c1'); P.c1 = P.c1(~complexC); end
    if isfield(P, 'neuron_sn'); P.neuron_sn = P.neuron_sn(~complexC); end
    
    if exist('YrA', 'var')
        YrA = YrA(~complexC, :);
    end
    
    if exist('AY', 'var')
        AY = AY(~complexC, :);
    end
    
    if exist('AA', 'var')
        AA = AA(~complexC, :);
    end
    
    %   Y_res(unsaturatedPix,:) = Y_res;
    %   Y_res(saturatedPix,:) = Ysat - Asat*C - bsat*f;
    %
    %   % MTK (to return A)
    %   A = A(:, 1:nr - nComplex);
    %   A(unsaturatedPix, :) = A;
    %   A(saturatedPix, :) = Asat;
    
    fprintf('Removed %d units with complex inferred activity\n', nComplex);
end


%% Repair neurons with small complex values

% MTK
if ~isreal(C)
    A = real(A);
    C = real(C);
    S = real(S);
    
    if exist('YrA', 'var')
        YrA = real(YrA);
    end
    
    if exist('AY', 'var')
        AY = real(AY);
    end
    
    if exist('AA', 'var')
        AA = real(AA);
    end    
    %   Yres = real(Y_res);
    fprintf('Repairing small complex values\n');
end
