function [ P, A ] = runDR( D, Pg, params )
% Runs dual regression on the data set to give a set of subject specific maps
% and time courses

nComps = size(Pg, 2);

iPg = pinv(Pg - mean(Pg,1));

% Regression for the time courses
A = cell(params.S,1);
for s = 1:params.S
    A{s} = cell(params.R(s),1);
    for r = 1:params.R(s)
        % Extract time courses
        Dsr = D{s}{r} - mean(D{s}{r},1);
        A{s}{r} = iPg * Dsr;
    end
end

% And for the subject maps
P = cell(params.S,1);
for s = 1:params.S
    % Extract data and time courses for this subject
    As = zeros(nComps, params.T*params.R(s));
    Ds = zeros(params.V, params.T*params.R(s));
    for r = 1:params.R(s)
        As(:, (r-1)*params.T+(1:params.T)) = A{s}{r} - mean(A{s}{r},2);
        Ds(:, (r-1)*params.T+(1:params.T)) = D{s}{r} - mean(D{s}{r},2);
    end
    % Pseudo-inverse to get P
    P{s} = Ds * pinv(As);
end

end
