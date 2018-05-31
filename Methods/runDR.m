function [ Ps, As ] = runDR( D, A, params )
%Runs dual regression on the data set given a set of subject specific time
%courses

nComps = size(A{1}{1}, 1);

%Regression to get subject maps
Ps = cell(params.S,1);
for s = 1:params.S
    %Extract data and time courses for this subject
    As = zeros(nComps, params.T*params.R(s));
    Ds = zeros(params.V, params.T*params.R(s));
    for r = 1:params.R(s)
        As(:, (r-1)*params.T+(1:params.T)) = A{s}{r};
        Ds(:, (r-1)*params.T+(1:params.T)) = D{s}{r};
    end
    %Pseudo-inverse to get Ps
    Ps{s} = Ds * As' / (As * As');
end

%And repeat for the time courses
As = cell(params.S,1);
for s = 1:params.S
    As{s} = cell(params.R(s),1);
    %Pseudo-inverse
    iPs = (Ps{s}' * Ps{s}) \ Ps{s}';
    for r = 1:params.R(s)
        %Extract time courses
        As{s}{r} = iPs * D{s}{r};
    end
end

end