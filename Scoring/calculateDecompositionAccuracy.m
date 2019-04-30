function [ Pscore, Ascore, cPscore, cAscore ] ...
    = calculateDecompositionAccuracy( P, infP, A, infA, params )
%CALCULATESCORES Returns accuracy of inferred maps and temporal
%correlations
%   Returns the correlation between the subject specific maps and time
%   courses, as well as the error between the spatial and temporal
%   correlation matrices

%%-----------------------------------------------------------------------------
% First calculate the raw correlations of maps / time courses to ground
% truth, and use that to pair components

%Calculate the map accuracies first
%Take cosine similarity between true and inferred maps
cP = 0;
for s = 1:params.S
    cP = cP + [P{s} infP{s}]' * [P{s} infP{s}];
end
cP = corrcov(cP);
%Just the scores between the two
cP = cP( 1:params.N, params.N+(1:params.iN) );

%Repeat for the temporal accuracy
%Take correlation between true and inferred time courses
cA = 0;
for s = 1:params.S
    for r = 1:params.R(s)
        cA = cA + cov([A{s}{r}; infA{s}{r}]');
    end
end
%Normalise
cA = corrcov(cA);
%Just the scores between the two
cA = cA( 1:params.N, params.N+(1:params.iN) );

%Match with the Hungarian algorithm and return scores sorted
matching = Hungarian( -abs( cP + cA ) );

%First find inds that match components
[iGT, iInf] = ind2sub(size(matching), find(matching));
signs = diag(sign( cP(iGT,iInf) + cA(iGT,iInf) ));

%%-----------------------------------------------------------------------------

% Spatial accuracy - cosine similarity between subject maps
Pscore = NaN(params.S, params.iN);
for s = 1:params.S
    cP = [P{s} infP{s}]' * [P{s} infP{s}];
    cP = corrcov(cP);
    % Just the scores between the two
    cP = cP( 1:params.N, params.N+(1:params.iN) );
    
    % Return scores in order of inferred components
    Pscore(s, iInf) = diag(cP(iGT, iInf)) .* signs;
end

%%-----------------------------------------------------------------------------

% Repeat for the temporal accuracy
% Take correlation between true and inferred time courses
Ascore = NaN(sum(params.R), params.iN); sr = 1;
for s = 1:params.S
    for r = 1:params.R(s)
        cA = cov([A{s}{r}; infA{s}{r}]');
        cA = corrcov(cA);
        % Just the scores between the two
        cA = cA( 1:params.N, params.N+(1:params.iN) );
        
        % Return scores in order of inferred components
        Ascore(sr, iInf) = diag(cA(iGT, iInf)) .* signs;
        
        sr = sr + 1;
    end
end

%%-----------------------------------------------------------------------------

%Now find how well the mode interactions have been calculated

%First for the spatial correlations
N = min(params.N, params.iN);
cPscore = NaN(params.S, N*(N-1)/2);
%Find the scores for each subjects corrcoef accuracy
for s = 1:params.S
    
    %Find z-scored cosine similarity of real and observed maps
    cP = P{s}(:,iGT)' * P{s}(:,iGT);
    cPz = r2z(corrcov(cP));
    
    infcP = (signs * signs') .* (infP{s}(:,iInf)' * infP{s}(:,iInf));
    infcPz = r2z(corrcov(infcP));
    
    %Now look at how similar the matrices are
    cPz = cPz(triu(ones(N),1)==1);
    infcPz = infcPz(triu(ones(N),1)==1);
    
    %Save the RMS error between the correlation matrices
    cPscore(s,:) = (cPz - infcPz).^2;
    
end

%cPscore = sqrt( mean(cPscore)' );

%%-----------------------------------------------------------------------------

%Now repeat for the temporal correlations
N = min(params.N, params.iN);
cAscore = NaN(sum(params.R), N*(N-1)/2); sr = 1;
%Find the scores for each subjects corrcoef accuracy
for s = 1:params.S
    
    %Find z-scored correlations of real and observed timecourses
    for r = 1:params.R(s)
        cA = cov(A{s}{r}(iGT,:)');
        infcA = (signs * signs') .* cov(infA{s}{r}(iInf,:)');
        
        %cAz = r2z(corrcov(cA));
        %infcAz = r2z(corrcov(infcA));
        % Tikhonov regularised partials
        cA = inv( corrcov(cA) + 0.01 * eye(N) );
        cAz = r2z( -corrcov(cA) );
        infcA = inv( corrcov(infcA) + 0.01 * eye(N) );
        infcAz = r2z( -corrcov(infcA) );
        
        %Now look at how similar the matrices are
        cAz = cAz(triu(ones(N),1)==1);
        infcAz = infcAz(triu(ones(N),1)==1);
        
        %Save the RMS error between the correlation matrices
        cAscore(sr,:) = (cAz - infcAz).^2;
        
        sr = sr + 1;
    end
end

%cAscore = sqrt( mean(cAscore)' );

%%-----------------------------------------------------------------------------

end
