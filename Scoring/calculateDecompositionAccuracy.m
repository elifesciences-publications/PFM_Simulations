function [ Pscore, Ascore, cPscore, cAscore ] ...
    = calculateDecompositionAccuracy( P, infP, A, infA, params )
%CALCULATESCORES Returns accuracy of inferred maps and temporal
%correlations
%   Returns the correlation between the subject specific maps and time
%   courses, as well as the error between the spatial and temporal
%   correlation matrices


% First calculate the raw correlations of maps / time courses to ground
% truth

%Calculate the map accuracies first
cP = 0;
for s = 1:params.S
    %Take correlations between true and inferred maps
    cP = cP + [P{s} infP{s}]' * [P{s} infP{s}];
end
%Normalise
dcP = sqrt( diag(cP) ); cP = cP ./ (dcP * dcP');
%Just the scores between the two
cP = cP( 1:params.N, params.N+(1:params.iN) );

%Repeat for the temporal accuracy
cA = 0;
for s = 1:params.S
    for r = 1:params.R(s)
        %Take correlations between true and inferred time courses
        cA = cA + [A{s}{r}; infA{s}{r}] * [A{s}{r}; infA{s}{r}]';
    end
end
%Normalise
dcA = sqrt( diag(cA) ); cA = cA ./ (dcA * dcA');
%Just the scores between the two
cA = cA( 1:params.N, params.N+(1:params.iN) );


%Match with the Hungarian algorithm and return scores sorted
matching = Hungarian( -abs( cP + cA ) );

%First find inds that match components
[iGT, iInf] = ind2sub(size(matching), find(matching));
signs = diag(diag(sign( cP(iGT,iInf) + cA(iGT,iInf) )));

%Return scores in order of inferred components
Pscore = NaN(params.iN, 1); Ascore = NaN(params.iN, 1);
Pscore(iInf) = diag(cP(iGT, iInf)) .* diag(signs);
Ascore(iInf) = diag(cA(iGT, iInf)) .* diag(signs);
%Pscore(iInf) = cP(matching==1) .* sign(cP(matching==1) + cA(matching==1));
%Ascore(iInf) = cA(matching==1) .* sign(cP(matching==1) + cA(matching==1));


%Now find how well the mode interactions have been calculated

%First for the spatial correlations
N = min(params.N, params.iN);
cPscore = NaN(params.S, N*(N-1)/2);
%Find the scores for each subjects corrcoef accuracy
for s = 1:params.S
    
    %Find correlations of real and observed maps
    cP = P{s}(:,iGT)' * P{s}(:,iGT);
    infcP = signs' * (infP{s}(:,iInf)' * infP{s}(:,iInf)) * signs;
    
    dcP = sqrt( diag(cP) ); cP = cP ./ (dcP * dcP');
    dinfcP = sqrt( diag(infcP) ); infcP = infcP ./ (dinfcP * dinfcP');
    
    %Now look at how similar the matrices are
    cP = cP(triu(ones(N),1)==1);
    infcP = infcP(triu(ones(N),1)==1);
    
    %Save the RMS error between the correlation matrices
    cPscore(s,:) = (cP-infcP).^2;
    
end

%Now repeat for the temporal correlations
N = min(params.N, params.iN);
cAscore = NaN(params.S, N*(N-1)/2);
%Find the scores for each subjects corrcoef accuracy
for s = 1:params.S
    
    cA = 0; infcA = 0;
    for r = 1:params.R(s)
        %Find correlations of real and observed timecourses
        cA = cA + A{s}{r}(iGT,:) * A{s}{r}(iGT,:)';
        infcA = infcA + signs * (infA{s}{r}(iInf,:) * infA{s}{r}(iInf,:)') * signs';
    end
    dcA = sqrt( diag(cA) ); cA = cA ./ (dcA * dcA');
    dinfcA = sqrt( diag(infcA) ); infcA = infcA ./ (dinfcA * dinfcA');
    
    %Now look at how similar the matrices are
    cA = cA(triu(ones(N),1)==1);
    infcA = infcA(triu(ones(N),1)==1);
    
    %Save the RMS error between the correlation matrices
    cAscore(s,:) = (cA-infcA).^2;
    
end

cPscore = sqrt( mean(cPscore)' );
cAscore = sqrt( mean(cAscore)' );

end