% linear regression of the signed error for each mean
% (to use as a biased policy)

function [e,b] = analysis_errsregression(sdata)
    
    % load data
    load_data
    
    % get errors
    e = mean(mmerr,1);
    
    % calculate the regression
    b = polyfit(e,u_mean,1);
end