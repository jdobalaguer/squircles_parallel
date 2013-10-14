% NOTE
%{
    only works if means are uniformly distributed
%}


function sdata = model_spectrum(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % load numbers
    load_numbers;
    
    % get closest values --------------------------------------------------
    % calculate limit between
    d_mean = diff(u_mean);
    b_mean = [-inf,u_mean(1:end-1)+.5*d_mean,+inf];
    % find position of responses between boundaries
    [h_respcat,i_respcat] = histc(sdata.resp_cat,b_mean);
    l_bmean = length(b_mean);
    i_respcat(i_respcat==l_bmean) = l_bmean-1;
    % replace with possible means
    sdata.resp_cat = u_mean(i_respcat);
    
    % response ----------------------------------------------------
    sdata.resp_err = sdata.resp_cat - sdata.vb_mean;
    sdata.resp_cor = ~sdata.resp_err;
    if isfield(sdata,'resp_code')
        sdata = rmfield(sdata,'resp_code');
    end
end
