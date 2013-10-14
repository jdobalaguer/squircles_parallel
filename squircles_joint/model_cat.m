function sdata = model_cat(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % load numbers
    load_numbers;
    
    % response ----------------------------------------------------
    % response
    sdata.resp_cat = sign(sdata.resp_cat);
    % fix 0cat trials
    i_trials = (sdata.resp_cat==0);
    sdata.resp_cat(i_trials) = sign(rand(1,sum(i_trials))-.5);
    % correct
    sdata.resp_cor = (sdata.resp_cat==sdata.sample_Vcat);
    % error
    sdata.resp_err = ~sdata.resp_cor;
    % code
    if isfield(sdata,'resp_code')
        sdata = rmfield(sdata,'resp_code');
    end
end