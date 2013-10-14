function sdata = model_abs(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % load numbers
    load_numbers;
    
    % response ----------------------------------------------------
    sdata.model_absmean = abs(sdata.vb_mean);
    % response
    sdata.resp_cat = abs(sdata.resp_cat);
    % correct
    sdata.resp_err = (sdata.resp_cat-sdata.model_absmean);
    % error
    sdata.resp_cor = ~sdata.resp_err;
    % code
    if isfield(sdata,'resp_code')
        sdata = rmfield(sdata,'resp_code');
    end
end