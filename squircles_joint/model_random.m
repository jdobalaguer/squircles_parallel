function sdata = model_random(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % load numbers
    load_numbers;
    
    % response ----------------------------------------------------
    sdata.resp_cat = u_mean(randi(nb_means,1,length(sdata.resp_err)));
    sdata.resp_err = sdata.resp_cat - sdata.vb_mean;
    sdata.resp_cor = ~sdata.resp_err;
    sdata = rmfield(sdata,'resp_code');
end