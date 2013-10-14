
function sdata = model_unbiaserr(sdata,biaserr)

    % load data
    load_data;
    
    if ~exist('biaserr','var')
        % load error bias
        [sdata.model_biaserr , ~] = analysis_errsregression(sdata);
    else
        sdata.model_biaserr = biaserr;
    end
    
    % get error bias for all trials
    sdata.model_respbias = nan(1,nb_trials*nb_subjects);
    for i_mean = 1:nb_means
        ii_trials = (sdata.vb_mean==u_mean(i_mean));
        sdata.model_respbias(ii_trials) = sdata.model_biaserr(i_mean);
    end
    
    % response
    sdata.resp_cat = sdata.resp_cat - sdata.model_respbias; % add error to the current response
    if sdata.id == 4
        sdata = model_spectrum(sdata);
    else
        sdata.resp_err = sdata.resp_cat - sdata.vb_mean;
        sdata.resp_cor = ~sdata.resp_err;
    end
end
