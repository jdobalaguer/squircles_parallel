% LOAD THE DATA. CORRECT FOR BIASES.
% NOTE. ONLY FOR SPECTRUM!!!!
function sdata = load_corrected()

    % load sdata
    sdata = load_badgui;
    
    % exit if non-spectrum
    if length(unique(sdata.resp_cat))<4
        error('load_corrected: error: categorical (non-spectrum) data');
    end

    % load fits
    fits = load_fits;

    u_subject   = unique(sdata.exp_sub);
    nb_subjects = length(u_subject);
    nb_psych    = max(sdata.vb_psych);

    u_means     = unique(sdata.vb_mean);
    nb_means    = length(u_means);

    % fit linear regression to each subject
    for i_subject = 1:nb_subjects
        i_trials = (sdata.exp_sub==u_subject(i_subject));
        sdata.resp_cat(i_trials) = (sdata.resp_cat(i_trials)-fits(i_subject,2))/fits(i_subject,1);
    end
    
    % correct other resp fields
    sdata.resp_err = sdata.resp_cat - sdata.vb_mean;
    sdata.resp_cor = ~sdata.resp_err;
    sdata = rmfield(sdata,'resp_code');

end