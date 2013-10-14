function sdata = model_mean(sdata,capacity,errp,erra)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % set capacity
    if ~exist('capacity','var')
        capacity = 1;
    end
    
    % set perceptual error
    if ~exist('errp','var')
        errp = 0;
    end
    
    % set averaging error
    if ~exist('erra','var')
        erra = 0;
    end
    
    % load numbers
    load_numbers;

    % set samples
    sdata.model_isample = nan(capacity,nb_trials*nb_subjects);
    sdata.model_vsample = nan(capacity,nb_trials*nb_subjects);
    for i_setsize = 1:nb_setsizes
        setsize = u_setsize(i_setsize);
        nb_x = nb_trials*nb_subjects/nb_setsizes;
        nb_y = min(capacity,setsize);
        % find trials
        ii_trial = (sdata.vb_setsize==u_setsize(i_setsize));
        
        % isamples --------------------------------------------------------
        % random permutation of samples
        [~,i_sample] = sort(rand(setsize,nb_x),1); % multiple randperm
        % keep only the capacity ones
        % nan if capacity>setsize
        j_sample = nan(capacity,nb_x);
        j_sample(1:nb_y,:) = i_sample(1:nb_y,:);
        
        sdata.model_isample(:,ii_trial) = j_sample;
        
        % vsamples --------------------------------------------------------
        f_trials = find(ii_trial);
        for f_trial = f_trials
            i_sample = sdata.model_isample(1:nb_y,f_trial);
            v_sample = sdata.sample_V(i_sample,f_trial);
            sdata.model_vsample(1:nb_y,f_trial) = v_sample;
        end
    end
    
    % error ---------------------------------------------------------------
    sdata.model_esample = errp*randn(capacity,nb_trials*nb_subjects);
    
    % noisy sample --------------------------------------------------------
    sdata.model_nsample = sdata.model_vsample + sdata.model_esample;
    
    % response ------------------------------------------------------------
    if capacity>1
        sdata.resp_cat = erra*randn(1,nb_trials*nb_subjects) + nanmean(sdata.model_nsample);
    else
        sdata.resp_cat = erra*randn(1,nb_trials*nb_subjects) + sdata.model_nsample;
    end
    sdata.resp_err = sdata.resp_cat - sdata.vb_mean;
    sdata.resp_cor = ~sdata.resp_err;
    if isfield(sdata,'resp_code')
        sdata = rmfield(sdata,'resp_code');
    end
end
