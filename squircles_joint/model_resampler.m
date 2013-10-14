function sdata = model_resampler(sdata,capacity,sampling_time,errp,erra)

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
    max_timing = max(u_timing);
    max_resampling = ceil(max_timing/sampling_time);

    % set samples
    sdata.model_isample = nan(capacity*max_resampling,nb_trials*nb_subjects);
    sdata.model_vsample = nan(capacity*max_resampling,nb_trials*nb_subjects);
    sdata.model_esample = nan(capacity*max_resampling,nb_trials*nb_subjects);
    
    % setsizes ------------------------------------------------------------
    for i_setsize = 1:nb_setsizes
        setsize = u_setsize(i_setsize);
        nb_x = nb_trials*nb_subjects/(nb_setsizes*nb_timings);
        nb_y = min(capacity,setsize);

        % timings ---------------------------------------------------------
        for i_timing = 1:nb_timings
            timing  = u_timing(i_timing);
            
            % find trials
            ii_trial = (sdata.vb_setsize==setsize) & (sdata.vb_timing==timing);
            
            % resampling
            nb_resampling = ceil(timing/sampling_time);
            
            % isamples ------------------------------------------------
            iii_y = [];
            for i_resampling = 1:nb_resampling
                ii_y = (i_resampling-1)*capacity + (1:nb_y);
                % random permutation of samples
                [~,i_sample] = sort(rand(setsize,nb_x),1); % multiple randperm
                sdata.model_isample(ii_y,ii_trial) = i_sample(1:nb_y,:);
                % store all indexes in a vector
                iii_y = [iii_y,ii_y];
            end
            
            % vsamples ------------------------------------------------
            f_trials = find(ii_trial);
            for f_trial = f_trials
                i_sample = sdata.model_isample(iii_y,f_trial);
                v_sample = sdata.sample_V(i_sample,f_trial);
                sdata.model_vsample(iii_y,f_trial) = v_sample;
            end
        end
    end
    
    % error ---------------------------------------------------------------
    sdata.model_esample = errp.*randn(capacity*max_resampling,nb_trials*nb_subjects);
    
    % noisy sample --------------------------------------------------------
    sdata.model_nsample = sdata.model_vsample + sdata.model_esample;
    
    % response ------------------------------------------------------------
    if capacity*max_resampling>1
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
