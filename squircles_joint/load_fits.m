function fits = load_fits(sdata)

    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % load numbers
    load_numbers;
    
    % fit linear regression to each psych/subject
    fits = nan(nb_subjects,2);
    for i_subject = 1:nb_subjects
        m_resp = nan(1,nb_means);
        % average response for each psych/subject
        for i_mean = 1:nb_means
            i_trials = ((sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_mean==u_mean(i_mean)));
            m_resp(i_mean) = mean(sdata.resp_cat(i_trials));
        end
        % store regression
        fits(i_subject,:) = polyfit(u_mean,m_resp,1);
    end
end