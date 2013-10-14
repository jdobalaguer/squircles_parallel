% fit model_time error ratios

function ret = analysis_timefit(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % grid of starting points
    c  = [1,2,3,4,6,9,12];
    tp = 0:.2:3;
    e1 = 0.0:0.1:0.5;
    e2 = 0.0:0.02:0.1;
    nb_c  = length(c);
    nb_tp = length(tp);
    nb_e1 = length(e1);
    nb_e2 = length(e2);
    z = nan(nb_c,nb_tp,nb_e1,nb_e2);
    
    % build parameter vectors ---------------------------------------------
    x_c = nan( nb_c,nb_tp,nb_e1,nb_e2);
    x_tp = nan(nb_c,nb_tp,nb_e1,nb_e2);
    x_e1 = nan(nb_c,nb_tp,nb_e1,nb_e2);
    x_e2 = nan(nb_c,nb_tp,nb_e1,nb_e2);
    for i_c = 1:nb_c
        x_c(i_c,:,:,:) = c(i_c);
        for i_tp = 1:nb_tp
            x_tp(:,i_tp,:,:) = tp(i_tp);
            for i_e1 = 1:nb_e1
                x_e1(:,:,i_e1,:) = e1(i_e1);
                for i_e2 = 1:nb_e2
                    x_e2(:,:,:,i_e2) = e2(i_e2);
                end
            end
        end
    end
    x_c  = x_c(:);
    x_tp = x_tp(:);
    x_e1 = x_e1(:);
    x_e2 = x_e2(:);

    % find error for model_mean -------------------------------------------
    nb_cases = nb_c*nb_tp*nb_e1*nb_e2;
    matlabpool open;
    tools_parforprogress(nb_cases);
    parfor i_case = 1:nb_cases
        % parameters
        this_c  = x_c(i_case);
        this_tp = x_tp(i_case);
        this_e1 = x_e1(i_case);
        this_e2 = x_e2(i_case);
        %fprintf(['analysis_resamplerfit: ',num2str(i_case),'/',num2str(nb_cases),': (c,tp,errp,erra) = (',num2str(this_c),';',num2str(this_tp),';',num2str(this_e1),';',num2str(this_e2),')\n']);
        z(i_case) = cost_function(sdata,this_c,this_tp,this_e1,this_e2);
        tools_parforprogress();
    end
    tools_parforprogress(0);
    matlabpool close;
    
    % return --------------------------------------------------------------
    ret = struct();
    ret.c = c;
    ret.tp = tp;
    ret.e1 = e1;
    ret.e2 = e2;
    ret.x_c = x_c;
    ret.x_tp = x_tp;
    ret.x_e1 = x_e1;
    ret.x_e2 = x_e2;
    ret.z = z;
end

% cost function of the fitting
function out = cost_function(sdata,capacity,timepower,errp,erra)

    % new model estimation
    mdata = model_cat(model_time(sdata,capacity,timepower,errp,erra));
    %mdata = model_spectrum(model_time(sdata,capacity,timepower,errp,erra));
    
    % load numbers
    load_numbers;
    
    % get samples of resp_cat(subject,psych,mean)
    nb_times = 5;
    s_respcat = nan(nb_times,nb_subjects,nb_psychs,nb_means);
    m_respcat = nan(nb_times,nb_subjects,nb_psychs,nb_means);
    for i_times = 1:nb_times
        for i_subject = 1:nb_subjects
            for i_psych = 1:nb_psychs
                for i_mean = 1:nb_means
                    i_trial = (sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_psych==u_psych(i_psych)) & (sdata.vb_mean==u_mean(i_mean));
                    s_respcat(i_times,i_subject,i_psych,i_mean) = mean(sdata.resp_cat(i_trial));
                    m_respcat(i_times,i_subject,i_psych,i_mean) = mean(mdata.resp_cat(i_trial));
                end
            end
        end
    end
    % estimate the average resp_cat
    mean_srespcat = mean(s_respcat,1);
    mean_mrespcat = mean(m_respcat,1);
    % error
    out = sqrt(mean(power(mean_srespcat(:) - mean_mrespcat(:),2)));
end