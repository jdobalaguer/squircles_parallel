% fit model_time error ratios

function ret = analysis_resamplerfit(sdata)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % grid of starting points
    c  = [1,2,3,4,6,8,10,12];
    st = [.05,.1,.2,.4,.8];
    e1 = [0.01,0.02,0.05,0.07,0.1,0.12,0.15,0.17,0.2];
    e2 = [0.17:0.01:0.23];
    ti = 1:10;
    
    nb_c  = length(c);
    nb_st = length(st);
    nb_e1 = length(e1);
    nb_e2 = length(e2);
    nb_ti = length(ti);
    
    % build parameter vectors ---------------------------------------------
    x_c = nan( nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    x_st = nan(nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    x_e1 = nan(nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    x_e2 = nan(nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    x_ti = nan(nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    for i_c = 1:nb_c
        x_c(i_c,:,:,:,:) = c(i_c);
        for i_st = 1:nb_st
            x_st(:,i_st,:,:,:) = st(i_st);
            for i_e1 = 1:nb_e1
                x_e1(:,:,i_e1,:,:) = e1(i_e1);
                for i_e2 = 1:nb_e2
                    x_e2(:,:,:,i_e2,:) = e2(i_e2);
                    for i_ti = 1:nb_ti
                        x_ti(:,:,:,:,i_ti) = ti(i_ti);
                    end
                end
            end
        end
    end
    
    % return --------------------------------------------------------------
    ret = struct();
    ret.c = c;
    ret.st = st;
    ret.e1 = e1;
    ret.e2 = e2;
    ret.e2 = ti;
    ret.x_c = x_c;
    ret.x_st = x_st;
    ret.x_e1 = x_e1;
    ret.x_e2 = x_e2;
    ret.x_ti = x_ti;
    
    % find error for model_mean -------------------------------------------
    nb_cases = nb_c*nb_st*nb_e1*nb_e2*nb_ti;
    z = nan(nb_c,nb_st,nb_e1,nb_e2,nb_ti);
    
    if ~matlabpool('size'); matlabpool('open'); end
    tools_parforprogress(nb_cases);
    parfor i_case = 1:nb_cases
        % parameters
        this_c  = ret.x_c(i_case);
        this_st = ret.x_st(i_case);
        this_e1 = ret.x_e1(i_case);
        this_e2 = ret.x_e2(i_case);
        this_ti = ret.x_ti(i_case);
        %fprintf(['analysis_resamplerfit: ',num2str(i_case),'/',num2str(nb_cases),': (c,st,errp,erra,ti) = (',num2str(this_c),';',num2str(this_st),';',num2str(this_e1),';',num2str(this_e2),';',num2str(this_ti),')\n']);
        z(i_case) = cost_function(sdata,this_c,this_st,this_e1,this_e2);
        tools_parforprogress();
    end
    tools_parforprogress(0);
    matlabpool('close');
    
    % return --------------------------------------------------------------
    ret = struct();
    ret.c = c;
    ret.st = st;
    ret.e1 = e1;
    ret.e2 = e2;
    ret.e2 = ti;
    ret.x_c = x_c;
    ret.x_st = x_st;
    ret.x_e1 = x_e1;
    ret.x_e2 = x_e2;
    ret.x_ti = x_ti;
    ret.z = z;
end

% cost function of the fitting
function out = cost_function(sdata,capacity,sampling_time,errp,erra)

    % new model estimation
    mdata = model_cat(model_resampler(sdata,capacity,sampling_time,errp,erra));
    %mdata = model_spectrum(model_resampler(sdata,capacity,sampling_time,errp,erra));
    
    % load numbers
    load_numbers;
    
    % get samples of resp_cat(subject,psych,mean)
    s_respcat = nan(nb_subjects,nb_psychs,nb_means);
    m_respcat = nan(nb_subjects,nb_psychs,nb_means);
    for i_subject = 1:nb_subjects
        for i_psych = 1:nb_psychs
            for i_mean = 1:nb_means
                i_trial = (sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_psych==u_psych(i_psych)) & (sdata.vb_mean==u_mean(i_mean));
                s_respcat(i_subject,i_psych,i_mean) = mean(sdata.resp_cat(i_trial));
                m_respcat(i_subject,i_psych,i_mean) = mean(mdata.resp_cat(i_trial));
            end
        end
    end
    % error
    out = sqrt(mean(power(s_respcat(:) - m_respcat(:),2)));
end