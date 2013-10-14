% fit model_mean error ratios

function ret = analysis_meanfit(sdata,ret)
    
    % load sdata
    if ~exist('sdata','var') || isempty(sdata)
        sdata = load_badgui();
    end
    
    if ~exist('ret','var') || isempty(ret)
        % grid of starting points
        c  = 1:2:12;
        e1 = 0.0:0.05:0.5;
        e2 = 0; %0.0:0.1:0.2;
        nb_c  = length(c);
        nb_e1 = length(e1);
        nb_e2 = length(e2);
        z = nan(nb_c,nb_e1,nb_e2);
        x_c = nan(nb_c,nb_e1,nb_e2);
        x_e1 = nan(nb_c,nb_e1,nb_e2);
        x_e2 = nan(nb_c,nb_e1,nb_e2);

        % precalculate the bias error -----------------------------------------
        biaserr = analysis_errsregression(sdata);
        biaserr(:) = 0;
        
        % build parameter vectors ---------------------------------------------
        x_c = nan( nb_c,nb_e1,nb_e2);
        x_e1 = nan(nb_c,nb_e1,nb_e2);
        x_e2 = nan(nb_c,nb_e1,nb_e2);
        for i_c = 1:nb_c
            x_c(i_c,:,:) = c(i_c);
            for i_e1 = 1:nb_e1
                x_e1(:,i_e1,:) = e1(i_e1);
                for i_e2 = 1:nb_e2
                    x_e2(:,:,i_e2) = e2(i_e2);
                end
            end
        end
        
        % find error for model_mean -------------------------------------------
        nb_cases = nb_c*nb_e1*nb_e2;
        matlabpool open;
        tools_parforprogress(nb_cases);
        parfor i_case = 1:nb_cases
            % parameters
            this_c  = x_c(i_case);
            this_e1 = x_e1(i_case);
            this_e2 = x_e2(i_case);
            %fprintf(['analysis_resamplerfit: ',num2str(i_case),'/',num2str(nb_cases),': (c,errp,erra) = (',num2str(this_c),';',num2str(this_e1),';',num2str(this_e2),')\n']);
            z(i_case) = cost_function(sdata,biaserr,this_c,this_e1,this_e2);
            tools_parforprogress();
        end
        tools_parforprogress(0);
        matlabpool close;
        
        % polynomial regression -----------------------------------------------
        if length(e2)==1
            nb_degrees = 3;  % deg(polynome)
            nb_regressors = sum(1:(nb_degrees+1));
            % x values
            l_x = nb_c*nb_e1;
            xx_c = x_c(:);
            xx_e1 = x_e1(:);
            
            % y value
            y = reshape(z,[1,nb_c*nb_e1]);
            % matrix of regressors
            X = nan(l_x,nb_regressors);
            j_regressor = 0;
            for i_degree = 0:nb_degrees
                for i_regressor = 0:i_degree
                    % counter
                    j_regressor = j_regressor + 1;
                    % x_c
                    this_xc  = power(xx_c,i_regressor);
                    % x_e1
                    this_xe1 = power(xx_e1,i_degree-i_regressor);
                    % regressor
                    this_regressor = this_xc.*this_xe1;
                    % store in X
                    X(:,j_regressor) = this_regressor;
               end
            end
            % find betas
            b = regress(y',X);
            % estimation
            est_y = (X*b)';
            % reshape
            est_z = reshape(est_y,nb_c,nb_e1);
        end
        
        % return
        ret = struct();
        ret.c = c;
        ret.e1 = e1;
        ret.e2 = e2;
        ret.z = z;
        if length(e2)==1
            ret.est_z = est_z;
        end
        ret.x_c = x_c;
        ret.x_e1 = x_e1;
        ret.x_e2 = x_e2;
    end

    if length(ret.e2)==1
        % find minimum in data
        [~,i_z] = min(ret.z(:));
        fprintf(['analysis_meanfit: data\n']);
        fprintf(['analysis_meanfit: min x_c  =',num2str(ret.x_c(i_z)),'\n']);
        fprintf(['analysis_meanfit: min x_e1 =',num2str(ret.x_e1(i_z)),'\n']);
        % find minimum in regression
        [~,i_estz] = min(ret.est_z(:));
        fprintf(['analysis_meanfit: regression\n']);
        fprintf(['analysis_meanfit: min x_c  =',num2str(ret.x_c(i_estz)),'\n']);
        fprintf(['analysis_meanfit: min x_e1 =',num2str(ret.x_e1(i_estz)),'\n']);

        % plot regression
        figure; surf(ret.c,ret.e1,ret.est_z');

        % plot model
        figure; surface(ret.c,ret.e1,ret.z');
    else
        fprintf(['analysis_meanfit: warning. length(ret.e2)==1, nothing to do.\n']);
    end
    
end

% cost function of the fitting
function out = cost_function(sdata,biaserr,capacity,errp,erra)

    % new model estimation
    %mdata = model_cat(model_biaserr(model_mean(sdata,capacity,errp,erra));
    mdata = model_spectrum(model_biaserr(model_mean(sdata,capacity,errp,erra),biaserr));
    
    % load numbers
    load_numbers;
    
    
%    criteria = 'setsize_error';
    criteria = 'response';
    
    switch criteria
        case 'response'
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
        case 'setsize_error'
            % get setsize's error
            nb_times = 50;
            s_resperr = nan(nb_times,nb_subjects,nb_setsizes);
            m_resperr = nan(nb_times,nb_subjects,nb_setsizes);
            for i_times = 1:nb_times
                for i_subject = 1:nb_subjects
                    for i_setsize = 1:nb_setsizes
                        i_trial = (sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_setsize==u_setsize(i_setsize));
                        s_resperr(i_times,i_subject,i_setsize) = mean(sdata.resp_err(i_trial));
                        m_resperr(i_times,i_subject,i_setsize) = mean(mdata.resp_err(i_trial));
                    end
                end
            end
            mean_sresperr = mean(s_resperr,1);
            mean_mresperr = mean(m_resperr,1);
            out = sqrt(mean(power(mean_sresperr(:) - mean_mresperr(:),2)));
    end
end