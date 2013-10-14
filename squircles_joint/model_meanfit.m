% fit model_mean error ratios

function sdata = model_meanfit(sdata,capacity)
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end
    
    % set capacity
    if ~exist('capacity','var')
        capacity = 1;
    end
    
    % fit the error ratios for model_mean
    fit = 1; % use method
    switch fit
        case 1
            % using multiple fminsearch (derivative-free method) ----------
            costf = @(e) cost_function(sdata,capacity,e(1),e(2));
            opt = optimset(...
                            'PlotFcns',{@optimplotx},...
                            'Display','iter',...
                            'MaxIter',20);
            % grid of starting points
            %e1 = 0:0.05:0.5;
            %e2 = 0:0.05:0.5;
            e1 = .175;
            e2 = .275;
            % do it multiple times
            min_e = inf;
            for i_e1 = 1:length(e1)
                for i_e2 = 1:length(e2)
                    % fminsearch
                    this_e = fminsearch(costf,[e1(i_e1),e2(i_e2)],opt);
                    % update min_e
                    if min_e > this_e
                        min_e = this_e;
                    end
                end
            end
            e = min_e;
        case 2
            % using genetic algorithms ------------------------------------
            costf = @(e) cost_function(sdata,capacity,e(1),e(2));
            lb = [0 , 0]; % lower boundaries
            ub = [1 , 1]; % upper boundaries
            opt = gaoptimset(...
                            'PlotFcns',{@gaplotbestf,@gaplotstopping},...
                            'Display','iter',...
                            'Generations',20);
            e = ga(costf,2,[],[],[],[],lb,ub,[],opt);
        otherwise
            error('model_meanfit: error. fit method doesn''t exist.');
    end
    
    % store the fit
    sdata = model_cat(model_mean(sdata,capacity,e(1),e(2)));
    sdata.model_efit = e;
end

function out = cost_function(sdata,capacity,errp,erra)
    % new model estimation
    mdata = model_cat(model_mean(sdata,capacity,errp,erra));
    
    % load numbers
    load_numbers;
    
    % get samples of resp_cat(subject,psych,mean)
    nb_times = 30;
    s_respcat = cell(1,nb_times);
    m_respcat = cell(1,nb_times);
    for i_times = 1:nb_times
        s_respcat{i_times} = nan(nb_subjects,nb_psychs,nb_means);
        m_respcat{i_times} = nan(nb_subjects,nb_psychs,nb_means);
        for i_subject = 1:nb_subjects
            for i_psych = 1:nb_psychs
                for i_mean = 1:nb_means
                    i_trial = (sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_psych==u_psych(i_psych)) & (sdata.vb_mean==u_mean(i_mean));
                    s_respcat{i_times}(i_subject,i_psych,i_mean) = mean(sdata.resp_cat(i_trial));
                    m_respcat{i_times}(i_subject,i_psych,i_mean) = mean(mdata.resp_cat(i_trial));
                end
            end
        end
    end
    
    % estimate the average resp_cat
    sum_srespcat = zeros(nb_subjects,nb_psychs,nb_means);
    sum_mrespcat = zeros(nb_subjects,nb_psychs,nb_means);
    for i_times = 1:nb_times
        sum_srespcat = sum_srespcat + s_respcat{i_times};
        sum_mrespcat = sum_mrespcat + m_respcat{i_times};
    end
    mean_srespcat = sum_srespcat / nb_times;
    mean_mrespcat = sum_mrespcat / nb_times;
    
    % error
    out = sqrt(mean(power(mean_srespcat(:) - mean_mrespcat(:),2)));
end