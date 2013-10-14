%% LOAD

% subjects' file names
lsdata = regexp(ls('olddata_2'),'\s','split');
i = 1;
while i<=length(lsdata)
    if isempty(lsdata{i})
        lsdata(i) = [];
    else
        i = i+1;
    end
end

% initialize variables
load(['olddata_2/',lsdata{1}]);

% create directory
if ~exist('data_2','dir')
    mkdir('data_2');
end

% loading cor
nb_subjects = length(lsdata);
u_setsize = unique(variables.setsizes);
nb_setsizes = length(u_setsize);
for i_subjects = 1:nb_subjects
    % load ----------------------------------------------------------------
    load(['olddata_2/',lsdata{i_subjects}]);
    
    nb_trials = data.i_trial;
    length_block = nb_trials / experiment.nb_blocks;
    o_trials = ones(1,nb_trials);
    % rename fields -------------------------------------------------------
    new_data = struct();
    
    % exp
    new_data.i_trial = data.i_trial;
    new_data.i_block = data.i_block;
    new_data.exp_sub = data.sub;
    new_data.exp_SR  = data.SR;
    new_data.exp_trial = 1:nb_trials;
    new_data.exp_block = floor((0:(nb_trials-1))/length_block)+1;
    
    % screen
    new_data.screen_theta        = data.theta;
    new_data.screen_fix_terror   = data.fix_terror;
    new_data.screen_stim_terror  = data.stim_terror;
    new_data.screen_mask_terror  = data.mask_terror;
    new_data.screen_pause_terror = data.pause_terror;
    
    % resp
    new_data.resp_RT   = data.RT;
    new_data.resp_code = data.respcode;
    new_data.resp_cat  = data.respcat;
    new_data.resp_cor  = data.cor;
    new_data.resp_err  = data.err;
    
    % vb
    new_data.vb_condition  = data.condition;
    new_data.vb_setsize    = data.setsize;
    new_data.vb_radi       = stimulus.radi         * o_trials;
    new_data.vb_var        = data.Vv;
    new_data.vb_constraint = variables.constraints * o_trials;
    new_data.vb_timing     = variables.timings     * o_trials;
    new_data.vb_mean       = data.Mv;
    new_data.vb_psych      = data.psych;
    
    % vb_counter
    new_data.vb_counter   = nan(1,nb_trials);
    u_mean = unique(new_data.vb_mean);
    nb_means = length(u_mean);
    u_psych = unique(new_data.vb_psych);
    nb_psychs = length(u_psych);
    for i_mean = 1:nb_means
        for i_psych = 1:nb_psychs
            i_trials = (new_data.vb_psych==u_psych(i_psych)) & (new_data.vb_mean==u_mean(i_mean));
            new_data.vb_counter(i_trials) = 1:sum(i_trials);
        end
    end
    
    % sample
    new_data.sample_V      = data.V;
    new_data.sample_C      = data.C;
    new_data.sample_S      = data.S;
    for i_setsize = 1:nb_setsizes
        setsize = u_setsize(i_setsize);
        i_trials = (new_data.vb_setsize==setsize);
        new_data.sample_mv(i_trials) = abs(mean(new_data.sample_V(1:setsize,i_trials)));
        new_data.sample_mc(i_trials) = mean(new_data.sample_C(1:setsize,i_trials));
        new_data.sample_ms(i_trials) = mean(new_data.sample_S(1:setsize,i_trials));
        new_data.sample_Vcat(i_trials) = sign(mean(new_data.sample_V(1:setsize,i_trials)));
        new_data.sample_Ccat(i_trials) = sign(new_data.sample_mc(i_trials) - 0.5);
        new_data.sample_Scat(i_trials) = sign(new_data.sample_ms(i_trials) - 1.5);
    end
    
    % CONDITIONS ----------------------------------------------------------
    % conditions_exp1 = [i_condition,setsize,       variance,constraint,timing,i_mean,i_psych,counter]
    % conditions_exp4 = [i_condition,setsize,i_radi,variance,constraint,timing,i_mean,i_psych,counter]
    new_conditions = nan(size(conditions,1),9);
    new_conditions(:,[1:2,4:9]) = conditions;
    new_conditions(:,3) = stimulus.radi;
    
    % FIX RESP --------------------------------------------------------
    new_data.resp_cat = -sign(new_data.exp_SR-.5).*(new_data.resp_code-2);
    new_data.resp_cor = (new_data.resp_cat == new_data.sample_Vcat);
    
    % save ----------------------------------------------------------------
    % rename
    data = new_data;
    conditions = new_conditions;
    % fix filename spaces
    filename = participant.filename{2};
    filename(filename==' ') = [];
    participant.filename{2} = filename;
    % save
    save(['data_2',filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'conditions','data','experiment','participant','stimulus','variables');
end
