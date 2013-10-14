%% LOAD

% subjects' file names
lsdata = regexp(ls('olddata_3'),'\s','split');
i = 1;
while i<=length(lsdata)
    if isempty(lsdata{i})
        lsdata(i) = [];
    else
        i = i+1;
    end
end

% initialize variables
load(['olddata_3/',lsdata{1}]);
nb_psych = max(data.vb_psych);
nb_means = length(variables.means);
nb_subjects = length(lsdata);
nb_conditions = size(conditions,1);
nb_trials = nb_conditions/(nb_means*nb_psych);
smeans = unique([-variables.means,+variables.means]);
nb_smeans = length(smeans);

% loading cor
for i_subjects = 1:nb_subjects
    % load ----------------------------------------------------------------
    load(['olddata_3/',lsdata{i_subjects}]);
    
    % SR ------------------------------------------------------------------
    if experiment.SR
        data.resp_cat = -data.resp_cat;
    end
    
    % Vcat ----------------------------------------------------------------
    % find no Vcat
    i_trials = ((data.vb_mean==0) & (data.vb_setsize==1));
    % decide new categories
    Scat = sign(randi(2,1,sum(i_trials))-1.5);
    % fix Scat
    data.sample_Scat(i_trials) = Scat;
    data.resp_cor(i_trials)    = (Scat == data.resp_cat(i_trials));
    data.resp_err(i_trials)    = ~data.resp_cor(i_trials);
    % fix Vcat
    data.sample_V    = data.sample_S;
    data.sample_Vcat = data.sample_Scat;

    % save ----------------------------------------------------------------
    mkdir('data_3');
    save(['data_3',filesep,participant.filename{2},'_',num2str(participant.filename{3}),'.mat'],'conditions','data','experiment','participant','stimulus','variables');
end
