clear all
close all

%% LOAD

% subjects' file names
lsdata = regexp(ls('data'),'\s','split');
i = 1;
while i<=length(lsdata)
    if isempty(lsdata{i})
        lsdata(i) = [];
    else
        i = i+1;
    end
end

% initialize variables
load(['data/',lsdata{1}]);
nb_psych = max(data.psych);
nb_means = length(variables.means);
nb_subjects = length(lsdata);
nb_conditions = size(conditions,1);
nb_trials = nb_conditions/(nb_means*nb_psych);
max_ss = max(variables.setsizes);
smeans = unique([-variables.means,+variables.means]);
nb_smeans = length(smeans);

%% LOADING
alldata = struct();
% split
alldata.cor         = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
alldata.respcat     = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
alldata.rt          = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
alldata.S           = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
alldata.C           = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
alldata.srespcat = zeros(nb_psych,nb_smeans,nb_subjects,nb_trials);
alldata.srespcat(:) = NaN;
% in a row
alldata.Mv_iar      = zeros(1,nb_conditions*nb_subjects);
alldata.Vv_iar      = zeros(1,nb_conditions*nb_subjects);
alldata.C_iar       = zeros(max_ss,nb_conditions*nb_subjects);
alldata.S_iar       = zeros(max_ss,nb_conditions*nb_subjects);
alldata.V_iar       = zeros(max_ss,nb_conditions*nb_subjects);
alldata.Vcat_iar    = zeros(1,nb_conditions*nb_subjects);
alldata.rt_iar      = zeros(1,nb_conditions*nb_subjects);
alldata.cor_iar     = zeros(1,nb_conditions*nb_subjects);
alldata.respcat_iar = zeros(1,nb_conditions*nb_subjects);
alldata.ss_iar      = zeros(1,nb_conditions*nb_subjects);
alldata.Q_iar       = zeros(1,nb_conditions*nb_subjects);
alldata.R_iar       = zeros(1,nb_conditions*nb_subjects);
alldata.T_iar       = zeros(1,nb_conditions*nb_subjects);


for i_subjects = 1:nb_subjects
    % save
    %save(['data/',lsdata{i_subjects}],'-append','psych');
    % load
    load(['data/',lsdata{i_subjects}]);
    % store cor
    for i_psych = 1:nb_psych
        for i_mean = 1:nb_means
            % split
            alldata.cor(i_psych,i_mean,i_subjects,:)        = data.cor(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.respcat(i_psych,i_mean,i_subjects,:)    = data.respcat(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.rt(i_psych,i_mean,i_subjects,:)         = data.RT(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.S(i_psych,i_mean,i_subjects,:)          = data.S(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.C(i_psych,i_mean,i_subjects,:)          = data.C(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.V(i_psych,i_mean,i_subjects,:)          = data.V(data.psych==i_psych & data.Mv==variables.means(i_mean));
            alldata.Vcat(i_psych,i_mean,i_subjects,:)       = data.Vcat(data.psych==i_psych & data.Mv==variables.means(i_mean));
            % in a row
            alldata.Mv_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.Mv;
            alldata.Vv_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.Vv;
            alldata.C_iar(:,((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.C;
            alldata.S_iar(:,((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.S;
            alldata.V_iar(:,((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.V;
            alldata.Vcat_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.Vcat;
            alldata.rt_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.RT;
            alldata.cor_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.cor;
            alldata.respcat_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = data.respcat;
            alldata.ss_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = conditions(:,2);
            Q = abs(data.V-(ones(max_ss,1)*mean(data.V)));
            Qmax = max(Q);
            Qmax(Qmax==0) = 1;
            Q = mean(Q ./ (ones(max_ss,1)*Qmax));
            alldata.Q_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = Q;
            alldata.R_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = mean(data.V>0);
            alldata.T_iar(((i_subjects-1)*nb_conditions+1):(i_subjects*nb_conditions)) = mean((data.V-(ones(max_ss,1)*mean(data.V)))>0);
        end
        for i_smean = 1:nb_smeans
            iii_smean = data.psych==i_psych & (data.Mv.*data.Vcat)==smeans(i_smean);
            a = -1 .* sign(data.SR(iii_smean)-.5) .* data.respcat(iii_smean);
            alldata.srespcat(i_psych,i_smean,i_subjects,1:length(a)) = a;
        end
    end
end
alldata.mcor = nanmean(alldata.cor,4);
alldata.mmcor = nanmean(alldata.mcor,3);
alldata.mrespcat = squeeze(nanmean(alldata.respcat,4));
alldata.mmrespcat = mean(alldata.mrespcat,3);
alldata.msrespcat = squeeze(nanmean(alldata.srespcat,4));
alldata.mmsrespcat = mean(alldata.msrespcat,3);

%% ANOVA performance
disp('ANOVA performance');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = alldata.mcor(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);

%% ANOVA response
disp('ANOVA response');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = alldata.mrespcat(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);


%% ANOVA reaction times
disp('ANOVA reaction times');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = alldata.rt(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);

%% psych plots
psych = load_psych();
u_setsize       = unique(psych(:,1));
u_variance      = unique(psych(:,2));
u_constraint    = unique(psych(:,3));
u_timing        = unique(psych(:,4));

x_means         = variables.means;
xx_means        = min(x_means):.001:max(x_means);
shifted_sigmoid =  @(p,x) p(2)+(p(4)-p(2)).*(1./(1+exp(-p(3)*(x-p(1)))));

% setsize
psych_setsize = [];
fits_setsize = [];
for iu_setsize = 1:length(u_setsize)
    % average of all psych curves with this setsize
    tmp_mmcor = alldata.mmcor(u_setsize(iu_setsize)==psych(:,1),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_setsize(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,0]);
    % store fit
    fits_setsize(end+1,:) = shifted_sigmoid(p,xx_means);
end
figure;
hold on;
plot(x_means,psych_setsize','+')
plot(xx_means,fits_setsize);

% variance
psych_variance = [];
fits_variance = [];
for iu_variance = 1:length(u_variance)
    % average of all psych curves with this variance
    tmp_mmcor = alldata.mmcor(u_variance(iu_variance)==psych(:,2),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_variance(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,0]);
    % store fit
    fits_variance(end+1,:) = shifted_sigmoid(p,xx_means);
end
figure;
hold on;
plot(x_means,psych_variance,'+');
plot(xx_means,fits_variance);

%{
% constraint
psych_constraint = [];
for iu_constraint = 1:length(u_constraint)
    % average of all psych curves with this constraint
    tmp_mmcor = alldata.mmcor(u_constraint(iu_constraint)==psych(:,3),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_constraint(end+1,:) = tmp_mmcor;
end
figure;
plot(variables.means,psych_constraint)

% timing
psych_timing = [];
for iu_timing = 1:length(u_timing)
    % average of all psych curves with this timing
    tmp_mmcor = alldata.mmcor(u_timing(iu_timing)==psych(:,4),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_timing(end+1,:) = tmp_mmcor;
end
figure;
plot(variables.means,psych_timing)
%}

%% LOG-LINEAR REGRESSION
%{
alldata.C_pc = {};
alldata.S_pc = {};
alldata.V_pc = {};
alldata.D_pc = {};
alldata.sortV_pc = {};
alldata.sortD_pc = {};
alldata.respcat_pc = {};
alldata.weights = zeros(max_ss,length(variables.setsizes));
for ss = variables.setsizes
    % values
    i_ss = (alldata.ss_iar==ss);
    alldata.C_pc{end+1}         = alldata.C_iar(1:ss,i_ss);
    alldata.S_pc{end+1}         = alldata.S_iar(1:ss,i_ss);
    alldata.V_pc{end+1}         = alldata.V_iar(1:ss,i_ss);
    alldata.sortV_pc{end+1}     = sort(alldata.V_pc{end});
    for i = ss
        p = 2;
        alldata.D_pc{end+1}         = norm(power(power(abs(alldata.C_iar(i,i_ss)),p)+power(abs(alldata.S_iar(i,i_ss)),p),1/p));
    end
    alldata.sortD_pc{end+1}     = sort(alldata.D_pc{end});
    alldata.respcat_pc{end+1}   = alldata.respcat_iar(i_ss);
    % log-linear regression
    w = glmfit(alldata.sortV_pc{end}',alldata.respcat_pc{end}==+1,'binomial','link','probit');
    alldata.weights(max_ss+1-(1:ss),(ss==variables.setsizes)) = w(2:end);
end
%}

%% CLEAN
clearvars -except alldata;
