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
smeans = unique([-variables.means,+variables.means]);
nb_smeans = length(smeans);


% loading cor
cor = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
respcat = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
rt = zeros(nb_psych,nb_means,nb_subjects,nb_trials);
srespcat = zeros(nb_psych,nb_smeans,nb_subjects,nb_trials);
srespcat(:) = NaN;
for i_subjects = 1:nb_subjects
    % load
    load(['data/',lsdata{i_subjects}]);
    % store cor
    for i_psych = 1:nb_psych
        for i_mean = 1:nb_means
            % pcor
            cor(i_psych,i_mean,i_subjects,:)     = data.cor(data.psych==i_psych & data.Ms==variables.means(i_mean));
            respcat(i_psych,i_mean,i_subjects,:) = data.respcat(data.psych==i_psych & data.Ms==variables.means(i_mean));
            rt(i_psych,i_mean,i_subjects,:)      = data.RT(data.psych==i_psych & data.Ms==variables.means(i_mean));
        end
        for i_smean = 1:nb_smeans
            iii_smean = data.psych==i_psych & (data.Ms.*data.Scat)==smeans(i_smean);
            a = -1 .* sign(data.SR(iii_smean)-.5) .* data.respcat(iii_smean);
            srespcat(i_psych,i_smean,i_subjects,1:length(a)) = a;
        end
    end
end
mcor = mean(cor,4);
mmcor = mean(mcor,3);
mrespcat = squeeze(nanmean(respcat,4));
mmrespcat = mean(mrespcat,3);
msrespcat = squeeze(nanmean(srespcat,4));
mmsrespcat = mean(msrespcat,3);

%% ANOVA performance
disp('ANOVA performance');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = mcor(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);

%% ANOVA response
disp('ANOVA response');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = mrespcat(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);

%% ANOVA reaction times
disp('ANOVA reaction times');
d = zeros(nb_subjects,nb_means*nb_psych);
for i_subjects = 1:nb_subjects
    d_subject = rt(:,:,i_subjects);
    d(i_subjects,:) = d_subject(:);
end
D = [5,3,2];
repanova(d,D);

%% psych plots (response)
u_setsize       = unique(psych(:,1));
u_variance      = unique(psych(:,2));
u_constraint    = unique(psych(:,3));
u_timing        = unique(psych(:,4));

x_means         = smeans;
xx_means        = min(x_means):.001:max(x_means);
shifted_sigmoid =  @(p,x) p(2)+(p(4)-p(2)).*(1./(1+exp(-p(3)*(x-p(1)))));

% setsize
psych_setsize = [];
fits_setsize = [];
for iu_setsize = 1:length(u_setsize)
    % average of all psych curves with this setsize
    tmp_mmcor = mmsrespcat(u_setsize(iu_setsize)==psych(:,1),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_setsize(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,1]);
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
    tmp_mmcor = mmsrespcat(u_variance(iu_variance)==psych(:,2),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_variance(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,1]);
    % store fit
    fits_variance(end+1,:) = shifted_sigmoid(p,xx_means);
end
figure;
hold on;
plot(x_means,psych_variance,'+');
plot(xx_means,fits_variance);


%% psych plots (performance)
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
    tmp_mmcor = mmcor(u_setsize(iu_setsize)==psych(:,1),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_setsize(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,1]);
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
    tmp_mmcor = mmcor(u_variance(iu_variance)==psych(:,2),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_variance(end+1,:) = tmp_mmcor;
    % fit
    p = nlinfit(x_means,tmp_mmcor,shifted_sigmoid,[0,0,1,1]);
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
    tmp_mmcor = mmcor(u_constraint(iu_constraint)==psych(:,3),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_constraint(end+1,:) = tmp_mmcor;
end
figure;
plot(variables.means,psych_constraint,'+')

% timing
psych_timing = [];
for iu_timing = 1:length(u_timing)
    % average of all psych curves with this timing
    tmp_mmcor = mmcor(u_timing(iu_timing)==psych(:,4),:);
    if size(tmp_mmcor,1)>1
        tmp_mmcor = mean(tmp_mmcor);
    end
    % store mcor
    psych_timing(end+1,:) = tmp_mmcor;
end
figure;
plot(variables.means,psych_timing,'+')
%}
