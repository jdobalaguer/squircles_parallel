% CALCULATE LOGISTIC REGRESSIONS FOR ALL SUBJECTS
% by subject
% by average

% ---

clc
clear all
close all

% load all data 
load_badgui;

% create variable
nb_setsizes = length(unique(sdata.vb_setsize));
nb_subjects = length(unique(sdata.exp_sub));
max_setsize = max(unique(sdata.vb_setsize));
beta = zeros(nb_setsizes*nb_subjects,max_setsize);

% by setsize
i_beta = 1;
for setsize = unique(sdata.vb_setsize)
    fprintf(['load_log: setsize',num2str(setsize),'\n']);
    i_setsizes = (sdata.vb_setsize==setsize);
    % by subject
    for subject = unique(sdata.exp_sub)
        fprintf(['load_log:   subject',num2str(subject),'\n']);
        i_subjects = (sdata.exp_sub==subject);
        % find trials
        i_trials = i_subjects & i_setsizes;
        V_trials = sdata.sample_V(1:setsize,i_trials);
        y_trials = (sdata.resp_cat(i_trials)>0);
        % sort values
        if setsize>1
            V_trials = sort(V_trials,'ascend');
        end
        % logistic regression
        [b,X,yest] = tools_logregression(y_trials',V_trials');
        beta(i_beta,1:setsize) = b';
        
        % plot average functions
        % create a histogram of X*b, etc...
        figure; hold on;
        plot(X*b,y_trials,'bx');
        plot(X*b,yest,'g.');
        drawnow;
        
        % beta index
        i_beta = i_beta+1;
    end
end