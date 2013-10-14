function analysis_anova(sdata)

    % load data
    load_data;
    
    % set effects
    switch sdata.id
        case 1
            % mean/ss/var
            D = [9,3,2];
        case 2
            % mean/ss/var
            D = [9,3,2];
        case 3
            % mean/ss/timing
            D = [9,4,4];
        case 4
            % mean/ss/timing
            D = [9,4,4];
        otherwise
            error('analysis_anova: experiment id not specified');
    end            

    %% ANOVA performance
    disp('ANOVA performance');
    d = zeros(nb_subjects,nb_means*nb_psychs);
    for i_subjects = 1:nb_subjects
        d_subject = mcor(:,:,i_subjects);
        d(i_subjects,:) = d_subject(:);
    end
    tools_repanova(d,D);

    %% ANOVA error
    disp('ANOVA error');
    d = zeros(nb_subjects,nb_means*nb_psychs);
    for i_subjects = 1:nb_subjects
        d_subject = merr(:,:,i_subjects);
        d(i_subjects,:) = d_subject(:);
    end
    tools_repanova(d,D);

    %% ANOVA abs error
    disp('ANOVA absolute error');
    d = zeros(nb_subjects,nb_means*nb_psychs);
    for i_subjects = 1:nb_subjects
        d_subject = maerr(:,:,i_subjects);
        d(i_subjects,:) = d_subject(:);
    end
    tools_repanova(d,D);

    %% ANOVA response
    disp('ANOVA response');
    d = zeros(nb_subjects,nb_means*nb_psychs);
    for i_subjects = 1:nb_subjects
        d_subject = mcat(:,:,i_subjects);
        d(i_subjects,:) = d_subject(:);
    end
    tools_repanova(d,D);

    %% ANOVA reaction times
    disp('ANOVA reaction times');
    d = zeros(nb_subjects,nb_means*nb_psychs);
    for i_subjects = 1:nb_subjects
        d_subject = mrt(:,:,i_subjects);
        d(i_subjects,:) = d_subject(:);
    end
    tools_repanova(d,D);
    
end