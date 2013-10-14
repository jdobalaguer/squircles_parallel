function all = load_subject(i_subject)
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
    load(['data/',lsdata{i_subject}]);

    all.conditions  = conditions;
    all.data        = data;
    all.experiment  = experiment;
    all.participant = participant;
    all.stimulus    = stimulus;
    all.variables   = variables;
end
