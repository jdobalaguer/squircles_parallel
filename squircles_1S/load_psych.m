function psych = load_psych()

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

    % create psych
    psych = [];
    for i_psych = 1:nb_psych
        i_conditions = find(conditions(:,7)==i_psych);
        i_conditions = i_conditions(1);
        psych = [psych; conditions(i_conditions,2),conditions(i_conditions,3),conditions(i_conditions,4),conditions(i_conditions,5) ];
    end
end