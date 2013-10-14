function psych = load_psych(id_data)
    % select folder
    if ~exist('id_data','var') || isnan(id_data)
        datafolder = 'data';
        id_data = nan;
    else
        datafolder = ['data_',num2str(id_data)];
    end

    % subjects' file names
    lsdata = regexp(ls(datafolder),'\s','split');
    i = 1;
    while i<=length(lsdata)
        if isempty(lsdata{i})
            lsdata(i) = [];
        else
            i = i+1;
        end
    end

    % initialize variables
    load([datafolder,filesep,lsdata{1}]);
    nb_psych = max(data.vb_psych);

    % create psych
    psych = [];
    for i_psych = 1:nb_psych
        i_conditions = find(conditions(:,8)==i_psych);
        i_conditions = i_conditions(1);
        psych = [psych; conditions(i_conditions,2),conditions(i_conditions,4),conditions(i_conditions,5),conditions(i_conditions,6) ];
    end
end