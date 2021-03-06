% LOAD EVERYTHING INTO A SINGLE STRUCT
function sdata = load_badgui(id_data)
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

    % subjects
    nb_subjects = length(lsdata);
    for i_subjects = 1:nb_subjects
        % load
        load([datafolder,filesep,lsdata{i_subjects}]);
        % concatenate
        if i_subjects==1
            sdata = data;
        else
            sdata = tools_catstruct(sdata,data);
        end
    end
    
    % resign first experiments
    if any(id_data == [1,2])
        sdata = model_resign(sdata);
    end
    
    % id
    sdata.id = id_data;
end
