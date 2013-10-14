function allsdata = load_some(id_data,nb_sdata)

    % select folder
    datafolder = ['data_',num2str(id_data)];
    
    % some
    if ~exist('nb_sdata','var')
        nb_sdata = 1;
    end
    
    sdata = load_badgui(id_data);
    load_numbers;
    
    allsdata = struct;
    
    for i_sdata = 0:(nb_sdata-1)
        % change subjects' number
        new_subjects = i_sdata*nb_subjects + (1:nb_subjects);
        this_sdata = change_participants(sdata,new_subjects);
        % concatenate
        if ~i_sdata
            allsdata = this_sdata;
        else
            allsdata = tools_catstruct(allsdata,this_sdata);
        end
    end
    
    % id data
    allsdata.id = id_data;
end

function sdata = change_participants(sdata,new_subjects)
    load_numbers;
    for i_subject = 1:nb_subjects
        i_trials = (sdata.exp_sub==u_subject(i_subject));
        sdata.exp_sub(i_trials) = new_subjects(i_subject);
    end
end