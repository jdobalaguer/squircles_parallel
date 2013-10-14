function ret = load_participants(id_data)

    datafolder = ['data_',num2str(id_data)];

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
    nb_subjects = length(lsdata);
    participants = cell(1,nb_subjects);
    datas = cell(1,nb_subjects);
    % loading data
    for i_subjects = 1:nb_subjects
        % load
        data = load([datafolder,filesep,lsdata{i_subjects}]);
        % store participant
        participants{i_subjects} = data.participant;
        % store data
        datas{i_subjects} = data.data;
    end

    % AGE, SEX, AERR
    age = {};
    sex = {};
    aerr = {};
    for i_subjects = 1:nb_subjects
        age{end+1} = participants{i_subjects}.age;
        sex{end+1} = participants{i_subjects}.sex;
        aerr{end+1} = mean(abs(datas{i_subjects}.resp_err));
    end

    ret = struct();
    ret.participants = participants;
    ret.datas = datas;
    ret.age = age;
    ret.sex = sex;
    ret.aerr = aerr;
end