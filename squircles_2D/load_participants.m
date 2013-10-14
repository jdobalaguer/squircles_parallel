clear all

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
nb_subjects = length(lsdata);
participants = cell(1,nb_subjects);

% loading data
for i_subjects = 1:nb_subjects
    % load
    data = load(['data/',lsdata{i_subjects}]);
    % store participant
    participants{i_subjects} = data.participant;
end

% AGE and SEX
age = {};
sex = {};
for i_subjects = 1:nb_subjects
    age{end+1} = participants{i_subjects}.age;
    sex{end+1} = participants{i_subjects}.sex;
end

% SEX
clearvars -except participants age sex
