% load sdata
if ~exist('sdata','var')
    sdata = load_badgui();
end

% load numbers
load_numbers;

% loading cor
cor = zeros(nb_psychs,nb_means,nb_subjects,nb_counters);
err = zeros(nb_psychs,nb_means,nb_subjects,nb_counters);
cat = zeros(nb_psychs,nb_means,nb_subjects,nb_counters);
rt  = zeros(nb_psychs,nb_means,nb_subjects,nb_counters);
for i_subjects = 1:nb_subjects
    % store cor
    for i_psych = 1:nb_psychs
        for i_mean = 1:nb_means
            cor(i_psych,i_mean,i_subjects,:) = sdata.cor(sdata.sub==u_subject(i_subjects) & sdata.psych==u_psych(i_psych) & sdata.Ms==u_mean(i_mean));
            err(i_psych,i_mean,i_subjects,:) = sdata.err(sdata.sub==u_subject(i_subjects) & sdata.psych==u_psych(i_psych) & sdata.Ms==u_mean(i_mean));
            cat(i_psych,i_mean,i_subjects,:) = sdata.respcat(sdata.sub==u_subject(i_subjects) & sdata.psych==u_psych(i_psych) & sdata.Ms==u_mean(i_mean));
            rt(i_psych,i_mean,i_subjects,:)  = sdata.RT( sdata.sub==u_subject(i_subjects) & sdata.psych==u_psych(i_psych) & sdata.Ms==u_mean(i_mean));
        end
    end
end

mcor = mean(cor,4);
mmcor = mean(mcor,3);
mmmcor = squeeze(mean(mmcor,2));

merr = mean(err,4);
mmerr = mean(merr,3);
mmmerr = squeeze(mean(mmerr,2));

aerr = abs(err);
maerr = mean(aerr,4);
mmaerr = mean(maerr,3);
mmmaerr = squeeze(mean(mmaerr,2));

mcat = mean(cat,4);
mmcat = mean(mcat,3);
mmmcat = squeeze(mean(mmcat,2));

mrt = mean(rt,4);
mmrt = mean(mrt,3);
mmmrt = squeeze(mean(mmrt,2));
