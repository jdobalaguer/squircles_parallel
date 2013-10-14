% load sdata
if ~exist('sdata','var')
    sdata = load_badgui();
end

% load numbers
load_numbers;

% loading cor
mcor  = zeros(nb_psychs,nb_means,nb_subjects);
merr  = zeros(nb_psychs,nb_means,nb_subjects);
maerr = zeros(nb_psychs,nb_means,nb_subjects);
mcat  = zeros(nb_psychs,nb_means,nb_subjects);
mrt   = zeros(nb_psychs,nb_means,nb_subjects);
for i_subjects = 1:nb_subjects
    % store cor
    for i_psych = 1:nb_psychs
        for i_mean = 1:nb_means
            mcor(i_psych,i_mean,i_subjects)  = mean(sdata.resp_cor(sdata.exp_sub==u_subject(i_subjects) & sdata.vb_psych==u_psych(i_psych) & sdata.vb_mean==u_mean(i_mean)));
            merr(i_psych,i_mean,i_subjects)  = mean(sdata.resp_err(sdata.exp_sub==u_subject(i_subjects) & sdata.vb_psych==u_psych(i_psych) & sdata.vb_mean==u_mean(i_mean)));
            maerr(i_psych,i_mean,i_subjects) = mean(abs(sdata.resp_err(sdata.exp_sub==u_subject(i_subjects) & sdata.vb_psych==u_psych(i_psych) & sdata.vb_mean==u_mean(i_mean))));
            mcat(i_psych,i_mean,i_subjects)  = mean(sdata.resp_cat(sdata.exp_sub==u_subject(i_subjects) & sdata.vb_psych==u_psych(i_psych) & sdata.vb_mean==u_mean(i_mean)));
            mrt(i_psych,i_mean,i_subjects)   = mean(sdata.resp_RT( sdata.exp_sub==u_subject(i_subjects) & sdata.vb_psych==u_psych(i_psych) & sdata.vb_mean==u_mean(i_mean)));
        end
    end
end

mmcor = mean(mcor,3);
mmmcor = squeeze(mean(mmcor,2));

mmerr = mean(merr,3);
mmmerr = squeeze(mean(mmerr,2));

mmaerr = mean(maerr,3);
mmmaerr = squeeze(mean(mmaerr,2));

mmcat = mean(mcat,3);
mmmcat = squeeze(mean(mmcat,2));

mmrt = mean(mrt,3);
mmmrt = squeeze(mean(mmrt,2));
