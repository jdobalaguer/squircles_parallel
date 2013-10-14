% variables
u_condition     = unique(sdata.vb_condition);
u_setsize       = unique(sdata.vb_setsize);
u_var           = unique(sdata.vb_var);
u_radi          = unique(sdata.vb_radi);
u_constraint    = unique(sdata.vb_constraint);
u_timing        = unique(sdata.vb_timing);
u_mean          = unique(sdata.vb_mean);
u_psych         = unique(sdata.vb_psych);
u_counter       = unique(sdata.vb_counter);

nb_conditions   = length(u_condition);
nb_setsizes     = length(u_setsize);
nb_vars         = length(u_var);
nb_radis        = length(u_radi);
nb_constraints  = length(u_constraint);
nb_timings      = length(u_timing);
nb_means        = length(u_mean);
nb_psychs       = length(u_psych);
nb_counters     = length(u_counter);

% experiment
u_subject   = unique(sdata.exp_sub);
u_block     = unique(sdata.exp_block);
u_trial     = unique(sdata.exp_trial);

nb_subjects = length(u_subject);
nb_blocks   = length(u_block);
nb_trials   = length(u_trial);