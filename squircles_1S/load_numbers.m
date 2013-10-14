% variables
u_condition    = unique(sdata.condition);
u_setsize       = unique(sdata.setsize);
u_mean          = unique(sdata.Ms);
u_psych         = unique(sdata.psych);

nb_conditions   = length(u_condition);
nb_setsizes     = length(u_setsize);
nb_means        = length(u_mean);
nb_psychs       = length(u_psych);

% experiment
u_subject   = unique(sdata.sub);
nb_subjects = length(u_subject);

nb_blocks   = unique(sdata.i_block);
nb_trials   = unique(sdata.i_trial);
u_block     = 1:nb_blocks;
u_trial     = 1:nb_trials;

u_counter = 1:25;
nb_counters = 25;