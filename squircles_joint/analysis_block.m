
function analysis_block(sdata)

    nb_blocks = max(sdata.exp_block);
    nb_psychs = max(sdata.vb_psych);

    % matrix with conditional performance
    c = cell(nb_psychs,nb_psychs);

    for sub = unique(sdata.exp_sub)
        % get data
        blocked_cor = zeros(1,nb_blocks);
        blocked_psy = zeros(1,nb_blocks);
        for blo = 1:nb_blocks
            i = ((sdata.exp_sub==sub)&(sdata.exp_block==blo));
            blocked_cor(blo) = mean(sdata.resp_cor(i));
            blocked_psy(blo) = unique(sdata.vb_psych(i));
        end
        % process data
        for blo = 2:nb_blocks
            c{blocked_psy(blo-1),blocked_psy(blo)}(end+1) = blocked_cor(blo);
        end     
    end

    % transform cell to matrix
    m = nan(nb_psychs,nb_psychs);
    for i_p = 1:nb_psychs
        for j_p = 1:nb_psychs
            m(i_p,j_p) = mean(c{i_p,j_p});
        end
    end

    % load psych
    psych = load_psych(sdata.id);

    % group by setsizes
    setsizes = unique(sdata.vb_setsize);
    nb_setsizes = length(setsizes);
    ms = nan(nb_setsizes,nb_setsizes);
    cs = cell(nb_setsizes,nb_setsizes);
    for i_ss1 = 1:nb_setsizes
        i1_m = (psych(:,1)==setsizes(i_ss1));
        for i_ss2 = 1:nb_setsizes
            i2_m = (psych(:,1)==setsizes(i_ss2));
            m_ss = m(i1_m,i2_m);
            c_ss = c(i1_m,i2_m);
            c_ss = [c_ss{:}];
            cs{i_ss1,i_ss2} = c_ss;
            ms(i_ss1,i_ss2) = mean(c_ss);
        end
    end

    % run ANOVAs
    ss0 = [];
    ss1 = [];
    y      = [];
    for i_ss1 = 1:nb_setsizes
        for i_ss2 = 1:nb_setsizes
            for i_cs = 1:length(cs{i_ss1,i_ss2})
                ss1(end+1) = setsizes(i_ss1);
                ss0(end+1) = setsizes(i_ss2);
                y(end+1)      = cs{i_ss1,i_ss2}(i_cs);
            end
        end
    end
    p = anovan(y,{ss1,ss0})
    
end