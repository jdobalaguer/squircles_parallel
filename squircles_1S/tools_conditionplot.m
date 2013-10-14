function tools_conditionplot(mvalues,vb)

    % mvalues = [nb_psychs,nb_means,nb_subjects]

    % load numbers
    sdata = load_badgui();
    load_numbers;
    
    % load psych
    psych = load_psych();

    % get variable values
    u_vb = unique(psych(:,vb));
    nb_vb = length(u_vb);
    
    % average each condition individualy
    % m_mvalues = [nb_vb,nb_subjects]
    m_mvalues = nan(nb_vb,nb_subjects);
    for i_vb = 1:nb_vb
        ii_vb = (psych(:,vb)==u_vb(i_vb));
        m_mvalues(i_vb,:) = mean(mean(mvalues(ii_vb,:,:),2),1);
    end

    % figure
    figure('color',[1 1 1]);
    hold on;

    % x-axis
    x_axis = cell(1,nb_vb);
    for i_xaxis = 1:nb_vb
        x_axis{i_xaxis} = num2str(u_vb(i_xaxis));
    end
    set(gca,'xtick',1:nb_vb);
    set(gca,'xticklabel',x_axis);

    % all circles
    os = cell(1,nb_vb);
    for i_os = 1:nb_vb
        os{i_os} = 'o';
    end
    
    % plot errorbars
    % p_mmvalues = [nb_subjects,nb_vb]
    p_mmvalues = permute(m_mvalues,[2,1]);
    tools_dotplot(p_mmvalues,{''},os,{[0,0,0]});
    
    % plot lines
    % m_pmmvalues = [nb_vb]
    m_pmmvalues = squeeze(mean(p_mmvalues,1));
    line(1:nb_vb,m_pmmvalues,'color',[0,0,0],'linewidth',3);
end