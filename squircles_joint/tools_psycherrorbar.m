function tools_psycherrorbar(sdata,mvalues,vbs,fit)

    use_color = 1;

    % load numbers
    load_numbers;
    
    % load psych
    psych = load_psych(sdata.id);

    % get variable values
    vb1 = vbs(1);
    vb2 = vbs(2);
    u_vb1 = unique(psych(:,vb1));
    u_vb2 = unique(psych(:,vb2));
    nb_vb1 = length(u_vb1);
    nb_vb2 = length(u_vb2);
    
    % average each condition individualy
    % m_mvalues = [nb_vb1,nb_vb2,nb_means,nb_subjects]
    m_mvalues = nan(nb_vb1,nb_vb2,nb_means,nb_subjects);
    for i_vb1 = 1:nb_vb1
        for i_vb2 = 1:nb_vb2
            ii_vb = (psych(:,vb1)==u_vb1(i_vb1)) & (psych(:,vb2)==u_vb2(i_vb2));
            m_mvalues(i_vb1,i_vb2,:,:) = mean(mvalues(ii_vb,:,:),1);
        end
    end
    
    % select the legend for all variables
    pos_legend = {' items',' variance',' contrained',' seconds'};
    c_legend = {};
    for i_vb = 1:size(psych,2)
        c_legend{i_vb} = {};
        u_vb = unique(psych(:,i_vb));
        for j_vb = 1:length(u_vb)
            c_legend{i_vb}{j_vb} = [num2str(u_vb(j_vb)),pos_legend{i_vb}];
        end
    end

    % colours
    if use_color
        colz = cell(1,5);
        colz{1} = 'b';
        colz{2} = 'g';
        colz{3} = 'r';
        colz{4} = 'c';
        colz{5} = 'm';
        colz{6} = 'y';
    else
        colz = {};
        for i_vb2 = 1:nb_vb2
            colz{i_vb2} = ((i_vb2-1)/nb_vb2)*ones(1,3);
        end
    end

    % figure
    figure('color',[1 1 1]);
    hold on;

    % x-axis
    x_axis = cell(1,nb_vb1);
    for i_xaxis = 1:nb_vb1
        x_axis{i_xaxis} = num2str(u_vb1(i_xaxis));
    end
    set(gca,'xtick',1:nb_vb1);
    set(gca,'xticklabel',x_axis);

    if fit
        % plot fitting
        % 1 = lower bound
        % 2 = upper bound
        % 3 = inflection pointclo
        % 4 = hills' slope
        f = @(p,x) p(1)*x + p(2);
        beta0 = [1 0];
        options  = statset('MaxIter',1000);
        for i_vb2 = 1:nb_vb2
            [p,~,~] = nlinfit(1:nb_means,squeeze(mean(mean(m_mvalues(:,i_vb2,:,:),3),1)),f,beta0,options);
            line(1:nb_means,f(p,1:nb_means),'color',colz{i_vb2},'linewidth',3);
        end
    end
        
    % plot values
    % pm_mvalues = [nb_subjects,nb_vb2,nb_vb1,nb_means]
    pm_mvalues = permute(m_mvalues,[4,2,1,3]);
    
    % all circles
    os = cell(1,nb_vb2);
    for i_os = 1:nb_vb2
        os{i_os} = 'o';
    end
    
    % plot errorbars
    % m_pmmvalues = [nb_subjects,nb_vb2,nb_vb1]
    m_pmmvalues = mean(pm_mvalues,4);
    tools_dotplot(m_pmmvalues,{''},os,colz);
    
    % plot lines
    for i_vb2 = 1:nb_vb2
        % mi_pmmvalues = [nb_vb1]
        mi_mpmmvalues = squeeze(mean(m_pmmvalues(:,i_vb2,:),1));
        line(1:nb_vb1,mi_mpmmvalues,'color',colz{i_vb2},'linewidth',3);
    end
    
    % x label
    pos_xlabel = {'setsize','variance','contrained','duration'};
    set(get(gca,'xlabel'),'string',pos_xlabel{vb1})
    
    % legend
    legend(c_legend{vb2},'Location','SouthEast')
    legend show;
end