function tools_plotcondition(sdata,mvalues,vb,toplot)
    
    plot_type = 'error bars';

    % mvalues = [nb_psychs,nb_means,nb_subjects]

    % load numbers
    load_numbers;
    
    % load psych
    psych = load_psych(sdata.id);

    % get variable values
    u_vb = unique(psych(:,vb));
    nb_vb = length(u_vb);
    
    % toplot label
    switch toplot
        case 1
            ylabel = 'Performance (% correct)';
        case 2
            ylabel = 'Error';
        case 3
            ylabel = 'Absolute error';
        case 4
            if sdata.id == 4
                ylabel = 'Response';
            else
                ylabel = 'Probability (Circle)';
            end
        case 5
            ylabel = 'Reaction Times';
        otherwise
            error(['tools_psychplot: toplot=',num2str(toplot),'is not valid.']);
    end
    
    % average each condition individualy
    % m_mvalues = [nb_vb,nb_subjects]
    m_mvalues = nan(nb_vb,nb_subjects);
    for i_vb = 1:nb_vb
        ii_vb = (psych(:,vb)==u_vb(i_vb));
        m_mvalues(i_vb,:) = mean(mean(mvalues(ii_vb,:,:),2),1);
    end
    
    % figure
    hold on;

    % transform values
    % p_mmvalues = [nb_subjects,nb_vb]
    p_mmvalues = permute(m_mvalues,[2,1]);
    % m_pmmvalues = [nb_vb]
    m_pmmvalues = squeeze(mean(p_mmvalues,1));
    
    switch plot_type
        case 'error bars'
            % all circles
            os = cell(1,nb_vb);
            for i_os = 1:nb_vb
                os{i_os} = 'o';
            end
            % plot errorbars
            tools_dotplot(p_mmvalues,{''},os,{[0,0,0]});
            % plot lines
            line(1:nb_vb,m_pmmvalues,'color',[0,0,0],'linewidth',3);
        case 'normal plot'
            % normal plot
            col = .8*ones(1,3);
            plot(1:length(m_pmmvalues),m_pmmvalues,'o-','color',col,'linewidth',3,'markeredgecolor',col,'markeredgecolor',col,'markerfacecolor',col,'markersize',10);
    end
    
    % x legend
    pos_xlabel = {'Setsize','Variance','Contrained','Duration'};
    switch vb
        case 3
        otherwise
            set(get(gca,'xlabel'),'string',pos_xlabel{vb},'fontsize',16,'fontname','Arial')
    end
    
    % y legend
    if exist('ylabel','var')
        set(get(gca,'ylabel'),'string',ylabel,'fontsize',16,'fontname','Arial');
    end

    % x-axis range
    switch vb
        case 3
            all_xaxis = {'Unconstrained','Constrained'};
            x_axis = cell(1,nb_vb);
            for i_xaxis = 1:nb_vb
                x_axis{i_xaxis} = all_xaxis{u_vb(i_xaxis)+1};
            end
            set(gca,'xtick',1:nb_vb);
            set(gca,'xticklabel',x_axis);
        otherwise
            x_axis = cell(1,nb_vb);
            for i_xaxis = 1:nb_vb
                x_axis{i_xaxis} = num2str(u_vb(i_xaxis));
            end
            set(gca,'xtick',1:nb_vb);
            set(gca,'xticklabel',x_axis);
    end
    
    % y-axis range
    switch toplot
        % COR
        case 1
            % EXPERIMENT 4 - SX+TSs
            if sdata.id == 4
                % at random
                plot(1:nb_vb,ones(1,nb_vb)/nb_means,'k--');
                % y-axis
                ymax = max(m_mvalues(:));
                v_ymax = 0:0.1:1;
                [~,i_vymax] = min(abs(ymax-v_ymax));
                round_ymax = v_ymax(i_vymax);
                ylim([.1,round_ymax]);
                set(gca,'ytick',[.1:.05:round_ymax]);
            else
                % y-axis
                ymax = max(m_mvalues(:));
                v_ymax = 0:0.05:1;
                [~,i_vymax] = min(abs(ymax-v_ymax));
                round_ymax = v_ymax(i_vymax);
                ylim([.5,round_ymax])
                set(gca,'ytick',[.5:.1:round_ymax]);
            end
        % AER
        case 3
            % EXPERIMENT 4 - SX+TSs
            if sdata.id == 4
                % y-axis
                ylim([0,.15])
                set(gca,'ytick',[0:.05:.15]);
            end
        % CAT
        case 4
            % EXPERIMENT 4 - SX+TSs
            if sdata.id == 4
                % y-axis
                ylim([-.2,.2])
                set(gca,'ytick',[-.2:.2:+.2]);
            else
                % y-axis
                ylim([0,1])
                set(gca,'ytick',[0:.5:1]);
            end
    end

end