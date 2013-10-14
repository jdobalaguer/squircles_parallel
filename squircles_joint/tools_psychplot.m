function tools_psychplot(sdata,mvalues,vbs,fit,toplot)
    use_color = 1;

    % load numbers
    load_numbers;
    nb_vbs = length(vbs);
    
    % load psych
    psych = load_psych(sdata.id);

    % figure
    figure('color',[1 1 1]);
    % change figure size
    f_xy = get(gcf,'position');
    set(gcf,'position',[f_xy(1),f_xy(2),nb_vbs*f_xy(3),f_xy(4)]);
    
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
        
    % select the legend for variables
    pos_legend = {' items',' variance',' contrained',' seconds'};
    c_legend = {};
    for j_vb = 1:size(psych,2)
        c_legend{j_vb} = {};
        u_vb = unique(psych(:,j_vb));
        for k_vb = 1:length(u_vb)
            c_legend{j_vb}{k_vb} = [num2str(u_vb(k_vb)),pos_legend{j_vb}];
        end
    end
        
    for i_vb = 1:nb_vbs
        vb = vbs(i_vb);
        
        % average each condition individualy
        u_vb = unique(psych(:,vb));
        nb_vb = length(u_vb);
        m_mvalues = nan(nb_vb,nb_means,nb_subjects);
        for j_vb = 1:nb_vb
            ii_vb = psych(:,vb)==u_vb(j_vb);
            m_mvalues(j_vb,:,:) = tools_nanmean(mvalues(ii_vb,:,:),1);
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
            colz = cell(1,nb_vb2);
            for i_vb2 = 1:nb_vb2
                colz{i_vb2} = ((i_vb2-1)/nb_vb2)*ones(1,3);
            end
        end
        
        % subplot
        subplot(1,nb_vbs,i_vb);
        hold on;
        
        if fit
            if sdata.id == 4
                % linear
                f = @(p,x) p(1)*x + p(2);
                beta0 = [1,0];
            else
                % plot fitting
                % 1 = lower bound
                % 2 = upper bound
                % 3 = inflection pointclo
                % 4 = hills' slope
                f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
                beta0 = [0 1 4.5 0.5];
            end
            options  = tools_statset('MaxIter',1000);
            for j_vb = 1:nb_vb
                [p,~,~] = tools_nlinfit(1:nb_means,squeeze(tools_nanmean(m_mvalues(j_vb,:,:),3)),f,beta0,options);
                line(1:nb_means,f(p,1:nb_means),'color',colz{j_vb},'linewidth',3);
            end
        end

        % plot values
        % m_mvalues = [nb_vb,nb_means,nb_subjects]
        % pm_mvalues = [nb_subjects,nb_vb,nb_means]
        pm_mvalues = permute(m_mvalues,[3,1,2]);
        tools_dotplot(pm_mvalues,{''},{'o','o','o','o'},colz);

        % x legend
        set(get(gca,'xlabel'),'string','Mean','fontsize',16,'fontname','Arial');
        
        % y legend
        if exist('ylabel','var')
            set(get(gca,'ylabel'),'string',ylabel,'fontsize',16,'fontname','Arial');
        end
        
        % x-axis range
        x_axis = {};
        for i_xaxis = 1:2:nb_means
            x_axis{end+1} = num2str(u_mean(i_xaxis));
        end
        set(gca,'xtick',1:2:nb_means);
        set(gca,'xticklabel',x_axis);
        
        % force y-axis range
        switch toplot
            case 1
            case 2
            case 3
            case 4
                if sdata.id == 4
                    ylim([-.2,.2])
                    set(gca,'ytick',[-.2:.1:+.2]);
                else
                    ylim([0,1])
                    set(gca,'ytick',[0:.25:1]);
                end
            case 5
        end
        
        % conditions legend
        if nb_vb>1
            legend(c_legend{vb},'Location','SouthEast');
            legend show;
        end
    end
end