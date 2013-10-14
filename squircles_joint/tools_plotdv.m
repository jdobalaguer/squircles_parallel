function tools_plotdv(sdata,vb,toplot)
    
    use_color = 1;

    if sdata.id~=4
        fprintf('tools_plotdv: error. only for spectrum task.\n');
        return;
    end

    % load ----------------------------------------------------------------
    
    % load numbers
    load_numbers;
    
    % load psych
    psych = load_psych(sdata.id);

    % err numbers
    u_err = unique(round(sdata.resp_err*100)./100);
    nb_errs = length(u_err);
    
    % variable numbers
    u_vb = unique(psych(:,vb));
    nb_vb = length(u_vb);
    
    % fitting function
    f = @(p,x) p(1)*exp(-p(2)*(power(x-p(3),2)));
    
    % calculate dv --------------------------------------------------------
    
    % fit linear regression to each psych/subject
    % mvalues = [nb_psychs,nb_means,nb_subjects]
    mvalues = nan(nb_psychs,nb_means,nb_subjects);
    for i_mean = 1:nb_means
        %figure;
        for i_psych = 1:nb_psychs
            for i_subject = 1:nb_subjects
                % fecth response error
                i_trials = ((sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_mean==u_mean(i_mean)) & (sdata.vb_psych==u_psych(i_psych)));
                errs     = sdata.resp_cat(i_trials);
                % distribution
                h_errs   = hist(errs,u_err);
                h_errs   = h_errs/sum(h_errs);
                % fit
                p0 = [0.1,10,u_mean(i_mean)];
                [p,~,~] = nlinfit(u_err,h_errs,f,p0);
                p(2) = 1./p(2);
                mvalues(i_psych,i_mean,i_subject) = p(toplot);
                % plot
                %{
                subplot(nb_psychs,nb_subjects,nb_subjects*(i_psych-1)+i_subject);
                hold on
                bar(u_err,h_errs);
                plot(u_err,f(p,u_err),'r');
                %}
            end
        end
    end
    
    % average each condition individualy
    % m_mvalues = [nb_vb,nb_means,nb_subjects]
    m_mvalues = nan(nb_vb,nb_means,nb_subjects);
    for j_vb = 1:nb_vb
        ii_vb = (psych(:,vb)==u_vb(j_vb));
        m_mvalues(j_vb,:,:) = nanmean(mvalues(ii_vb,:,:),1);
    end
    
    % plot ----------------------------------------------------------------
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

    % plot values
    % pm_mvalues = [nb_subjects,nb_vb,nb_means]
    pm_mvalues = permute(m_mvalues,[3,1,2]);
    hold on;
    for i_vb = 1:nb_vb
        % mi_pmmvalues = [nb_vb]
        mi_mpmmvalues = squeeze(mean(pm_mvalues(:,i_vb,:),1));
        line(1:nb_means,mi_mpmmvalues,'color',colz{i_vb},'linewidth',3);
    end
    tools_dotplot(pm_mvalues,{''},{'o','o','o','o'},colz);
    
    % x legend
    pos_xlabel = 'Mean';
    set(get(gca,'xlabel'),'string',pos_xlabel,'fontsize',16,'fontname','Arial')
    
    % toplot label
    switch toplot
    % 1 = y-scale
    % 2 = variance
    % 3 = x-mean
        case 1
            ylabel = 'Y-scale';
        case 2
            ylabel = 'Variance of the estimation';
        case 3
            ylabel = 'Average response';
        otherwise
            error(['tools_plotdv: toplot=',num2str(toplot),'is not valid.']);
    end
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
    
    % y-axis range
    ymin = nanmin(m_mvalues(:));
    ymax = nanmax(m_mvalues(:));
    switch toplot
        case 2
            v_ymin = 0:.2:2;
            v_ymax = 0:.2:2;
            [~,i_vymin] = min(abs(ymin-v_ymin));
            [~,i_vymax] = min(abs(ymax-v_ymax));
            round_ymin = v_ymin(i_vymin);
            round_ymax = v_ymax(i_vymax);
            ylim([round_ymin,round_ymax])
            set(gca,'ytick',[round_ymin:.2:round_ymax]);
        case 2
            v_ymin = 0:.02:.1;
            v_ymax = 0:.02:.1;
            [~,i_vymin] = min(abs(ymin-v_ymin));
            [~,i_vymax] = min(abs(ymax-v_ymax));
            round_ymin = v_ymin(i_vymin);
            round_ymax = v_ymax(i_vymax);
            ylim([round_ymin,round_ymax])
            set(gca,'ytick',[round_ymin:.02:round_ymax]);
        case 3
            v_ymin = u_mean;
            v_ymax = u_mean;
            [~,i_vymin] = min(abs(ymin-v_ymin));
            [~,i_vymax] = min(abs(ymax-v_ymax));
            round_ymin = v_ymin(i_vymin);
            round_ymax = v_ymax(i_vymax);
            ylim([round_ymin,round_ymax])
            set(gca,'ytick',[round_ymin:.1:round_ymax]);
    end
    
    % legend
    pos_legend = {' items',' variance',' contrained',' seconds'};
    c_legend = {};
    for j_vb = 1:size(psych,2)
        c_legend{j_vb} = {};
        u_vb = unique(psych(:,j_vb));
        for k_vb = 1:length(u_vb)
            c_legend{j_vb}{k_vb} = [num2str(u_vb(k_vb)),pos_legend{j_vb}];
        end
    end
    legend(c_legend{vb},'Location','SouthEast')
    legend show;
end