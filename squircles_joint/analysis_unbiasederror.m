function errfits = analysis_unbiasederror(sdata,pleaseplot)
    if ~exist('pleaseplot','var')
        pleaseplot = 1;
    end

    u_errs = -.4:.05:.4;

    % load corrected sdata
    if ~exist('sdata','var')
        sdata = load_badgui(sdata.id);
        %sdata = load_corrected();
    end

    % load psych
    psych = load_psych(sdata.id);
    nb_psych    = size(psych,1);

    u_subject   = unique(sdata.exp_sub);
    nb_subjects = length(u_subject);
    
    % fit linear regression to each psych/subject
    errfits = nan(nb_psych,nb_subjects,4);
    f = @(p,x) p(1)+p(2)*exp(-p(3)*(power(x-p(4),2)));
    % figure
    if pleaseplot
        figure;
    end
    for i_psych = 1:nb_psych
        for i_subject = 1:nb_subjects
            % fecth response error
            i_trials = ((sdata.exp_sub==u_subject(i_subject)) & (sdata.vb_psych==i_psych));
            errs     = sdata.resp_err(i_trials);
            % distribution
            h_errs   = hist(errs,u_errs);
            h_errs   = h_errs/sum(h_errs);
            % plot histogram
            if pleaseplot
                subplot(nb_psych,nb_subjects,nb_subjects*(i_psych-1)+i_subject);
                hold on;
                bar(u_errs,h_errs)
            end
            % fit
            p0 = [0,1,1,0];
            [p,~,~] = nlinfit(u_errs,h_errs,f,p0);
            errfits(i_psych,i_subject,:) = p;
            % plot fit
            if pleaseplot
                plot(u_errs,f(p,u_errs),'r');
            end
        end
    end
    
    if ~pleaseplot
        return;
    end
    
    % plots
    % 1 = y-bias
    % 2 = y-scale
    % 3 = 1/variance
    % 4 = x-mean
    toplots = 1;
    
    % do the analysis over... (see 'psych' variable)
    % 1 = setsize
    % 2 = variance
    % 3 = constraint
    % 4 = duration
    vbs = [1,4];
    
    for toplot = toplots
        for vb = vbs
            % average each condition individualy
            u_vb = unique(psych(:,vb));
            nb_vb = length(u_vb);
            m_errfits = nan(nb_vb,nb_subjects);
            for i_vb = 1:nb_vb
                ii_vb = psych(:,vb)==u_vb(i_vb);
                m_errfits(i_vb,:) = mean(errfits(ii_vb,:,toplot),1);
            end

            % figure
            figure('color',[1 1 1]);

            % x-axis
            x_axis = cell(1,nb_vb);
            for i_xaxis = 1:nb_vb
                x_axis{i_xaxis} = num2str(u_vb(i_xaxis));
            end
            set(gca,'xtick',1:nb_vb);
            set(gca,'xticklabel',x_axis);

            % plot values
            p_errfits = permute(m_errfits,[2,1]);
            tools_dotplot(p_errfits,{''},{'o'},{[0,0,0]});
        end
    end
end