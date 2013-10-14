function analysis_slopes(sdata)
    
    % load sdata
    if ~exist('sdata','var')
        sdata = load_badgui();
    end

    % load fits
    fits = load_psychfits(sdata);
    
    % load psychc
    psych = load_psych(sdata.id);
    
    nb_subjects = size(fits,2);

    for i_slopes = 1:2
        
        % extract slopes
        slopes = fits(:,:,i_slopes);

        % do the analysis over... (see 'psych' variable)
        % 1 = setsize
        % 2 = variance
        % 3 = constraint
        % 4 = duration
        vbs = [1,4];
        for vb = vbs
            % average each condition individualy
            u_vb = unique(psych(:,vb));
            nb_vb = length(u_vb);
            m_slopes = nan(nb_vb,nb_subjects);
            for i_vb = 1:nb_vb
                ii_vb = psych(:,vb)==u_vb(i_vb);
                m_slopes(i_vb,:) = mean(slopes(ii_vb,:),1);
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
            % m_slopes = [nb_vb,nb_subjects]
            % p_slopes = [nb_subjects,nb_vb]
            p_slopes = permute(m_slopes,[2,1]);
            tools_dotplot(p_slopes,{''},{'o'},{[0,0,0]});
        end

        % do the analysis over two variables at the same time
        % 1 = setsize
        % 2 = variance
        % 3 = constraint
        % 4 = duration
        vbs = [1,4];
        if length(vbs)~=2
            return
        end
        u_vb1 = unique(psych(:,vbs(1)));
        u_vb2 = unique(psych(:,vbs(2)));
        nb_vb1 = length(u_vb1);
        nb_vb2 = length(u_vb2);
        m_slopes = nan(nb_vb1,nb_vb2,nb_subjects);
        for i_vb1 = 1:nb_vb1
            for i_vb2 = 1:nb_vb2
                i_psych = ((psych(:,vbs(1))==u_vb1(i_vb1)) & (psych(:,vbs(2))==u_vb2(i_vb2)));
                if size(slopes(i_psych,:),1)==1
                    m_slopes(i_vb1,i_vb2,:) = slopes(i_psych,:);
                else
                    m_slopes(i_vb1,i_vb2,:) = mean(slopes(i_psych,:),1);
                end
            end
        end

        % color
        colz = {};
        for i_vb1 = 1:nb_vb1
            colz{i_vb1} = ((i_vb1-1)/nb_vb1)*ones(1,3);
        end

        % figure
        figure('color',[1 1 1]);

        % x-axis
        x_axis = cell(1,nb_vb2);
        for i_xaxis = 1:nb_vb2
            x_axis{i_xaxis} = num2str(u_vb2(i_xaxis));
        end
        set(gca,'xtick',1:nb_vb2);
        set(gca,'xticklabel',x_axis);

        % plot values
        % m_slopes  = [nb_vb1,nb_v2,nb_subjects]
        % pm_slopes = [nb_subjects,nb_vb1,nb_v2]
        % so
        %{
            x = u_v2
            nb_v1 curves
            std driven by nb_subjects
        %}
        pm_slopes = permute(m_slopes,[3,1,2]);
        tools_dotplot(pm_slopes,{''},{'o','o','o','o'},colz);

        % select the legend for variables
        pos_legend = {' items',' variance',' contrained',' seconds'};
        c_legend = {};
        for i_vb1 = 1:nb_vb1
            c_legend{i_vb1} = [num2str(u_vb1(i_vb1)),pos_legend{vbs(1)}];
        end
        legend(c_legend,'Location','SouthEast')
        % show legend
        legend show;
    end

end