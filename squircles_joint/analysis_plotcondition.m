function analysis_plotcondition(sdata)

    % load data
    load_data;

    % do the analysis over... (see 'psych' variable)
    % 1 = setsize
    % 2 = variance
    % 3 = constraint
    % 4 = duration
    vbs = [1,4];
    nb_vbs = length(vbs);

    % plot ...
    % 1 = performance
    % 2 = error
    % 3 = absolute error
    % 4 = response
    % 5 = reaction times
    toplots = [1,3];
    
    
    for toplot = toplots
        % figure
        figure('color',[1 1 1]);

        % change figure size
        f_xy = get(gcf,'position');
        set(gcf,'position',[f_xy(1),f_xy(2),nb_vbs*f_xy(3),f_xy(4)]);
        for i_vb = 1:nb_vbs
            vb = vbs(i_vb);

            % subplot
            subplot(1,nb_vbs,i_vb);

            switch toplot
                case 1
                    mvalues = mcor;
                case 2
                    mvalues = merr;
                case 3
                    mvalues = maerr;
                case 4
                    mvalues = mcat;
                case 5
                    mvalues = mrt;
                otherwise
                    error(['analysis_plots: toplot=',num2str(toplot),'is not valid.']);
            end
            tools_plotcondition(sdata,mvalues,vb,toplot);
        end
    end
end
