function analysis_errorbar(sdata)

    % load data
    load_data;

    % do the analysis over... (see 'psych' variable)
    % 1 = setsize
    % 2 = variance
    % 3 = constraint
    % 4 = duration
    vbs = [4,1];

    % plot ...
    % 1 = performance
    % 2 = error
    % 3 = absolute error
    % 3 = response
    % 4 = reaction times
    toplots = [1,3];

    % check vbs
    if length(vbs)~=2 || vbs(1)==vbs(2)
        error('analysis_errorbar: error. vbs needs to values.');
    end
    
    for toplot = toplots
        % select variable to plot
        fit = 0;
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
        tools_psycherrorbar(sdata,mvalues,vbs,fit);
    end
end
