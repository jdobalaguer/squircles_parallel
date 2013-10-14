function [zevalues] = tools_mypseudorandrange(array_mean,array_sd,nitms,constraint,critmean,critsd,range,nsamples)
% critmean is the maximal distance between outcome mean and expected mean
%       (0 = outcome very close to expected value)
% critsd is the maximal distance between outcome variance and expected variance
%       (0 = outcome very close to expected value)
% range is the range, of course

    % check for values
    if array_mean<=range(1) || array_mean>=range(2)
        error(['mypseudorange: error: array_mean(',num2str(array_mean),') out of range(',num2str(range),')\n']);
    end

    if array_sd==0
        zevalues = array_mean*ones(nsamples,nitms);
    elseif nitms==1
        zevalues = array_mean*ones(nsamples,nitms);
    else
        pass = 0;
        zevalues = zeros(0,nitms);
        counter = 0;
        while pass < nsamples
            counter = counter+1;
            if counter > 1000
                error(['mypseudorange: exit: 1000 loops(',num2str(array_mean),...
                                                      ',',num2str(array_sd),...
                                                      ',',num2str(nitms),...
                                                      ',',num2str(constraint),...
                                                      ',',num2str(critmean),...
                                                      ',',num2str(critsd),...
                                                      ',',num2str(range),...
                                                      ',',num2str(nsamples),')\n']);
            end
            values   = array_mean + array_sd*randn(nitms,10000); % generate 10000 samples
            out      = any( (values<range(1)) + (values>range(2)) )>0;
            if any(~out) % take samples that are not out of range
                values   = values(:,~out);
                distsd   = abs(std(values) - array_sd); % distance between the sample statistics and expected statistics
                distmean = abs(mean(values) - array_mean);
                zecrit   = distsd<critsd & distmean<critmean; % criterion
                if any(zecrit) || ~constraint % add the samples that pass the criterion
                    pass     = pass + sum(zecrit);
                    zevalues = [zevalues; values(:,zecrit)'];
                    if size(zevalues,1) > nsamples
                        zevalues((nsamples+1):end,:) = [];
                    end
                end
                % save for distribution statisticts
                %     n = n+1; zestds = [zestds std(zevalues)]; zemeans = [zemeans,mean(zevalues)];
            end
        end
    end

return





