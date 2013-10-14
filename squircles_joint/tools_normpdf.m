function y = tools_normpdf(x,p)
    if length(p)~=3
        error('tools_normpdf: length(p)~=3');
    end
    
    y = p(1) + p(2)*exp(-(x.*x)/p(3));
end