function gammacor = tools_RGBcor(raw,lumin)
    % raw   = RGB in 1x3, RGB values defined by design
    % lumin = luminance we want to have between 0-1
    
    if size(raw,1)==1 && size(raw,2)==1
        raw=ones(1,3)*raw;
    end
    RGB = [0.300 0.590 0.110];
    
    % cor are the values observed by the subject so that the true luminance is "lumin"
    cor = raw./(RGB*raw')*lumin; % now cor is in 0-1 range
    gammacor = cor.^(1/2.2);
    
    % now gammacor is values we give to the graphics card, with a gamma correction, gamma=2.2
    gammacor = gammacor*255;% put the gammacor values back in the 0-255 range
end