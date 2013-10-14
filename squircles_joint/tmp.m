
clear all;

sdata = load_badgui(1);
load_numbers;

cor_S = nan(1,nb_setsizes);
cor_A = nan(1,nb_setsizes);

for i_setsize = 1:nb_setsizes
    
    setsize = u_setsize(i_setsize);
    ii_trials = (sdata.vb_setsize==setsize);
    
    x = sdata.sample_V(1:setsize,ii_trials);
    if setsize>1
        mx = mean(x);
    else
        mx = x;
    end
    sx = (mx > 0);
    
    % response S
    rx_S = (x(1,:) > 0);
    cor_S(i_setsize) = mean(rx_S == sx);
    
    % response A
    rx_A = (mx > 0);
    cor_A(i_setsize) = mean(rx_A == sx);
    
end

plot(u_setsize,cor_A,u_setsize,cor_S,'+-');
ylim([0,1]);