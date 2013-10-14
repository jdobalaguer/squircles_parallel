% model for squircles parallel

% load data
clear all;close all;

% switch here
sdata=load_badgui(1); % 1=exp 1; 2 = exp 4; 3 = exp 2; 4= exp3;

% new variables
sdata.vb_absmean=abs(sdata.vb_mean);

% free parameters
noisemat=[0:0.1:5];  % perceptual noise
smat=[0.2:0.2:4];        % slope of sigmoid
capmat=[2 4 8];         % capacity

% derived parameters
setsizes=unique(sdata.vb_setsize(sdata.vb_setsize>1));
variances=unique(sdata.vb_var);
durations=unique(sdata.vb_timing);
meanz=unique(sdata.vb_absmean);
submat=unique(sdata.exp_sub);

% fixed parameters
st = 50; % sampling time, 50ms.
plotornot=1;
f1=figure('color',[1 1 1]);
colz={'r','g','b'};

% sigmoid function
f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4))); % probability functions

% find participant data
for m=1:length(meanz);
    for v=1:length(variances);
        for z=1:length(setsizes);
            for d=1:length(durations);
                for ss=1:length(submat);
                    sindx=find(sdata.vb_absmean==meanz(m) & sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d) & sdata.exp_sub==submat(ss));
                    realcor(m,z,v,d,ss)=mean(sdata.resp_cor(sindx));
                end
            end
        end
    end
end
mrealcor=squeeze(mean(realcor,5));


h=waitbar(0,'yo');

%%
counter=0;
fit=nan(length(noisemat),length(smat),length(capmat));
for n=1:length(noisemat);
    for s=1:length(smat);
        counter=counter+1;
        
        if ~plotornot;
            waitbar(counter./(length(noisemat)*length(smat)),h);
        end
        
        for m=1:length(meanz);
            for v=1:length(variances);
                for z=1:length(setsizes);
                    for d=1:length(durations)
                        
                        samples=(durations(d)*1000)./st;
                        
                        indx=find(sdata.vb_absmean==meanz(m) & sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d));
                        x=(-1.5+sdata.sample_S(1:setsizes(z),indx)')*4;;  % input data, normalised to mean zero
                        
                        paramz=[-1 2 0 smat(s)];
                        R=f(paramz,x);    % R is now a value between -1 and 1
                        R=R+randn(length(indx),setsizes(z))*noisemat(n);  % add noise

                        for cp=1:length(capmat); % up to capacity
                            clear dv;
                            for p=1:samples
                                sR=tools_shuffle(R,2); % shuffle
                                dv(:,p)=mean(sR(:,1:min(capmat(cp),setsizes(z))),2);
                            end
                            tdv=mean(dv,2);
                            cor{cp}(m,z,v,d)=mean(sign(tdv')==sign(sdata.sample_Scat(indx)));
                        end
                        
                    end
                end
            end
        end
        
        %% plotting
        if plotornot
            figure(f1);
            for cp=1:length(capmat);
                for z=1:length(setsizes);
                    subplot(1,length(setsizes),z);
                    cla;
                    for r=1:size(mrealcor,3);
                        plot(squeeze(mrealcor(:,z,r)),'o','markersize',10,'markerfacecolor',colz{r},'markeredgecolor',colz{r});
                        hold on;
                        plot(squeeze(cor{cp}(:,z,r))','-','linewidth',2,'color',colz{r});
                    end
                    
                    xlim([0.5 size(cor{cp},1)+0.5]);
                    ylim([0.4 1]);
                    title(['noise=',num2str(noisemat(n)),' slope=',num2str(smat(s)),' cp=',num2str(capmat(cp))]);
                    % drawnow;
                end
            end
        end
        for cp=1:length(capmat);
            fit(n,s,cp)=sum(abs(squelch(squeeze(cor{cp})-mrealcor)));
        end
        
    end
end


% plot exhaustively
viz(fit);

i=find4(fit==min(fit(:)));
n=i(1);
s=i(2);
cp=i(3);

for m=1:length(meanz);
    for v=1:length(variances);
        for z=1:length(setsizes);
            for d=1:length(durations)
                
                samples=(durations(d)*1000)./st;
                
                indx=find(sdata.vb_absmean==meanz(m) & sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d));
                x=(-1.5+sdata.sample_S(1:setsizes(z),indx)')*4;;  % input data, normalised to mean zero
                
                paramz=[-1 2 0 smat(s)];
                R=f(paramz,x);    % R is now a value between -1 and 1
                R=R+randn(length(indx),setsizes(z))*noisemat(n);  % add noise
               
                clear dv;
                for p=1:samples
                    sR=tools_shuffle(R,2); % shuffle
                    dv(:,p)=mean(sR(:,1:min(capmat(cp),setsizes(z))),2);
                end
                tdv=mean(dv,2);
                finalcor(m,z,v,d)=mean(sign(tdv')==sign(sdata.sample_Scat(indx)));
                
            end
        end
    end
end


%%
f2=figure('color',[1 1 1]);
for z=1:length(setsizes);
    subplot(1,length(setsizes),z);
    for r=1:size(mrealcor,3);
        plot(squeeze(mrealcor(:,z,r)),'o','markersize',10,'markerfacecolor',colz{r},'markeredgecolor',colz{r});
        hold on;
        plot(squeeze(finalcor(:,z,r))','-','linewidth',2,'color',colz{r});
    end
    xlim([0.5 size(finalcor,1)+0.5]);
    ylim([0.4 1]);
    if z==1;
        title(['noise=',num2str(noisemat(n)),' slope=',num2str(smat(s)),' cp=',num2str(capmat(cp))]);
    end
end


close (h);
