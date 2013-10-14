% model for squircles parallel

% load data
clear all;close all;

sdata=load_badgui(1); % 1=exp 1; 2 = exp 4; 3 = exp 2; 4= exp3;

% free parameters
noisemat=[0.2:0.2:2];  % perceptual noise
smat=[0.2:0.1:1];        % slope of sigmoid
capmat=[1:12];         % capacity

% derived parameters
setsizes=unique(sdata.vb_setsize(sdata.vb_setsize>1));
variances=unique(sdata.vb_var);
durations=unique(sdata.vb_timing);
submat=unique(sdata.exp_sub);

% fixed parameters
st = 50; % sampling time, 50ms.
plotornot=1;
f1=figure('color',[1 1 1]);
colz={'r','g','b'};

% sigmoid function
f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4))); % probability functions

% find participant data
for v=1:length(variances);
    for z=1:length(setsizes);
        for d=1:length(durations);
            for ss=1:length(submat);
                sindx=find(sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d) & sdata.exp_sub==submat(ss));
                realcor(z,v,d,ss)=mean(sdata.resp_cor(sindx));
            end
        end
    end
end
mrealcor=mean(realcor,4);

%cor=nan(length(variances),length(setsizes),length(durations));
fit=nan(length(noisemat),length(smat),length(capmat));
for n=1:length(noisemat);
    for s=1:length(smat);
        
        % waitbar here
        
        for v=1:length(variances);
            for z=1:length(setsizes);
                for d=1:length(durations)
                    
                    samples=(durations(d)*1000)./st;
                    
                    indx=find(sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d));
                    x=(-1.5+sdata.sample_S(1:setsizes(z),indx)')*4;;  % input data, normalised to mean zero
                    x=x+randn(length(indx),setsizes(z))*noisemat(n);  % add noise
                    
                    paramz=[-1 2 0 smat(s)];
                    R=f(paramz,x);    % R is now a value between -1 and 1
                    
                    for cp=1:length(capmat); % up to capacity
                        clear dv;
                        for p=1:samples
                            sR=tools_shuffle(R,2); % shuffle
                            dv(:,p)=mean(sR(:,1:min(capmat(cp),setsizes(z))),2);
                        end
                        tdv=mean(dv,2);
                        cor{cp}(z,v,d)=mean(sign(tdv')==sign(sdata.sample_Scat(indx)));
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
                  for r=1:size(mrealcor,2);
                        plot(squeeze(mrealcor(z,r,:)),'o','markersize',10,'markerfacecolor',colz{r},'markeredgecolor',colz{r});
                        hold on;
                        plot(squeeze(cor{cp}(z,r,:))','d-','linewidth',2,'color',colz{r});
                    end
                    xlim([0.5 length(durations)+0.5]);
                    ylim([0.5 1]);
                    title(['noise=',num2str(noisemat(n)),' slope=',num2str(smat(s)),' cp=',num2str(capmat(cp))]);
                    % drawnow;
                end
            end
        end
        for cp=1:length(capmat);
            fit(n,s,cp)=sum(abs(squelch(cor{cp}-mrealcor)));
        end
        
    end
end


% plot exhaustively
viz(fit);

i=find4(fit==min(fit(:)));
n=i(1);
s=i(2);
cp=i(3);

for v=1:length(variances);
    for z=1:length(setsizes);
        for d=1:length(durations)
            
            samples=(durations(d)*1000)./st;
            
            indx=find(sdata.vb_var==variances(v) & sdata.vb_setsize==setsizes(z) & sdata.vb_timing==durations(d));
            x=(-1.5+sdata.sample_S(1:setsizes(z),indx)')*4;;  % input data, normalised to mean zero
            x=x+randn(length(indx),setsizes(z))*noisemat(n);  % add noise
            
            paramz=[-1 2 0 smat(s)];
            R=f(paramz,x);    % R is now a value between -1 and 1
            
            clear dv;
            for p=1:samples
                sR=tools_shuffle(R,2); % shuffle
                dv(:,p)=mean(sR(:,1:min(capmat(cp),setsizes(z))),2);
            end
            tdv=mean(dv,2);
            finalcor(z,v,d)=mean(sign(tdv')==sign(sdata.sample_Scat(indx)));
            
        end
    end
end


%%
f2=figure('color',[1 1 1]);
for z=1:length(setsizes);
    subplot(1,length(setsizes),z);
    for r=1:size(mrealcor,2);
        plot(squeeze(mrealcor(z,r,:)),'o','markersize',10,'markerfacecolor',colz{r},'markeredgecolor',colz{r});
        hold on;
        plot(squeeze(finalcor(z,r,:))','d-','linewidth',2,'color',colz{r});
    end
    xlim([0.5 size(finalcor,3)+0.5]);
    ylim([0.5 1]);
    if z==1;
        title(['noise=',num2str(noisemat(n)),' slope=',num2str(smat(s)),' cp=',num2str(capmat(cp))]);
    end
end



