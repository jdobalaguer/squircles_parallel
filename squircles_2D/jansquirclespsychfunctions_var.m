close all

msrespcat = alldata.msrespcat;

dog=permute(msrespcat,[3 1 2]);
dog_ss = [];

% dog_ss(:,1,:)=mean(dog(:,1:2,:),2);
% dog_ss(:,2,:)=mean(dog(:,3:4,:),2);
% dog_ss(:,3,:)=mean(dog(:,5:6,:),2);

dog_ss(:,1,:)=mean(dog(:,[1 3 5],:),2);  % low variance
dog_ss(:,2,:)=mean(dog(:,[2 4 6],:),2);  % high variance

%
figure('color',[1 1 1]);
paramz=[0 1 4.5 0.5];
param=statset('MaxIter',1000);

colz={[0 0 0],[.4 .4 .4]};
hold on;
for j=1:size(dog_ss,2)

f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));  % 2=lower bound 2=upper bound, 3= inflection point 4=hills' slope 
[p r jj] = nlinfit(1:size(dog_ss,3),squeeze(mean(dog_ss(:,j,:)))',f,paramz,param);
hold on;

% for s=1:10;
% [pp(s,j,:) rr jj] = nlinfit(1:size(dog_ss,3),squeeze(dog_ss(s,j,:))',f,paramz,param);  %%1=lower bound 2=upper bound, 3= inflection point 4=hills' slope
% line(1:size(dog_ss,3),f(pp(s,j,:),1:size(dog_ss,3)),'color',colz{j},'linewidth',3);
% end


line(1:size(dog_ss,3),f(p,1:size(dog_ss,3)),'color',colz{j},'linewidth',3);

end

set(gca,'FontSize',16);
xlim([0.5 9.5]);
legend({'low variance','high variance'},'Location','SouthEast')
dotplot(dog_ss,{''},{'o','o','o'},colz);
set(gca,'xtick',1:9);
set(gca,'xticklabel',{'-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2'});
