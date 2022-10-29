 close all
 clc
 clear 
 iter=1;
 dim=7;
 format long

 Gm_o=500;
 AN=6;
%  Best_score=zeros(AN,Gm_o,iter);
 Best_score=zeros(AN,Gm_o);
 SearchAgents_no=50;
 FES=Gm_o*SearchAgents_no;
 %f = @funcAHT2019;
 f = @funcNs2019;
 Mean=zeros(1,AN);
 Std=zeros(1,AN);
 run = 1;
for k=1:run
 for i=1: SearchAgents_no

 L1 = rand*0.9+0.1;
 L2 = rand*0.9+0.1;
 H = rand*(0.01-0.002)+0.002;
 t = rand*0.0001+0.0001;
 n = rand*900+100;
 l = rand*(0.01-0.001)+0.001;
% Nh = rand*9+1;
 Nh = rand*199+1;
 x = [L1,L2,H,t,n,l,Nh];
 
 pop(i,:) = x;

 end
 
Lowerbound = [0.1 0.1 0.002 0.0001   100 0.001 1];
% Upperbound = [1    1   0.01  0.0002 1000 0.01  200];%2
Upperbound = [1    1   0.01  0.0002 1000 0.01  10];%1
lb = Lowerbound;
ub = Upperbound;
  
% for func_num=[1]
%     for count=1:iter
% disp(count)
tic;
% 
%      [Best_score(1,:,count)]=wPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
%      [Best_score(2,:,count)]=cfPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
%      [Best_score(3,:,count)]=cfwPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
%      [Best_score(4,:,count)]=ETLBO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
%      [Best_score(5,:,count),popSCA]=SCA(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
%      [Best_score(6,:,count),popDOLSCA]=EDOLSCA(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(1,:)]=wPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(2,:)]=cfPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(3,:)]=cfwPSO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(4,:)]=ETLBO(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(5,:),popSCA]=SCA(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
     [Best_score(6,:),popDOLSCA]=EDOLSCA(pop,Gm_o,dim,SearchAgents_no,lb,ub,f);
toc;
TotalTime = toc;
disp(['You ve wasted  ' num2str(TotalTime) '  seconds of your life']);
    
    %Best_score=Best_score-func_num*100;
%     if Best_score<1e-08
%         Best_score=0;
%     end
% 
%     for i=1:AN
%         Mean(i)=mean(Best_score(i,:));     
%         for j=1:iter
%             Std(i,j)=std(Best_score(i,:));
%         end
%     end
% 
% ii=linspace(1,Gm_o,Gm_o);   

%     end



%  for k = 1:AN
%  A=Best_score(k,:,:);
%  A=squeeze(A);
%  MIN=min(A(Gm_o,:));
%  MAX=max(A(Gm_o,:));
%  Min(:,k)=MIN;
%  Max(:,k)=MAX;
%  A=mean(A,2);
%  mean_score(:,k)=A;
%  end
%   for k = 1:AN
%  B=Best_score(k,:,:);
%  B=squeeze(B);
%  B=B(Gm_o,:);
%  STd(:,k)=std(B);
%  end
 
% mean=mean_score(Gm_o,:);
% mean=mean';
% STd=STd';
% Min=Min';
% Max=Max';
best=Best_score';
figure1 = figure('OuterPosition',[1133 509.8 521.6 464.8]);
axes1 = axes('Parent',figure1,...
    'Position',[0.0967996705313203 0.126586344633663 0.88157894736842 0.814410480349345]);
hold(axes1,'on');

% plot(log10(mean_score(:,1)),'-v','color','#377245','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
% plot(log10(mean_score(:,2)),'-^','color','#6977ab','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
% plot(log10(mean_score(:,3)),'-o','color','#9699e9','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
% plot(log10(mean_score(:,4)),'-x','color','#484ca3','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
% plot(log10(mean_score(:,5)),'-<','color','#4f54bf','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
% plot(log10(mean_score(:,6)),'-h','color','#8e9950','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,1)),'-v','color','#377245','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,2)),'-^','color','#6977ab','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,3)),'-o','color','#9699e9','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,4)),'-x','color','#484ca3','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,5)),'-<','color','#4f54bf','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;
plot(log10(best(:,6)),'-h','color','#8e9950','linewidth',0.7,'MarkerIndices',1:Gm_o/25:Gm_o,'Markersize',6);hold on;

legend('PSO','cfPSO','cfwPSO','ETLBO','SCA','EDOLSCA');
xlabel('FES');
%title('CF2');
% ylabel('log10(Number of Entropy Production)');
ylabel('log10(Number of Heat Transfer Area)');
box(axes1,'on');
set(axes1,'XTickLabel',{'0','5','10','15','20','25'});
legend(axes1,'show');
annotation(figure1,'textbox',...
    [0.90694006309148 0.0276595744680851 0.0767823343848568 0.0531914893617021],...
    'String',{'¡Á10^{4}'},...
    'LineStyle','none',...
    'FitBoxToText','off');
Filename = 'NS';
save(Filename)

end