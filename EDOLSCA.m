function[ge,pop]=EDOLSCA(pop,Gm_o,D,Np,lb,ub,fobj)
disp('                    EDOLSCA                    ')
%Elite Dynamic Opposite Learning Sin Cosin Algo
%Thierry Hu, 17-06-2020
elite_num=round(0.5*Np);
Lowerbound=lb;
Upperbound=ub;
Jr=0;
w=0;
Gm=Gm_o;
G=1;
%St=pop;
P=zeros(1,D);%Destination_position,最优值时的粒子坐标
Destination_position=zeros(1,D);
Destination_fitness=inf;
Convergence_curve=zeros(1,Gm);
Objective_values = zeros(1,Np);
%  for j=1:Np
%      Fit_Tr(j)=fobj(St(j,:)',func_num);%所有个体的适度值(1*Np)
%      if j==1
%          P=St(j,:);
%          Best=Fit_Tr(j);
%      elseif Fit_Tr(j)<P
%          P=St(j,:);
%          Best=Fit_Tr(j);
%      end
%  All_objective_values(j)=Fit_Tr(j);
for i=1:Np
    Objective_values(1,i)=fobj(pop(i,:),D);%所有个体的适度值(1*Np)
    if i==1
        Destination_position=pop(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=pop(i,:);
        Destination_fitness=Objective_values(1,i);%每代最优
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end
 


%SCA Phase%
while G<=Gm
    
 %Elite Phase 1
[elite_value,e_location]=min(Objective_values);%找到种群中最小的位置
elite=zeros(1,D);
elite=pop(e_location,:);%最小的种群

    a = 2;
    r1=a-G*((a)/Gm);
    
    for i=1:Np % in i-th solution
        for j=1:D % in j-th dimension
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();

            if r4<0.5
                % Eq. (3.1)
                pop_new(i,j)= pop(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-pop(i,j)));
%pop_new(i,j)= Checkbound(pop_new,Lowerbound,Upperbound,Np,D,G);
            else
                % Eq. (3.2)
                pop_new(i,j)= pop(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-pop(i,j)));
%pop_new(i,j)= Checkbound(pop_new,Lowerbound,Upperbound,Np,D,G);
            end
        end
    end
    %CheckBound
    for i=1:D
        for k=1:Np
            if pop_new(k,i)<lb(i)
                pop_new(k,i)=rand*(ub(i)-lb(i))+lb(i);%*(1+G/10);
            else if pop_new(k,i)>ub(i)
                pop_new(k,i)=ub(i)-rand*(ub(i)-lb(i));%/(1+G/10);
                end
            end
        end
    end
            fitnew=fobj(pop_new(i,:),D);
            if fitnew<Objective_values(i)
                pop(i,:)=pop_new(i,:);
                Objective_values(i)=fitnew;
                if fitnew<Destination_fitness
                    Destination_fitness=fitnew;
                end
            end

        ge(G)=Destination_fitness;
    G=G+1;
               %***************Dynamic opposite phase**************%  
    if rand<Jr
       for i=1:Np
              for j=1:D
                  Upperbound1(j)=max(pop(:,j));
                 Lowerbound1(j)=min(pop(:,j)); 
                  op(i,j)=Upperbound1(j)+Lowerbound1(j)-pop(i,j);
              end
            popO(i,:)=pop(i,:)+w*rand*(rand*op(i,:)-pop(i,:));
            popO_new(i,:)=pop(i,:);
            popO_new(Np+i,:)=popO(i,:);
            
                for i=1:D
        for k=1:Np
            if popO_new(k,i)<lb(i)
                popO_new(k,i)=rand*(ub(i)-lb(i))+lb(i);%*(1+G/10);
            else if popO_new(k,i)>ub(i)
                popO_new(k,i)=ub(i)-rand*(ub(i)-lb(i));%/(1+G/10);
                end
            end
        end
    end
end
%           for j=1:2*Np
%                Objective_values0(j)=fobj(popO_new(j,:),D);
%                if Objective_values0(j)<Destination_fitness
%                Destination_fitness=Objective_values0(j);
%                end
%           end
%              [Value,Index]=sort(Objective_values0);
%           for i=1:Np
%                  pop(i,:)=popO_new(Index(i),:);
%           end            

            fitnew=fobj(popO_new(i,:),D);
            if fitnew<Objective_values(i)
                pop(i,:)=popO_new(i,:);
                Objective_values(i)=fitnew;
                if fitnew<Destination_fitness
                    Destination_fitness=fitnew;
                end
            end
        
        if G==1
            ge(1)=Destination_fitness;
        end
        if G>1
            if Destination_fitness<ge(G-1)
                ge(G)=Destination_fitness;
             else 
                ge(G)=ge(G-1);   
            end
        end
        G=G+1;        
    end
%Elite, choice and save the best one, delete the worst one 
%Obj_values：已知的适应度，1*Np
%pop:种群，Np*Dim
[elite_value_new,e_location_new]=min(Objective_values);%找到种群中最小的位置
elite_new=zeros(1,D);
elite_new=pop(e_location_new,:);%最小的种群

for l=1:elite_num
[worst_value,worst_location]=max(Objective_values);%找到最差的
pop(worst_location,:)=[];%删除最差的
Objective_values(worst_location)=[];
end

if elite_value_new<elite_value
    elite=elite_new;
end
for m=1:elite_num
Objective_values(Np+1-m)=elite_value_new;
pop(Np+1-m,:)=elite;
end

end
    




[m , n]=size(ge);
% if n==5001
%     ge(:,1)=[];
% else
%     if n==5002
%         ge(:,1:2)=[];
%     end
% 
% end
% if n==10001
%     ge(:,1)=[];
% else
%     if n==10002
%         ge(:,1:2)=[];
%     end
% 
% end
% 
% if n==20001
%     ge(:,1)=[];
% else
%     if n==20002
%         ge(:,1:2)=[];
%     end
% 
% end
% if n==30001
%     ge(:,1)=[];
% else
%     if n==30002
%         ge(:,1:2)=[];
%     end
% 
% end
% if n==60001
%     ge(:,1)=[];
% else
%     if n==60002
%         ge(:,1:2)=[];
%     end
% 
% end
% if n==100001
%     ge(:,1)=[];
% else
%     if n==100002
%         ge(:,1:2)=[];
%     end
% 
% end
if n==500001
    ge(:,1)=[];
else
    if n==500002
        ge(:,1:2)=[];
    end

end
end

