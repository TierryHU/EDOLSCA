function[ge,pop]=DOLSCA(pop,Gm_o,D,Np,lb,ub,fobj)
disp('                    DOLSCA                    ')
Lowerbound=lb;
Upperbound=ub;
Jr=0.5;
w=12;
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
  
%                 [pop_new]=CheckboundSCA(pop_new,Lowerbound,Upperbound,Np,D,G);
%             fitnew=fobj(pop_new(i,:)',func_num);
%             if fitnew<Objective_values(i)
%                 pop(i,:)=pop_new(i,:);
%                 Objective_values(i)=fitnew;
%                 if fitnew<Destination_fitness
%                     Destination_fitness=fitnew;
%                 end
%             end
%         Flag4ub=pop_new(i,:)>ub;
%         Flag4lb=pop_new(i,:)<lb;
%         pop_new(i,:)=(pop_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%         fitnew=fobj(pop_new(i,:)',func_num);
%         if fitnew< Destination_fitness
%             Destination_position=pop_new(i,:);
%             Destination_fitness=Objective_values(1,i);
%     end
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
if n==30001
    ge(:,1)=[];
else
    if n==30002
        ge(:,1:2)=[];
    end

end
end

