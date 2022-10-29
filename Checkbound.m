function [St]=Checkbound(St,Stlow_bound,Stup_bound,Np,Dimension,G)

    for i=1:Dimension
        for k=1:Np
            if St(k,i)<Stlow_bound(i)
                St(k,i)=rand*(Stup_bound(i)-Stlow_bound(i))+Stlow_bound(i);%*(1+G/10);
            else if St(k,i)>Stup_bound(i)
                St(k,i)=Stup_bound(i)-rand*(Stup_bound(i)-Stlow_bound(i));%/(1+G/10);
                end
            end
        end
    end
end