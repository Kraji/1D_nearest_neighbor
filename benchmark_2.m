clear
close all
set(0,'DefaultFigureWindowStyle','docked')

%% Choice of the parameters
upper=200;
lower=20;
stepp=4;
stat_max=zeros(length(lower:stepp:upper),1);
stat_mean=zeros(length(lower:stepp:upper),1);
stat_std=zeros(length(lower:stepp:upper),1);
pow=20;
z=1;
%%
h = waitbar(0,'Quantum computing being done...');
for n=lower:stepp:upper
    waitbar(n / upper,h,strcat('Step ',num2str(n)))
    iteration=1000;
    time_step=ones(iteration,1);
    
    position=(1:n);
    linear_chordal=round((n/pi)^pow*asin(sin(pi/n*(position-1))).^pow);
    %% Algorithm
    parfor k=1:iteration
        interaction=1:n;
        for i=1:2:n-1
            interaction([i i+1])=interaction([i+1 i]);
        end;
        
        for i=1:2*n
            a=randi(n);
            b=randi(n);
            interaction([a b])=interaction([b a]);
            interaction([interaction(a) interaction(b)])=interaction([interaction(b) interaction(a)]);
        end;
        m=floor(n-sum(interaction==position)); %number of nodes to gather
        
        %% Vector formulation
        
        diag_distance=linear_chordal(abs(interaction-position)+1)';
        D2=sum(diag_distance);
        %% Permutation
        permutation=zeros(n^2/2,3);
        j=1;
        result=0;
        distance_vector=D2*ones(n,1);
        
        while round(D2)>m && j<n^2/2 && result==0
            for i=1:n-1
                interac_temp=interaction;
                interac_temp([i i+1])=interac_temp([i+1 i]);
                interac_temp([interac_temp(i) interac_temp(i+1)])=interac_temp([interac_temp(i+1) interac_temp(i)]);
                distance_vector(i)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
                if round(distance_vector(i))==m
                    result=1;
                    break;
                end;
            end;
            
            interac_temp=interaction;
            interac_temp([1 n])=interac_temp([n 1]);
            interac_temp([interac_temp(1) interac_temp(n)])=interac_temp([interac_temp(n) interac_temp(1)]);
            distance_vector(n)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
            if round(distance_vector(n))==m
                result=1;
            end;
            
            [mini,index]=min(distance_vector);
            if min(distance_vector)==max(distance_vector)
                index=randi(n);
            end;
            if index==n
                permutation(j,:)=[n,1,mini];
                interaction([n 1])=interaction([1 n]);
                interaction([interaction(1) interaction(n)])=interaction([interaction(n) interaction(1)]);
            else
                permutation(j,:)=[index,index+1,mini];
                interaction([index index+1])=interaction([index+1 index]);
                interaction([interaction(index) interaction(index+1)])=interaction([interaction(index+1) interaction(index)]);
            end;
            j=j+1;
            D2=mini;
        end;
        permutation(j:n^2/2,:)=[];
        
        %% Number of time step
        
        tr=permutation(:,1:2)';
        tr2=zeros(size(tr,2),n);
        
        for i=1:size(tr,2)
            for j=time_step(k):-1:1
                if tr2(j,tr(1,i))==0 && tr2(j,tr(2,i))==0 && j~=1
                    continue
                elseif j==1 && tr2(j,tr(1,i))==0 && tr2(j,tr(2,i))==0
                    tr2(j,tr(1,i))=tr(2,i);
                    tr2(j,tr(2,i))=tr(1,i);
                    break;
                elseif tr2(j,tr(1,i))~=0 || tr2(j,tr(2,i))~=0
                    tr2(j+1,tr(1,i))=tr(2,i);
                    tr2(j+1,tr(2,i))=tr(1,i);
                    time_step(k)=max(time_step(k),j+1);
                    break;
                end
            end
        end
        
        tr2(time_step(k)+1:size(tr,2),:)=[];
    end
    
    stat_max(z)=max(time_step);
    stat_mean(z)=mean(time_step);
    stat_std(z)=std(time_step);
    z=z+1;
end;
%% Analysis of the data
save benchmark.mat
fit1 = polyfit(lower:stepp:upper,stat_max',1);
fit2 = polyfit(lower:stepp:upper,stat_mean',1);

max_fit=fit1(1)*(lower:stepp:upper)+fit1(2);
mean_fit=fit2(1)*(lower:stepp:upper)+fit2(2);

figure();
plot(lower:stepp:upper,stat_max);
xlabel('Number of nodes'); ylabel('Maximum depth'); title('Maximum depth of 1D cyclic NN network');
hold on;
plot(lower:stepp:upper,max_fit,'r');
legend('Simulation of the maximum depth',['Linear fit y= ' num2str(fit1(1)) 'x +' num2str(fit1(2))]);
savefig(gcf,'max')
saveas(gcf,'max','pdf')

figure();
errorbar(lower:stepp:upper,stat_mean,stat_std);
xlabel('Number of nodes'); ylabel('Average depth with standard deviation'); title('Average depth of 1D cyclic NN network');
hold on;
plot(lower:stepp:upper,mean_fit,'r');
legend('Simulation of average depth & standard deviation',['Linear fit y= ' num2str(fit2(1)) 'x ' num2str(fit2(2))]);
savefig(gcf,'mean')
saveas(gcf,'mean','pdf')