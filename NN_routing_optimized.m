clear
close all
set(0,'DefaultFigureWindowStyle','docked')

%% Choice of the parameters
n=300; % number of nodes
pow=20;
position=1:n;

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

linear_chordal=round((n/pi)^pow*asin(sin(pi/n*(position-1))).^pow);
diag_distance=linear_chordal(abs(interaction-position)+1)';
D2=sum(diag_distance);

%% Permutation
permutation=zeros(n^2/2,3);
j=1;
result=0;
disp(['The quadratic node distance is ' num2str(D2)]);
distance_vector=D2*ones(n,1);
tic
while round(D2)>m && j<n^2/2 && result==0
    for i=1:n-1
        interac_temp=interaction;
        interac_temp([i i+1])=interac_temp([i+1 i]);
        interac_temp([interac_temp(i) interac_temp(i+1)])=interac_temp([interac_temp(i+1) interac_temp(i)]);
        distance_vector(i)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
        if round(distance_vector(i))==m
            disp('Optimal permutation has been found');
            result=1;
            break;
        end;
    end;
    
    interac_temp=interaction;
    interac_temp([1 n])=interac_temp([n 1]);
    interac_temp([interac_temp(1) interac_temp(n)])=interac_temp([interac_temp(n) interac_temp(1)]);
    distance_vector(n)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
    if round(distance_vector(n))==m
        disp('Optimal permutation has been found');
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
time=toc
permutation(j:n^2/2,:)=[];
%display(permutation);
%display(interaction);

%% Number of time step

tr=permutation(:,1:2)';
tr2=zeros(size(tr,2),n);
%disp(tr);
time_step=1;

for i=1:size(tr,2)
    for j=time_step:-1:1
        if tr2(j,tr(1,i))==0 && tr2(j,tr(2,i))==0 && j~=1
            continue
        elseif j==1 && tr2(j,tr(1,i))==0 && tr2(j,tr(2,i))==0
            tr2(j,tr(1,i))=tr(2,i);
            tr2(j,tr(2,i))=tr(1,i);
            break;
        elseif tr2(j,tr(1,i))~=0 || tr2(j,tr(2,i))~=0
            tr2(j+1,tr(1,i))=tr(2,i);
            tr2(j+1,tr(2,i))=tr(1,i);
            time_step=max(time_step,j+1);
            break;
        end
    end
end

tr2(time_step+1:size(tr,2),:)=[];
verification=sum(tr2(:))==sum(tr(:));

%disp(tr2);
disp(verification);
display(time_step);

