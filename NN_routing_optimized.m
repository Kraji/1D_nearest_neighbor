clear
close all
set(0,'DefaultFigureWindowStyle','docked')

% I present here an algorithm to sort a 1D cyclic nearest neighbor network.
% The idea is to gather pairs of nodes defined using only nearest neighbor
% permutations.

%% Choice of the parameters
% I define first the number of nodes our 1D cyclic NN network will have.

n=10; % Number of nodes in the network
position=1:n; % Labels of the physical nodes.

%% Definition of the example interaction
% Here I generate a random interaction graph called "interaction" that
% states that the i-th node should be gathered with the interaction(i)-th
% node of the network.

interaction=1:n;
for i=1:2:n-1
    interaction([i i+1])=interaction([i+1 i]);
end;

for i=1:2*n
    a=randi(n);
    b=randi(n);
    interaction([interaction(a) interaction(b)])=interaction([interaction(b) interaction(a)]);
    interaction([a b])=interaction([b a]);
end
temp=interaction;
display(interaction);

%% Distance definition
% Here I define a function called "linear_chordal" that defines the
% distance between two nodes i and j to be linear_chordal(abs(i-j)+1).


% This is a parameter that that can be tuned to optimize the algorithm.
pow2=4;

% This is the definition of the distance between two nodes in the network.
linear_chordal=pow2.^round((n/pi)*asin(sin(pi/n*(position-1)))); 

% The actual overall distance in the network is D2
diag_distance=linear_chordal(abs(interaction-position)+1)';
D2=sum(diag_distance);

% The number "m" defined below is the minimum distance one can achieve if
% all nodes are gathered. This will be the goal to reach when we look at
% optimal permutations.
m=pow2*round(n-sum(interaction==position))+sum(interaction==position); 

%% Permutation
% In this section we find the optimal permutation to sort the network. For
% this purpose, we apply all nearest neighbor transpositions and look how
% the overall distance evolves. We only keep at each iteration the
% transposition that lowers the most the overall distance. Then we create a
% tabu list stating that we cannot permute nodes that have been permuted
% recently. If the tabu list is full or if one cannot reduce the overall
% distance anymore we empty the tabu list.

permutation=zeros(n^2/2,3); % This matrix will tell us which transposition to do to sort the network, it also tells us the overall distance after each transposition
j=1;

result=0; %This will be declared to 1 when we reach an optimal permutation to break the while loop.
disp(['The overall distance is ' num2str(D2)]);
distance_vector=D2*ones(n,1); %In this vector we will store the overall distance obtained after applying each transposition.
tabu=1:n; % We introduce a tabu list.

while round(D2)>m && j<n^2/2 && result==0
    for i=tabu(tabu>0)
        if i<n
            interac_temp=interaction;
            
            interac_temp([interac_temp(i) interac_temp(i+1)])=interac_temp([interac_temp(i+1) interac_temp(i)]);
            interac_temp([i i+1])=interac_temp([i+1 i]);
            
            distance_vector(i)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
            if round(distance_vector(i))==m
                disp('Optimal permutation has been found');
                result=1;
                break;
            end;
        else
            interac_temp=interaction;
            
            interac_temp([interac_temp(1) interac_temp(n)])=interac_temp([interac_temp(n) interac_temp(1)]);
            interac_temp([1 n])=interac_temp([n 1]);
            
            distance_vector(n)=min(round(sum(linear_chordal(abs(interac_temp-position)+1))),round(D2));
            if round(distance_vector(n))==m
                disp('Optimal permutation has been found');
                result=1;
            end
        end
    end
    
    [mini,index]=min(distance_vector);
    
    tabu(index)=0;
    if index==1
        tabu(n)=0;
        tabu(2)=0;
    elseif index==n
        tabu(1)=0;
        tabu(n-1)=0;
    else
        tabu(index+1)=0;
        tabu(index-1)=0;
    end
    
    if sum(tabu)==0
        tabu=1:n;
    end
    
    if min(distance_vector)<D2
        if min(distance_vector)==max(distance_vector)
            index=randi(n);
        end
        if index==n
            permutation(j,:)=[n,1,mini];
            interaction([interaction(1) interaction(n)])=interaction([interaction(n) interaction(1)]);
            interaction([n 1])=interaction([1 n]);
        else
            permutation(j,:)=[index,index+1,mini];
            interaction([interaction(index) interaction(index+1)])=interaction([interaction(index+1) interaction(index)]);
            interaction([index index+1])=interaction([index+1 index]);
        end
        j=j+1;
        D2=mini;
    else
        tabu=1:n;
    end
end

permutation(j:n^2/2,:)=[];

%% Number of time step
% In this section, we have all the transpositions that we need to implement
% to sort the network and we order them to minimize the number of parallel
% operations we can do in one time step to minimize the overhead of the
% network.

tr=permutation(:,1:2)';
tr2=zeros(size(tr,2),n);
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

% This the list of permutations that one needs to do to sort the network. 
% One needs to read this matrix line per line. 
% On each line if you see for instance « 0 0 4 3 0 7 6 0 0 » it means that elements 3 and 4 have to be permuted, and elements 7 and 6 have to be permuted. 
% The number of lines gives you the number of time steps required.
display(tr2);
disp([ num2str(time_step) ' timesteps required']);

