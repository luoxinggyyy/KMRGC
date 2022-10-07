function [output] = mrgcdemo(input,k)  % input is the input data, k is the nearest neighbor coefficient, and output is the clustering result and saved in the form of tuple.
data = input;
[m,~] = size(data); 
pattern = data;
%% Calculate the Euclidean distance (calculate the distance of the points)
dist1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        dist1(i,j) = norm(data(i,:)-data(j,:));% Distance formula, (dist (i, j) represents the distance between the i-th row and the j-th row)
    end
end
%% Sorting the distance
k1 =k; 
[p,I] = sort(dist1);%p is the sorted array, I contains the position of the array before sorting
%% Find the coefficient of attraction and store the ones in the neighborhood in the structure (the attraction of the point with other points is counted)
yy = I(2:k1,1:m);% Adjacent sets
rr = struct([]);
for t = 1:m
    [r,o] = find(yy==t); 
    rr(t).a = t;   %rr is a structure, and the first column stores the column where t is located, which is the position of the point
    rr(t).b = r; %r represents the proximity of point t to other points
    rr(t).c = o; %o represents the location of other points in the neighboring domain with t as the center point
end
for t = 1:m
    rr(t).e = sum(exp(-(rr(t).b)/10));% represents the sum of S(n) in the article
end
for t = 1:m
    rr(t).f = (rr(t).e - min([rr.e]))/(max([rr.e]) - min([rr.e]));  %The NI(x) in the article, i.e., into the neighborhood index and after normalization, is stored in s.f
end

    
%% Calculating KRI
KRI = zeros(1,m);
for t = 1:m
    for j = 1:k1-1
        if rr(yy(j,t)).f > rr(t).f
            KRI(t) = rr(t).f * j * p(j+1,t);
            break
        end

    end
end
[ykri,~] = sort(KRI,'descend');  % Descending order
figure
bar(ykri(1:30));  %Draw bar chart
hold on;
xlabel('Cluster number'); ylabel('KRI');

%% Find cluster centroids
hh = (1:m);
ge =vertcat(hh, yy); %Attractive set
tt = struct([]);
for i =1:m
    [~,temp] = min([rr(ge(:,i)).f]);%Boundary Points
    tt(i).a = i;
    tt(i).b = ge(temp,i);
    [~,temp1] = max([rr(ge(:,i)).f]);% Center Point
    tt(i).c = ge(temp1,i);
    tt(i).d = ge(:,i);
    tt(i).e = max(KRI([tt(i).d]));
end

kkk = [];%Array of boundary points
attr = zeros(k1,k1);
for i = 1:m
    w = 0;
    for j = 1:k1
        for z = 1:k1
            if ~isempty(intersect(tt(i).d(j),tt(tt(i).d(z)).d))
                attr(j,z) = rr(tt(i).d(j)).f - rr(tt(i).d(z)).f; %The attraction equation
            else
                attr(j,z) = -rr(tt(i).d(z)).f;
            end
        end
        if length(find(attr(j,:) <= 0)) == k1
            w = w+1;
            kkk(i,w) = tt(i).d(j); %kkk is the set of boundary points
            
        end
    end
    tt(i).f = kkk(i,:); %Store the boundary points in column f of the structure tt 
    tt(i).f(tt(i).f == 0) = [];%Remove the zeros from the set of boundary points
end
% tt.f is the boundary point, tt.c is the cluster centroid, and tt.d is the neighborhood set of the column


%% Calculate channel 
vvv = struct([]); %Used to place channels
for i = 1:m
    ww = 0;
    nn = 0;
    for j = 1:m
        if intersect(tt(i).f,tt(j).d) > 0 %Determine whether the boundary point intersects with other neighboring sets
            ww = ww+1;
            vvv(i).a = i;
            vvv(i).b(ww) = j; %Store the columns corresponding to clusters i and j, and if there are more than one, put them together
            vvv(i).c(ww) = {intersect([tt(i).f],[tt(j).d])}; %Store the channels in a tuple way to distinguish the columns corresponding to each array, the intersection of the boundary points corresponding to the columns with the other neighboring and
        end
        if intersect(tt(i).d,tt(j).f) > 0
            nn = nn+1;
            vvv(i).e(nn) = j;
            vvv(i).d(nn) = {intersect(tt(i).d,tt(j).f)};% The intersection of the set of neighbors corresponding to this column with the other boundary points
        end
    end
    vvv(i).f = union(vvv(i).b,vvv(i).e);

    vvv(i).f(vvv(i).f == i) = [];
end



for i = 1:m
    d = tt(i).f; %The boundary point of the column
    f = vvv(i).f; % The column that intersects the column
    ppp = zeros(1,length(f));
    for j = 1:length(f)
        g = tt(f(j)).f;% The boundary point of the column that intersects the column
        ddd = zeros(1,length(d));
        for z = 1:length(d)
            disg =zeros(1,length(g));
            for zz = 1:length(g)
                disg(zz) =  norm(data(d(z),:)-data(g(zz),:));

            end
            a = (disg~=0);
            b = sum(a(:));
            if b > 0
                disg(disg==0) = [];%Removing the case where the channel is 0
                [~,temp2] = sort(disg);
                ddd(z) = disg(temp2(1));
            else
                ddd(z) = nan;
            end
        end
        ppp(1,j) = min(ddd);
    end
    vvv(i).h = ppp; %ppp is the minimum channel distance of the column corresponding to the other columns
end

for i = 1:m
    tt(i).g = vvv(i).f; %The set of intersections of the column with the boundary points of other columns 
    tt(i).h = vvv(i).h; %Corresponding channel distance 
end
        
% Sorting of channel distances
for i = 1:m
    xxx = zeros(1,length(tt(i).g));
    [~,xxx] = sort(tt(i).h,'MissingPlacement','last');
    tt(i).g = tt(i).g(xxx);
    tt(i).h = sort(tt(i).h,'MissingPlacement','last');
end
% Elimination of NaN
for i = 1:m
    for j = 1:length(tt(i).g)
        if ~isnan(tt(i).h(j))
            tt(i).i(j) = tt(i).g(j);
            tt(i).j(j) = tt(i).h(j);
        else
            break
        end
    end
end


for i = 1:m
    tt(i).k  = cat(2,tt(i).a,tt(i).i);
    tt(i).l = length(tt(i).k);
end

sss = struct([]); %Used to place channels
for i =1:m
    sss(i).a = tt(i).j;
    sss(i).b = tt(i).i;
    sss(i).c = {};
    for j = 1:length(sss(i).b)
        sss(i).c(j) = {[i,sss(i).b(j)]};
    end
end

ge = [sss.a];%Put is the clustering of the boundary points of the two clusters to be fused
go = [sss.c];% The two clusters placed are the ones to be fused
[ge1,ge2] = sort(ge);% Sorting of channel distances
go1 = go(ge2); % The corresponding channels are also listed in ascending order
ger = ge1';
go2 =go1';
biu = cell2mat(go2);
get = cat(2,ger,biu); %Sorted channel distance and channel

dk = unique(get(:,1)); 
%%  Reverse sorting using NI at channel boundaries
gett = get; 
for i = 1:length(gett)
    gett(i,4) = max(rr(tt(gett(i,2)).b).f,rr(tt(gett(i,3)).b).f);
end
[~,gev2] = sort(gett(:,4),"descend");
getx = gett(gev2,:);
dkk = unique(getx(:,4));
dkk = sort(dkk,"descend");%Here the getx variables, the first column is the channel distance, the second and third columns are the clusters to be fused, and the third column is the minimum value of the NI at the boundary point of the cluster where these two clusters are to be fused, because the NI needs to be fused sequentially from largest to smallest, so directional ordering is used.
%% 
Q = {[]};
for i = 1:length(dkk)
    [row,~] = find(getx==dkk(i));
    U = getx(row,2);
    Q(i) = {unique(U)'};
end
%% The neighborhood set corresponding to the fusion set
B = {};
for i = 1:length(Q)
    B(i) = {unique([tt(cell2mat(Q(i))).d])'};
end
%% Integration
P =Q;
for i = 1:length(dkk)-1
    M = P(1:i);
    if intersect(cell2mat(P(i+1)),cell2mat(P(1:i))) > 0
        for j = 1:i
            if intersect(cell2mat(P(i+1)),cell2mat(M(j))) > 0
                P(j) = {union(cell2mat(P(i+1)),cell2mat(M(j)))};
                P{i+1} =[];
                break
            end
        end
    end
    W = P(1:i);

    if length(unique(cell2mat(P(1:i)))) >= m 
        break
    end
end

K = {};
j = 1;
for i = 1:length(W)
    if ~isempty(cell2mat(W(i)))
        K(j) = W(i); %  K is the final clustering effect
        j = j+1;
    end
end
output = K;