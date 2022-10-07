function [output] = kmrgcdemo(input,k)  % input is the input data, k is the nearest neighbor coefficient, and output is the clustering result and saved in the form of tuple.
data = input;
[m,n] = size(data); 
pattern = data; 
%% Calculate the Euclidean distance (calculate the distance of the points)
dist1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        dist1(i,j) = norm(data(i,:)-data(j,:));% Distance formula, (dist (i, j) represents the distance between the i-th row and the j-th row)
    end
end
%% Sorting the distance
k1 = k;
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
    rr(t).e = sum(exp(-(rr(t).b)/10));% represents the sum of S(n)
end
for t = 1:m
    rr(t).f = (rr(t).e - min([rr.e]))/(max([rr.e]) - min([rr.e]));  %The NI(x) in the article, i.e., into the neighborhood index and after normalization, is stored in s.f
end


%% Calculating KRI
KRI = zeros(1,m);
for t = 1:m
    for j = 1:k1-1
        if rr(yy(j,t)).f >= rr(t).f
            KRI(t) = rr(t).f * j * p(j+1,t); %KRI is the core index, here is the set of neighbors of column t has NI greater than it, then execute this
            break
        end
    end
    if max([rr(yy(:,t)).f]) < rr(t).f
        KRI(t) = rr(t).f*0.0001; % 这里是为了避免KRI为0的情况，后续都为0，采用的是NI乘与一个很小的数，使得不对前面的KRI产生影响同时有大小差异。
    end
end
figure
ykri = sort(KRI,'descend');  % Descending order
bar(ykri(1:30));  %Draw bar chart
hold on;
xlabel('Cluster number'); ylabel('KRI');
%% Find cluster centroids
hh = (1:m);
ge =vertcat(hh, yy); %Attractive set
tt = struct([]);
for i =1:m
    [~,temp] = min([rr(ge(:,i)).f]);% Boundary Points
    tt(i).a = i;
    tt(i).b = ge(temp,i);
    [~,temp1] = max([rr(ge(:,i)).f]);% Center Point
    tt(i).c = ge(temp1,i);
    tt(i).d = ge(:,i);
    tt(i).e = KRI(tt(i).c);
    tt(i).m = rr(i).f;
end

%% Computational clusters
a = [1:m];
u = [];
hhh = {};
for i = 1:m
    uu = setdiff(a,u);
    [~,zz] = max ([tt(uu).m]);
    Ctmax = tt(uu(zz)).c;
    ct = zeros(1,length(uu));
    for j = 1:length(uu)
        ct(j) = norm(data(Ctmax,:)-data([tt(uu(j)).c],:));
    end
    [~,CT] = sort(ct);
    YY = [];% Empty array, used for initialization of each cluster
    for z = 1:m
        YY = union(YY,[tt(uu(CT(z))).d]); %denotes that u is merged with the neighboring set of the corresponding column and then assigned to u, playing the role of annexation and accumulation
        if length(u) + z < m
            if length(intersect(YY,[tt(uu(CT(z+1))).d])) == 0 %If there is no intersection, exit the loop. (meaning that the fusion of clusters is complete and the next larger cluster is executed for fusion)
                break
            end
        else
            break
        end
    end
    hhh(i) = {YY'};
    u = cell2mat(hhh);
    u = unique(u);
    if length(u) >= m % The loop ends when the sum of the columns contained in all clusters is the same as the total number of points
        break
    end
end
output = hhh;