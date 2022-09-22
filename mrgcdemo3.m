function [output] = mrgcdemo3(input,k)  % input为输入数据，k为最邻近系数，output为聚类结果，并以元胞的形式保存。
data = input;
[m,~] = size(data); % m = 300,n = 2
pattern = data; %数据制作完成
%% 计算欧式距离(计算所有点的距离)
dist1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        dist1(i,j) = norm(data(i,:)-data(j,:));% 距离公式,（dist（i，j）代表第i行与第j行的距离）
    end
end
%% 对距离进行排序
k1 =k; 
[p,I] = sort(dist1);%p为排序后的数组，I包含了排序前数组的位置
%% 求吸引系数，将邻近区域中的储存在结构体中  （算的是点与其他点的吸引）
yy = I(2:k1,1:m);% 邻近集
rr = struct([]);
for t = 1:m
    [r,o] = find(yy==t); 
    rr(t).a = t;   %rr为一个结构体，第一列储存的是t所在的列，就是点的位置，
    rr(t).b = r; %r代表t点与其他点的相邻情况
    rr(t).c = o; %o代表t为中心点组成的邻近域中其他点所在位置
end
for t = 1:m
    rr(t).e = sum(exp(-(rr(t).b)/10));% 代表了文章中S(n)的和
end
for t = 1:m
    rr(t).f = (rr(t).e - min([rr.e]))/(max([rr.e]) - min([rr.e]));  %文章中的NI（x），即进邻指数，且归一化之后的，储存在了s.f中
end

    
%% 计算KRI
KRI = zeros(1,m);
for t = 1:m
    for j = 1:k1-1
        if rr(yy(j,t)).f > rr(t).f
            KRI(t) = rr(t).f * j * p(j+1,t);
            break
        end

    end
end
[ykri,bbb] = sort(KRI,'descend');  % 降序排列
figure
bar(ykri(1:30));  %画柱状图
hold on;
xlabel('Cluster number'); ylabel('KRI');

%% 找到聚类中心点
hh = (1:m);
ge =vertcat(hh, yy); %吸引集
tt = struct([]);
for i =1:m
    [~,temp] = min([rr(ge(:,i)).f]);% 边界点
    tt(i).a = i;
    tt(i).b = ge(temp,i);
    [~,temp1] = max([rr(ge(:,i)).f]);% 中心点
    tt(i).c = ge(temp1,i);
    tt(i).d = ge(:,i);
    tt(i).e = max(KRI([tt(i).d]));
end

kkk = [];%放边界点的
attr = zeros(k1,k1);
for i = 1:m
    w = 0;
    for j = 1:k1
        for z = 1:k1
            if ~isempty(intersect(tt(i).d(j),tt(tt(i).d(z)).d))
                attr(j,z) = rr(tt(i).d(j)).f - rr(tt(i).d(z)).f; %吸引方程
            else
                attr(j,z) = -rr(tt(i).d(z)).f;
            end
        end
        if length(find(attr(j,:) <= 0)) == k1
            w = w+1;
            kkk(i,w) = tt(i).d(j); %kkk就是边界点集合
            
        end
    end
    tt(i).f = kkk(i,:); %把边界点存放在结构体tt中的f列中 
    tt(i).f(tt(i).f == 0) = [];%去除边界点集合中的0
end
% tt.f为边界点，tt.c为聚类中心点，tt.d为该列的邻近集


%% 计算通道test
vvv = struct([]); %用来放通道
for i = 1:m
    ww = 0;
    nn = 0;
    for j = 1:m
        if intersect(tt(i).f,tt(j).d) > 0 %判断边界点与其他邻近集是否有交集
            ww = ww+1;
            vvv(i).a = i;
            vvv(i).b(ww) = j; %将i簇与j簇所对应的列储存起来,有多个的时候放在一块
            vvv(i).c(ww) = {intersect([tt(i).f],[tt(j).d])}; %存放通道，用元胞的方式来区分每个数组对应的列，改列对应的边界点与其他邻近及的交集
        end
        if intersect(tt(i).d,tt(j).f) > 0
            nn = nn+1;
            vvv(i).e(nn) = j;
            vvv(i).d(nn) = {intersect(tt(i).d,tt(j).f)};% 该列对应的邻近集与其他边界点的交集
        end
    end
    vvv(i).f = union(vvv(i).b,vvv(i).e);

    vvv(i).f(vvv(i).f == i) = [];
end



for i = 1:m
    d = tt(i).f; %该列的边界点
    f = vvv(i).f; % 与该列相交的列
    ppp = zeros(1,length(f));
    for j = 1:length(f)
        g = tt(f(j)).f;% 与该列相交的列的边界点
        ddd = zeros(1,length(d));
        for z = 1:length(d)
            disg =zeros(1,length(g));
            for zz = 1:length(g)
                disg(zz) =  norm(data(d(z),:)-data(g(zz),:));

            end
            a = (disg~=0);
            b = sum(a(:));
            if b > 0
                disg(disg==0) = [];%去除通道为0的情况
                [~,temp2] = sort(disg);
                ddd(z) = disg(temp2(1));
            else
                ddd(z) = nan;
            end
        end
        ppp(1,j) = min(ddd);
    end
    vvv(i).h = ppp; %ppp为该列对应其他列的最小通道距离
end

for i = 1:m
    tt(i).g = vvv(i).f; %该列与其他列边界点相交的集合 
    tt(i).h = vvv(i).h; %所对应的通道距离 
end
        
% 对距离进行排序
for i = 1:m
    xxx = zeros(1,length(tt(i).g));
    [~,xxx] = sort(tt(i).h,'MissingPlacement','last');
    tt(i).g = tt(i).g(xxx);
    tt(i).h = sort(tt(i).h,'MissingPlacement','last');
end
% 消除NaN
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

sss = struct([]); %用来放通道
for i =1:m
    sss(i).a = tt(i).j;
    sss(i).b = tt(i).i;
    sss(i).c = {};
    for j = 1:length(sss(i).b)
        sss(i).c(j) = {[i,sss(i).b(j)]};
    end
end

ge = [sss.a];%放的是要融合的两个簇的边界点的聚类
go = [sss.c];% 放的是要融合的两个簇
[ge1,ge2] = sort(ge);% 对通道距离进行排序
go1 = go(ge2); % 对应的通道也按照升序排列
ger = ge1';
go2 =go1';
biu = cell2mat(go2);
get = cat(2,ger,biu); %排序后的通道距离及通道

dk = unique(get(:,1)); % 后面发现利用通道距离效果不好，不如文章中利用的NI来进行通道融合
%%  已通道边界的NI来反向排序
gett = get; 
for i = 1:length(gett)
    gett(i,4) = max(rr(tt(gett(i,2)).b).f,rr(tt(gett(i,3)).b).f);
end
[gev1,gev2] = sort(gett(:,4),"descend");
getx = gett(gev2,:);
dkk = unique(getx(:,4));
dkk = sort(dkk,"descend");%这里的getx变量，第一列为通道距离，第二、第三列为要融合的簇，第三列为这两个簇要融合的簇的边界点的NI的最小值，因为需要对NI从大到小进行顺序融合，因此用了方向排序。
%% 
Q = {[]};
for i = 1:length(dkk)
    [row,col] = find(getx==dkk(i));
    U = getx(row,2);
    Q(i) = {unique(U)'};
end
%% 融合集合对应的邻近集
B = {};
for i = 1:length(Q)
    B(i) = {unique([tt(cell2mat(Q(i))).d])'};
end
%% 融合
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
        K(j) = W(i); % K为最终聚类效果
        j = j+1;
    end
end
output = K;