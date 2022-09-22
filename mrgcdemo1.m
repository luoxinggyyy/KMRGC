function [output] = mrgcdemo(input,k)  % input为输入数据，k为最邻近系数，output为聚类结果，并以元胞的形式保存。
data = input;
[m,n] = size(data); % m = 300,n = 2
pattern = data; %数据制作完成
%% 计算欧式距离(计算所有点的距离)
dist1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        dist1(i,j) = norm(data(i,:)-data(j,:));% 距离公式,（dist（i，j）代表第i行与第j行的距离）
    end
end
%% 对距离进行排序
k1 = k;
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
        if rr(yy(j,t)).f >= rr(t).f
            KRI(t) = rr(t).f * j * p(j+1,t); %KRI为核心指数，这里是t列的邻近集有大于它的NI，就执行这个
            break
        end
    end
    if max([rr(yy(:,t)).f]) < rr(t).f
        KRI(t) = rr(t).f*0.0001; % 这里是为了避免KRI为0的情况，后续都为0，采用的是NI乘与一个很小的数，使得不对前面的KRI产生影响同时有大小差异。
    end
end
figure
ykri = sort(KRI,'descend');  % 降序排列
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
    tt(i).e = KRI(tt(i).c);
    tt(i).m = rr(i).f;
end

%% 计算簇
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
    YY = [];% 空数组，用来做每次簇的初始化
    for z = 1:m
        YY = union(YY,[tt(uu(CT(z))).d]); %表示u分别与列对应的邻近集求并集再赋值给u，起到吞并累加的作用
        if length(u) + z < m
            if length(intersect(YY,[tt(uu(CT(z+1))).d])) == 0 %如果没有交集，就退出循环。（意思就是簇的融合完毕，执行下一个较大的簇进行融合）
                break
            end
        else
            break
        end
    end
    hhh(i) = {YY'};
    u = cell2mat(hhh);
    u = unique(u);
    if length(u) >= m % 当所有簇包含的列的总和与点的总数相同时，循环结束
        break
    end
end
output = hhh;