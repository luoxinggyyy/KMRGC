# KMRGC
multi-resolution graph-based clustering and its impoving method（KMRGC）
Instructions for use
（1）. KMRGC.m is the code of KMRGC, and the specific principle is in the code comment, and MRGC.m is the code of MRGC, 
    and the principle is also in the code comment.
（2）. test.mlx is the running code to call K-means, MRGC and KMRGC algorithms, and the datasets are Dataset Ⅰ and Dataset Ⅱ in the article, 
    and the production of the datasets is also in the test.xlm file.
（3）. If you want to run a custom dataset, please open test1.mlx, change the address of the dataset in the xlread function inside, 
    and you can run it (Note: data normalization is used uniformly here, if you do not want to normalize, remember to delete the normalization code in the test1.mlx file)

tips:The algorithm is written in matlab platform, I use the genuine matlab 2022b version, if you run a function error please update the latest version.
 Luo Xin
