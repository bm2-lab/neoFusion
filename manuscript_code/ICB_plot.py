import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
plt.switch_backend('agg')


### 两个数据的比率图

# f, ax = plt.subplots(figsize=(16, 8))
# fig = sns.boxplot(x=var, y="SalePrice", data=data)
# fig.axis(ymin=0, ymax=800000)
# plt.xticks(rotation=90)


#### 画总体概览图
def f3():
    #file1 = '/home/wzt/project/GeneFusion/Data2/results/OS_91011_Score_CTL.tsv'
    #file2 = '/home/wzt/project/GeneFusion/Data2/results/OS_91011_Burden_CTL.tsv'
    file1 = '/home/wzt/project/GeneFusion/Data1/results/PFS_91011_Score2_CTL.tsv'
    file2 = '/home/wzt/project/GeneFusion/Data1/results/PFS_91011_Burden_CTL.tsv'
    dat = pd.read_csv(file1,sep='\t',index_col=0,header=0)
    dat1 = pd.read_csv(file2,sep='\t',index_col=0,header=0)
    dat.sort_values(by='SNVIndel',inplace=True,ascending=False)
    dat1 = dat1.loc[dat.index]
    width = 0.05
    n = dat.shape[0]
    ind = np.arange(0, width * n, width)  # the x locations for the groups
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=1, nrows=3, gridspec_kw={'height_ratios': [10, 1, 10]})
    fig.set_figwidth(5)
    ax1.bar(ind, np.log2(dat.SNVIndel + 1), width, edgecolor=['black'] * n)
    ax1.bar(ind, np.log2(dat.Fusion + 1), width, bottom=np.log2(dat.SNVIndel + 1), edgecolor=['black'] * n)
    ax1.set_xticks([])
    # ax1.set(yscale='log')
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    color = []
    #a = 'CR'
    #b = 'PR'
    a = 'long-survival'
    b = 'response'
    for i in dat.group:
        if i == a:
            color.append('maroon')
        elif i == b:
            color.append('red')
        else:
            color.append('black')
    ax2.bar(ind, [0.095] * n, width, color=color, edgecolor=['white'] * n)
    ax2.set_xticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_yticks([])
    #ax2.set_ylim(0, 0.1)
    Fusion = [np.log2(i + 1) for i in dat1.Fusion]
    Fusion = [i * -1 for i in Fusion]
    SNVIndel = [np.log2(i + 1) for i in dat1.SNVIndel]
    SNVIndel = [i * -1 for i in SNVIndel]
    ax3.bar(ind, SNVIndel, width, edgecolor=['black'] * n)
    ax3.bar(ind, Fusion, width, bottom=SNVIndel, edgecolor=['black'] * n)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_xticks([])
    fig.subplots_adjust(hspace=0)
    fig.savefig('overview.pdf',dpi=100, bbox_inches='tight')


#### 画箱线图
def f4():
    file1 = '/home/wzt/project/GeneFusion/Data1/results/PFS_91011_Score2_CTL.tsv'
    file2 = '/home/wzt/project/GeneFusion/Data1/results/PFS_91011_Burden_CTL.tsv'
    dat = pd.read_csv(file1, sep='\t', index_col=0, header=0)
    dat['Sum'] = dat.Fusion + dat.SNVIndel
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.ylim(0, 40)
    sns.boxplot(x=dat.group, y=dat.Sum, showfliers=False)
    fig.savefig('box_Score.pdf', dpi=100, bbox_inches='tight')

    dat = pd.read_csv(file2, sep='\t', index_col=0, header=0)
    dat['Sum'] = dat.Fusion + dat.SNVIndel
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.ylim(0, 900)
    sns.boxplot(x=dat.group, y=dat.Sum, showfliers=False)
    fig.savefig('box_burden.pdf', dpi=100, bbox_inches='tight')

    file1 = '/home/wzt/project/GeneFusion/Data2/results/OS_91011_Score_CTL.tsv'
    file2 = '/home/wzt/project/GeneFusion/Data2/results/OS_91011_Burden_CTL.tsv'
    dat = pd.read_csv(file1, sep='\t', index_col=0, header=0)
    dat['Sum'] = dat.Fusion + dat.SNVIndel
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.ylim(0, 30)
    sns.boxplot(x=dat.group, y=dat.Sum, showfliers=False)
    fig.savefig('box_Score.pdf', dpi=100, bbox_inches='tight')

    dat = pd.read_csv(file2, sep='\t', index_col=0, header=0)
    dat['Sum'] = dat.Fusion + dat.SNVIndel
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.ylim(0, 1000)
    sns.boxplot(x=dat.group, y=dat.Sum, showfliers=False)
    fig.savefig('box_burden.pdf', dpi=100, bbox_inches='tight')





