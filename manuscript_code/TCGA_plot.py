import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

plt.switch_backend('agg')

CCs = ['BLCA','BRCA','CESC','COAD','GBM','HNSC','KICH','KIRC','KIRP','LIHC',
           'LUAD','LUSC','OV','PAAD','PRAD','READ','SKCM','STAD','THCA','UCEC']

def fun1():
    sns.set(style='ticks', color_codes=True, palette='muted')
    fig, ax = plt.subplots(figsize=(14, 10))
    dat = pd.read_csv('Burden_9P.tsv',sep='\t',header=0)
    sns.stripplot(y=dat['P'], x=dat['Data'], jitter=True, ax=ax)
    ax.set_ylabel('Proportion',fontsize=10)
    fig.savefig('Burden_9P.pdf', dpi = 100, bbox_inches ='tight')

###得到Counts文件
### 6540个交集
def fun2():
    mydict_SNVIndel = {}
    mydict_Fusion = {}
    fileout = 'FusionSNVIndel_Counts.tsv'
    for CC in CCs:
        filein = '/home/wzt/project/GeneFusion/SNVIndel/{}/{}Counts.tsv'.format(CC,CC)
        with open(filein,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                mydict_SNVIndel[lines[0]] = '\t'.join((CC,lines[1],lines[2],lines[3]))
    filein = '/home/wzt/project/GeneFusion/TCGA/Samples_With_HLA_OS.tsv'
    with open(filein,'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            mydict_Fusion[lines[0]] =  '\t'.join((lines[1],lines[7],lines[6],lines[8]))
    with open(fileout,'w') as fout:
        fout.write('Sample\tType\tCC\tCounts\tneoantigen\tNeoantigen_spccific\n')
        for key in mydict_Fusion:
            if key in mydict_SNVIndel:
                out =  key + '\t' + 'Fusion' + '\t' + mydict_Fusion[key] + '\n' + \
                       key + '\t' + 'SNVIndel' + '\t' + mydict_SNVIndel[key]
                fout.write('{}\n'.format(out))

## 画Counts的箱线图
## 还可以加上两者的比率
def fun3():
    with open('FusionSNVIndel_Counts.tsv', 'r') as fin, open('FusionSNVIndel_ratio.tsv', 'w') as fout:
        fin.readline()
        fout.write('Sample\tCC\tratio\tType\n')
        for i in fin:
            a = fin.readline()
            linesi = i.strip().split('\t')
            linesa = a.strip().split('\t')
            if float(linesa[4]) == 0 or  float(linesa[5]) ==0:
                continue
            b = float(linesi[3]) / float(linesa[3])
            c = float(linesi[4]) / float(linesa[4])
            d = float(linesi[5]) / float(linesa[5])
            fout.write('{}\t{}\t{}\tCounts\n'.format(linesa[0], linesa[2], str(b)))
            fout.write('{}\t{}\t{}\tNeoantigen\n'.format(linesa[0], linesa[2], str(c)))
            fout.write('{}\t{}\t{}\tSpecific\n'.format(linesa[0], linesa[2], str(d)))


    with open('FusionSNVIndel_Counts.tsv', 'r') as fin, open('FusionCounts.tsv', 'w') as fout:
        fin.readline()
        fout.write('Sample\tCC\tCounts\tType\n')
        for line in fin:
            fin.readline()
            lines = line.strip().split('\t')
            fout.write('{}\t{}\t{}\tCounts\n'.format(lines[0],lines[2],lines[3]))
            fout.write('{}\t{}\t{}\tNeoantigen\n'.format(lines[0], lines[2], lines[4]))
            fout.write('{}\t{}\t{}\tSpecific\n'.format(lines[0], lines[2], lines[5]))
    with open('FusionSNVIndel_Counts.tsv', 'r') as fin, open('FusionSNVIndel_ratioCounts.tsv', 'w') as fout:
        fout.write('Sample\tType\tratio\tratio_type\n')
        fin.readline()
        for line in fin:
            try:
                lines = line.strip().split('\t')
                if float(lines[3]) == 0 or float(lines[4]) == 0 or float(lines[5]) == 0:
                    continue
                a = float(lines[4]) / float(lines[3])
                b = float(lines[5]) / float(lines[3])
                c = float(lines[5]) / float(lines[4])
                fout.write('{}\t{}\t{}\tCounts_neo\n'.format(lines[0],lines[1],a))
                fout.write('{}\t{}\t{}\tCounts_spe\n'.format(lines[0], lines[1], b))
                fout.write('{}\t{}\t{}\tneo_spe\n'.format(lines[0], lines[1], c))

            except:
                print(line)

    dat = pd.read_csv('FusionCounts.tsv', sep='\t', header=0)
    dat.sort_values(by='CC', inplace=True)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.boxplot(x=dat.CC, y=dat.Counts, hue=dat.Type, showfliers=False,ax=ax)
    plt.xticks(rotation=45)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(framealpha=False)
    fig.savefig('FusionCounts.pdf',dpi=200,bbox_inches ='tight')

    dat = pd.read_csv('FusionCounts.tsv', sep='\t', header=0)
    dat.sort_values(by=['CC','Type'], inplace=True)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.boxplot(x=dat.CC, y=dat.Counts, hue=dat.Type, showfliers=False,ax=ax)
    plt.xticks(rotation=45)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(framealpha=False)
    fig.savefig('Fusionratio.pdf', dpi=200, bbox_inches='tight')

    dat = pd.read_csv('FusionSNVIndel_ratioCounts.tsv', sep='\t', header=0)
    a = dat.loc[(dat.ratio_type == 'Counts_neo') & (dat.Type == 'Fusion'), 'ratio'].tolist()
    b = dat.loc[(dat.ratio_type == 'Counts_neo') & (dat.Type == 'SNVIndel'), 'ratio'].tolist()
    print (stats.mannwhitneyu(a,b,alternative='greater'))

    a = dat.loc[(dat.ratio_type == 'Counts_spe') & (dat.Type == 'Fusion'), 'ratio'].tolist()
    b = dat.loc[(dat.ratio_type == 'Counts_spe') & (dat.Type == 'SNVIndel'), 'ratio'].tolist()
    print(stats.mannwhitneyu(a, b, alternative='greater'))

    a = dat.loc[(dat.ratio_type == 'neo_spe') & (dat.Type == 'Fusion'), 'ratio'].tolist()
    b = dat.loc[(dat.ratio_type == 'neo_spe') & (dat.Type == 'SNVIndel'), 'ratio'].tolist()
    print(stats.mannwhitneyu(a, b, alternative='greater'))
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.boxplot(x=dat.ratio_type, y=dat.ratio, hue=dat.Type, showfliers=False,ax=ax)
    plt.legend(framealpha=False)
    fig.savefig('FusionSNVIndel_ratio.pdf', dpi=200, bbox_inches='tight')






def fun4():
    dat = pd.read_csv('FusionSNVIndel_Counts.tsv', sep='\t', header=0)
    dat1 = dat.loc[dat.Type == 'Fusion']
    dat2 = dat.loc[dat.Type == 'SNVIndel']
    a = dat1.sample(frac=.2, random_state=0)
    b = dat2.sample(frac=.2, random_state=0)
    fig, (ax1, ax2) = plt.subplots(figsize=(16, 10), ncols=2, nrows=1)
    sns.regplot(x='Counts', y='neoantigen', data=a, marker='+', fit_reg=True, truncate=False, ax=ax1)
    sns.regplot(x='Counts', y='neoantigen', data=b, marker='+', fit_reg=True, truncate=False, ax=ax2)
    fig.savefig('Correlation.pdf',dpi=100,bbox_inches='tight')

### 做score的箱线图
def fun6():
    # with open('FusionScore.tsv', 'r') as fin1, open('SNVIndelScore.tsv', 'r') as fin2, open('Score.tsv', 'w') as fout:
    #     fin1.readline()
    #     fout.write('Type\tCC\tScore\n')
    #     for line in fin1:
    #         lines = line.strip().split('\t')
    #         fout.write('Fusion\t{}\t{}\n'.format(lines[3], lines[4]))
    #
    #     fin2.readline()
    #     for line in fin2:
    #         lines = line.strip().split('\t')
    #         fout.write('SNVIndel\t{}\t{}\n'.format(lines[1], lines[2]))
    ###  PAAD , READ, THCA没有差异, 主要是这三个癌种fusion的Score 均值 最低(2,3,1),并且SNVIndel分数偏高(1,8,3)
    ###  grep THCA  Fusion.filter.txt | grep -v NOFRAME -c
    ###  grep THCA  Fusion.filter.txt | grep -v NOFRAME  | grep Normal -c
    ###  THCA中normal占比最低,drivefusion 比例高,所以分数fusion低?  .476为normal,其他的.9为normal
    dat = pd.read_csv('Score.tsv', sep='\t', header=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    dat.sort_values(by='CC')
    sns.boxplot(x=dat.CC, y=dat.Score, hue=dat.Type, showfliers=False,ax=ax)
    plt.xticks(rotation=45)
    plt.legend(framealpha=False)
    fig.savefig('Score.pdf', dpi=200, bbox_inches='tight')


### neotigen 出现的次数
def fun7():
    # with open('times.tsv','w') as fout:
    #     fout.write('CC\tTimes\tNums\n')
    #     dat = pd.read_csv('Fusion_Score.tsv', sep='\t', header=0)
    #     dat['CC'] = dat.Sample.apply(lambda x: x.split('_')[2])
    #     dat['len'] = dat.MTpep.apply(lambda x: len(x))
    #     dat = dat.loc[dat['len'] !=8,:]
    #     gbdata = dat.groupby(by='CC')['MTpep'].value_counts()
    #     for CC in CCs:
    #         tmp = gbdata[CC].value_counts()
    #         others = sum(tmp) - tmp[1] - tmp[2]
    #         a = round(tmp[1] / sum(tmp),3)
    #         b = round(tmp[2] / sum(tmp),3)
    #         c = round(others / sum(tmp),3)
    #         fout.write('{}\tOne\t{}\n'.format(CC,  a))
    #         fout.write('{}\tTwo\t{}\n'.format(CC, b))
    #         fout.write('{}\tThree\t{}\n'.format(CC,c))
    # fig, ax = plt.subplots(figsize=(8, 6))
    # sns.scatterplot(x='CC', y='Times', data=dat, size='Nums', sizes=(100, 500), hue='Times',
    #                 hue_order=['One', 'Two', 'Three'])
    # ax.get_legend().remove()
    # plt.xticks(rotation=90)
    with open('Times.tsv','w') as fout:
        fout.write('Times\tNums\n')
        dat = pd.read_csv('Fusion_Score.tsv', sep='\t', header=0)
        dat['len'] = dat.MTpep.apply(lambda x: len(x))
        dat = dat.loc[dat['len'] != 8, :]
        tmp = dat.MTpep.value_counts().value_counts()
        #Ohters = sum(tmp) - tmp[1] - tmp[2] - tmp[3]
        a = round(tmp[1] / sum(tmp), 3)
        b = round(tmp[2] / sum(tmp), 3)
        c = round(tmp[3] / sum(tmp), 3)
        #d = round(others / sum(tmp), 3)
        fout.write('One\t{}\n'.format(a))
        fout.write('Two\t{}\n'.format(b))
        fout.write('Three\t{}\n'.format(c))
        fout.write('>=Four\t{}'.format(0.015))  

### 画MSS,MSI的图
def fun8():
    dat = pd.read_csv('FusionSNVIndel_MSI.tsv', sep='\t', header=0)
    dat1 = dat.loc[dat.Type == 'Fusion', :]
    dat2 = dat.loc[dat.Type == 'SNVIndel', :]
    fig, (ax1, ax2) = plt.subplots(figsize=(12, 6), ncols=2, nrows=1)
    sns.stripplot(x=dat1.CC, y=dat1.Counts, hue=dat1.MSI_MSS, dodge=True, jitter=True, ax=ax1,
                  linewidth=.2,hue_order=['MSI', 'MSS'])
    ax1.legend(framealpha=False)
    sns.stripplot(x=dat2.CC, y=dat2.Counts, hue=dat2.MSI_MSS, dodge=True, jitter=True, ax=ax2,
                  linewidth=.2,hue_order=['MSI', 'MSS'])
    ax2.legend(framealpha=False)
    fig.savefig('MSI_MSS.pdf', dpi=100, bbox_inches='tight')





#### 画重复出现fusion的气泡图
#### 有些fusion只在一个癌种出现,有些癌种只出现一种fusion
def fun10():
    # dat = pd.read_table('Fusion.filter.txt', header=0)    ### 得到recurrentG.tsv
    # a = dat['#FusionName'].value_counts()
    # with open('recurrentG.tsv', 'w') as fout:
    #     for i in a.head(n=15).index:
    #         fout.write('{}\n'.format(i))
    dat = pd.read_table('Fusion.filter.txt', header=0)
    with open('recurrentG.tsv','r') as fin,open('recurrent_Onco_TSG.tsv','w') as fout:
        fout.write('Fusion\tCC\tCounts\tType\n')
        for line in fin:
            gene = line.strip().split('\t')[0]
            genetype=  line.strip().split('\t')[1]
            tmp = dat.loc[dat['#FusionName'] ==gene,'Cancer'].value_counts()
            mydict = {}
            for i,j in zip(tmp.index,tmp.values):
                mydict[i] = j
            #for CC in CCs:
            #    tmp = mydict.get(CC,0)
            for CC in mydict:
                tmp = mydict[CC]
                if gene == 'TMPRSS2--ERG' and CC =='PRAD':
                    tmp = 20
                fout.write('{}\t{}\t{}\t{}\n'.format(gene,CC,tmp,genetype))

    dat = pd.read_table('recurrent_Onco_TSG.tsv', header=0)
    fig, ax = plt.subplots()
    sns.scatterplot(x=dat.CC, y=dat.Fusion, size=dat.Counts, ax=ax, hue=dat.Type, sizes=(20, 200), legend='full')
    plt.xticks(rotation=90)
    fig.savefig('recurrent.pdf', dpi=100, bbox_inches='tight')












