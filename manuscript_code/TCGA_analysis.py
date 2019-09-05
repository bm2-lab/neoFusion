import numpy as np,pandas as pd
import os,sys,glob,re
from scipy import stats
from collections import defaultdict
CCs = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
        'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD',
        'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

kmers = [9,10,11]

def f1(Type):
    mydict = {}
    myScore= {'FRAMESHIFT':[],'INFRAME':[]}
    with open('Fusion_Score.tsv','r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if len(lines[3]) not in kmers:
                continue
            sample = lines[1].split('_')[1]
            if Type == 'Score':
                mydict[sample] = mydict.get(sample,0) + float(lines[-1])  ##其他指标,如neoantigen的个数
            elif Type == 'Counts':
                mydict[sample] = mydict.get(sample, 0) + 1  ##其他指标,如neoantigen的个数
    with open('Fusion.filter.txt','r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if lines[6] == 'INFRAME' or  lines[6] =='FRAMESHIFT':
                myScore[lines[6]].append(mydict.get(lines[1],0))
    print (stats.mannwhitneyu(myScore['FRAMESHIFT'], myScore['INFRAME'],alternative='greater'))
    return myScore


def f2():
    myScore = f1('Score')
    myCounts = f1('Counts')
    with open('frame.tsv', 'w') as fout:
        fout.write('Type1\tType2\tNums\n')
        for i in myScore['INFRAME']:
            fout.write('Score\tINFRAME\t{}\n'.format(i))
            #fout.write('Score\tINFRAME\t{}\n'.format(-np.log10(i+ 10**-20)))
        for i in myScore['FRAMESHIFT']:
            fout.write('Score\tFRAMESHIFT\t{}\n'.format(i))
            #fout.write('Score\tFRAMESHIFT\t{}\n'.format(-np.log10(i+10**-20)))
        for i in myCounts['INFRAME']:
            fout.write('Counts\tINFRAME\t{}\n'.format(i))
        for i in myCounts['FRAMESHIFT']:
            fout.write('Counts\tFRAMESHIFT\t{}\n'.format(i))

    dat = pd.read_csv('frame.tsv', sep='\t', header=0)
    fig, ax = plt.subplots()
    sns.boxplot(x='Type1', y='Nums', hue='Type2', data=dat, showfliers=False,ax=ax)
    plt.legend(framealpha=0)
    fig.savefig('frame.pdf', dpi=100, bbox_inches='tight')




def f2():
    os.chdir('/home/wzt/project/GeneFusion/TCGA')
    mydict = {}
    myScore = {'Drive': [], 'nonDrive': []}
    genetype = 'Onco'
    with open('Fusion_Score.tsv', 'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            sample = lines[1].split('_')[1]
            if len(lines[3]) not in kmers:
                continue
            mydict[sample] = mydict.get(sample, 0) + float(lines[-1])
    with open('Fusion.filter.txt', 'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if lines[-2] == 'NO':
                continue
            genetypes = lines[-1].split(',')
            if genetype in genetypes and lines[1] in mydict:
                myScore['Drive'].append(mydict.get(lines[1], 0))
            elif lines[-1] == 'Normal' and lines[1] in mydict:
                myScore['nonDrive'].append(mydict.get(lines[1], 0))
    print(stats.stats.mannwhitneyu(myScore['Drive'], myScore['nonDrive'],alternative='less'))
    print(np.mean(myScore['Drive']))
    print(np.mean(myScore['nonDrive']))



def f3(ratio=True):
    os.chdir('/home/wzt/project/GeneFusion/TCGA')
    mydict = defaultdict(dict)
    with open('Fusion.filter.txt') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            genetype = lines[-1].split(',')
            for i in genetype:
                mydict[i][lines[-4]] = mydict[i].get(lines[-4],0) + 1
    if ratio is True:
        for i in mydict:
            S = mydict[i]['FRAMESHIFT'] + mydict[i]['INFRAME'] + mydict[i]['NOFRAME']
            for j in ['FRAMESHIFT','INFRAME','NOFRAME']:
                mydict[i][j] = round((mydict[i][j]) / S,3)
    a = [mydict['Onco']['INFRAME'],mydict['TSG']['INFRAME'],mydict['Kinase']['INFRAME'],mydict['Normal']['INFRAME']]
    b = [mydict['Onco']['FRAMESHIFT'], mydict['TSG']['FRAMESHIFT'], mydict['Kinase']['FRAMESHIFT'], mydict['Normal']['FRAMESHIFT']]
    c = [mydict['Onco']['NOFRAME'], mydict['TSG']['NOFRAME'], mydict['Kinase']['NOFRAME'], mydict['Normal']['NOFRAME']]

    barWidth = 0.85
    names = ('Onco', 'TSG', 'Kinase', 'Normal')
    r = [0, 1, 2, 3]
    p1 = plt.bar(r, a, edgecolor='white', width=barWidth)
    p2 = plt.bar(r, b, bottom=a, edgecolor='white', width=barWidth)
    p3 = plt.bar(r, c, bottom=[i + j for i, j in zip(a, b)], width=barWidth)
    plt.xticks(r, names)
    plt.xlim(-1, 6)
    plt.ylim(0, 1.3)
    plt.legend((p1[0], p2[0], p3[0]), ('inframe', 'frameshift', 'noframe'))
    plt.savefig('frameratio.pdf', dpi=100, bbox_inches='tight')



def f4():
    def getScorelist(mydict, k):
        mylow = []
        myhigh = []
        for gene in mydict:
            if mydict[gene]['Counts'] >= k:
                myhigh.append(mydict[gene]['Score'])
            elif mydict[gene]['Counts'] == 1:
                mylow.append(mydict[gene]['Score'])
        return myhigh, mylow
    os.chdir('/home/wzt/project/GeneFusion/TCGA')
    myScore = {}
    mydict = defaultdict(dict)
    with open('Fusion_Score.tsv', 'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            sample = lines[1].split('_')[1]
            myScore[sample] = myScore.get(sample, 0) + float(lines[-1])
    with open('Fusion.filter.txt','r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if lines[-2] == 'NO':
                continue
            #genes = lines[0].split('--')    ####  gene没有score也要算
            # if lines[1] not in myScore:
            #     continue
            mydict[lines[0]]['Counts'] = mydict[lines[0]].get('Counts',0) + 1
            mydict[lines[0]]['Score'] = mydict[lines[0]].get('Score',0) + myScore.get(lines[1],0)
    for i in mydict:
        mydict[i]['Score'] = mydict[i]['Score'] / mydict[i]['Counts']


    myhigh, mylow = getScorelist(mydict, 10)
    c = []
    for i in range(1000):
        a = myhigh
        b = np.random.choice(mylow, size=2000, replace=False)
        if np.mean(a) < np.mean(b):
            c.append(1)
        else:
            c.append(0)
    print(sum(c))


def f5():

    def getScorelist(mydict, k):     ### 以pep计算,一个样品可以因为HLA的原因可以出现几次
        mylow = []
        myhigh = []
        for pep in mydict:
            if  mydict[pep]['Counts'] >= k:
                myhigh.append(mydict[pep]['Score'])
            elif mydict[pep]['Counts'] == 1:
                mylow.append(mydict[pep]['Score'])
        return myhigh, mylow
    def getScoreTime(mydict,k):
        mylist = []
        for pep in mydict:
            if mydict[pep]['Counts'] == k:
                mylist.append(mydict[pep]['Score'])
        return mylist
    mydict = defaultdict(dict)
    with open('Fusion_Score.tsv', 'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            sample  = lines[1][:7] + '_' + lines[1][13:]
            if lines[3] not in mydict:
                mydict[lines[3]]['Sample'] = set()
            mydict[lines[3]]['Score'] = mydict[lines[3]].get('Score', 0) + float(lines[-1])
            mydict[lines[3]]['Counts'] = mydict[lines[3]].get('Counts', 0) + 1
            mydict[lines[3]]['Sample'].add(sample)
        for i in mydict:
            mydict[i]['Score'] = mydict[i]['Score'] / mydict[i]['Counts']
    myhigh,mylow = getScorelist(mydict,20)
    mylist1 = getScoreTime(mydict, 1)
    mylist2 = getScoreTime(mydict, 2)
    mylist3 = getScoreTime(mydict, 3)
    mylist4, _ = getScorelist(mydict, 4)
    mylist11 = [- np.log10(i + 10 ** -20) for i in mylist1]
    mylist22 = [- np.log10(i + 10 ** -20) for i in mylist2]
    mylist33 = [- np.log10(i + 10 ** -20) for i in mylist3]
    mylist44 = [- np.log10(i + 10 ** -20) for i in mylist4]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, figsize=(8, 5), sharey=True)
    sns.boxplot(mylist11, ax=ax1, showfliers=True, orient='v')
    sns.boxplot(mylist22, ax=ax2, showfliers=True, orient='v')
    sns.boxplot(mylist33, ax=ax3, showfliers=True, orient='v')
    sns.boxplot(mylist44, ax=ax4, showfliers=True, orient='v')
    fig.savefig('TimesScore.pdf', dpi=100, bbox_inches='tight')
    return myhigh,mylow,mydict



def f7():
    fileout = '/home/wzt/project/GeneFusion/SNVIndel/FusionSNVIndel_MSI.tsv'
    with open(fileout,'w') as fout:
        fout.write('Sample\tCC\tMSI_MSS\tType\tCounts\n')
        CCs = ['COAD','STAD','UCEC']
        filemsi = '/home/wzt/project/GeneFusion/SNVIndel/MSI.txt'
        filein = '/home/wzt/project/GeneFusion/TCGA/Fusion.filter.txt'
        mydict = {}
        with open(filein,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                mydict[lines[3][5:12]] = mydict.get(lines[3][5:12],0) + 1
        with open(filemsi,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if lines[1][5:] in CCs and float(lines[2]) >= 0.4:
                    temp = mydict.get(lines[0][5:],0)
                    fout.write('{}\t{}\tMSI\tFusion\t{}\n'.format(lines[0][5:],lines[1][5:],temp))
                elif lines[1][5:] in CCs and float(lines[2]) < 0.4:
                    temp = mydict.get(lines[0][5:], 0)
                    fout.write('{}\t{}\tMSS\tFusion\t{}\n'.format(lines[0][5:], lines[1][5:], temp))
        filemsi = '/home/wzt/project/GeneFusion/SNVIndel/MSI.txt'
        filein = '/home/wzt/project/GeneFusion/SNVIndel/FusionSNVIndel_Counts.tsv'
        mydict = {}
        with open(filemsi, 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                mydict[lines[0][5:]] = float(lines[2])
        with open(filein, 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if lines[1] == 'SNVIndel' and lines[2] in CCs:
                    if lines[0] in mydict and mydict[lines[0]] >= 0.4:
                        fout.write('{}\t{}\tMSI\tSNVIndel\t{}\n'.format(lines[0], lines[2], lines[3]))
                    elif lines[0] in mydict and mydict[lines[0]] < 0.4:
                        fout.write('{}\t{}\tMSS\tSNVIndel\t{}\n'.format(lines[0], lines[2], lines[3]))

def f8():
    os.chdir('/home/wzt/project/GeneFusion/TCGA/shuffle_hla')
    Onco = pd.read_table('OncoScoreCompare.tsv', header=0)
    TSG = pd.read_table('TSGScoreCompare.tsv', header=0)
    TPG = pd.read_table('TPGScoreCompare.tsv', header=0)
    stats.mannwhitneyu(TPG.MeanScore, Onco.MeanScore, alternative='greater')   #### 显著性差异
    stats.mannwhitneyu(TPG.MeanScore, TSG.MeanScore, alternative='greater')    #### 无显著差异






def f9():
    fileout = '/home/wzt/project/GeneFusion/SNVIndel/FusionSNVIndel_Drive_MSI.tsv'
    with open(fileout,'w') as fout:
        fout.write('Sample\tCC\tMSI_MSS\tType\tCounts\n')
        CCs = ['COAD','STAD','UCEC']
        filemsi = '/home/wzt/project/GeneFusion/SNVIndel/MSI.txt'
        filein = '/home/wzt/project/GeneFusion/TCGA/Fusion.filter.txt'
        mydict = {}
        with open(filein,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                genetypes = lines[-1].split(',')
                if 'Drive' in genetypes:
                    mydict[lines[3][5:12]] = mydict.get(lines[3][5:12],0) + 1
        with open(filemsi,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if lines[1][5:] in CCs and float(lines[2]) >= 0.4:
                    temp = mydict.get(lines[0][5:],0)
                    fout.write('{}\t{}\tMSI\tFusion\t{}\n'.format(lines[0][5:],lines[1][5:],temp))
                elif lines[1][5:] in CCs and float(lines[2]) < 0.4:
                    temp = mydict.get(lines[0][5:], 0)
                    fout.write('{}\t{}\tMSS\tFusiont\t{}\n'.format(lines[0][5:], lines[1][5:], temp))
