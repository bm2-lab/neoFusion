### 查看高Fusion和低Fusion样品之间 NMD的不同
import os,re,sys,glob,subprocess



mydict = {}


def getExpress(mylist,filein):
    express_list = []
    for ensembl_gene in mylist:
        cmd = "grep {} {}".format(ensembl_gene,filein)
        output = subprocess.getoutput(cmd).strip().split('\t')
        if len(output) !=2:
            express = str(0)
        else:
            express = output[1]
        express_list.append(express)
    return express_list


def main1():
    mylist1 = []    ####  geneName
    mylist2 = []    ####  ensembl gene

    os.chdir('/home/wzt/project/GeneFusion/SNVIndel/NMD')
    with open('gene.txt','r') as fin:
        for line in fin:
            lines = line.strip().split('\t')
            mylist1.append(lines[0])
            mylist2.append(lines[1])

    with open('NMD.txt','r') as fin,open('express.txt','w') as fout:
        fout.write('Sample\t{}\n'.format('\t'.join(mylist1)))
        for line in fin:
            lines = line.strip().split()
            tmpfile = lines[13]
            sampleName = lines[7]
            express_list = getExpress(mylist2,tmpfile)
            fout.write('{}\t{}\n'.format(sampleName,'\t'.join(express_list)))

def main2():
    mydict = {}
    os.chdir('/home/wzt/project/GeneFusion/SNVIndel/NMD')
    with open('/home/wzt/project/GeneFusion/SNVIndel/FusionCounts.tsv','r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if lines[3] == 'Counts':
                mydict[lines[0]] = lines[1] + '\t' + lines[2]
    with open('express.txt','r') as fin,open('express1.txt','w') as fout:
        fout.write('{}\tCC\tCounts\n'.format(fin.readline().strip()))
        for line in fin:
            lines = line.strip().split('\t',maxsplit=1)
            if lines[0] in mydict:
                fout.write('{}\t{}\t{}\n'.format(lines[0],lines[1],mydict[lines[0]]))





##main1()   ### 获得表达量
main2()



