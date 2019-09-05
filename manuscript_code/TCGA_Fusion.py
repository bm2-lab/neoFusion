import os,re,sys,glob,subprocess,pickle,warnings
import pandas as pd,numpy as np
from subprocess import PIPE
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool
from math import log, exp
from getNatureA import Neoantigen
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from collections import defaultdict
warnings.filterwarnings('ignore')
hydro_score = dict(R=-4.5,K=-3.9,N=-3.5,D=-3.5,Q=-3.5,E=-3.5,H=-3.2,P=-1.6, Y=-1.3,W=-0.9,S=-.8,
                   T=-0.7,G=-0.4,A=1.8,M=1.9,C=2.5,F=2.8,L=3.8,V=4.2,I=4.5)

def getR(neo_seq,iedb_seq):     ### nature文章的相似性
    align_score = []
    a = 26
    k = 4.86936
    for seq in iedb_seq:
        aln_score = aligner(neo_seq,seq)
        if aln_score:
            localds_core = max([line[2] for line in aln_score])
            align_score.append(localds_core)

    bindingEnergies = list(map(lambda x: -k * (a - x), align_score))
    lZk = logSum(bindingEnergies + [0])
    lGb = logSum(bindingEnergies)
    R=exp(lGb-lZk)
    return R

def getwildpep(mut_pep):
    path = '/home/zhouchi/software/MuPeXI/bin/pepmatch_db_x86_64'
    with open('/tmp/tmp.fa','w') as fout:
        fout.write('{}\n'.format(mut_pep))
    ref_seq = '/home/zhouchi/database/Annotation/references/reference_peptide_{}.txt'.format(len(mut_pep))
    cmd = '{}  -thr 10  temp/{}.fa  {}  > temp/{}.tsv'.format(path,temp,ref_seq,temp)
    os.system(cmd)
    lines = os.popen("grep -v '#' temp/{}.tsv".format(temp)).read().strip().split()
    return lines[3]




def aligner(seq1,seq2):
    matrix = matlist.blosum62
    gap_open = -11
    gap_extend = -1
    aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
    return aln

def cal_similarity(mut_pep,wild_pep):    ### 序列相似性
    score_pair = aligner(mut_pep,wild_pep)[0][2]
    score_self = aligner(mut_pep, mut_pep)[0][2]
    score = score_pair / score_self
    return round(score,3)

def cal_netctlpan(hla,pep,index):
    temp = index
    with open('temp/{}.fa'.format(temp), 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'python2  /home/zhouchi/software/netchop/predict.py --method netctlpan --length {}  --allele {}  ' \
          '--epitope_threshold -99.9  --threshold -99.9  --noplot  temp/{}.fa > temp/{}.tsv'.format(len(pep),hla,temp,temp)
    os.system(cmd)
    BOOL = True
    while BOOL:
        lines = os.popen('grep {} temp/{}.tsv'.format(pep,temp)).read().strip().split()
        if len(lines) == 0:
            os.system(cmd)
        else:
            return lines[3],lines[4],lines[5]

def cal_netMHCpan(hla,pep,index):
    temp = index
    with open('temp/{}.fa'.format(temp), 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'netMHCpan  -a {}  -l {} -f temp/{}.fa -s  -BA > temp/{}.tsv'.format(hla,len(pep),temp,temp)
    os.system(cmd)
    results = os.popen('grep {} temp/{}.tsv'.format(pep,temp)).read().strip().split()
    return results[11],results[12],results[13]

def tranform(x):
    return 50000**(1-x)


def isWild(pep):
    ref_pep = '/home/wzt/database/STAR_Fusion/GRCh38/ctat/ref_annot_fusion.pep'
    fin = open(ref_pep, 'r')
    title, ref_seq = next(SimpleFastaParser(fin))
    fin.close()
    if ref_seq.find(pep) == -1:
        return False
    else:
        return True

def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x-ma),v))) + ma

def cal_mismatch(mtpep,wtpep):
    mis = 0
    for i,j in zip(mtpep,wtpep):
        if i != j:
            mis += 1
    return mis

def getiedbseq():
    iedb_seq = []
    with open('/home/wzt/project/GeneFusion/iedb.fasta', 'r') as fin:
        for t, seq in SimpleFastaParser(fin):
            iedb_seq.append(seq)
    return iedb_seq

def getmodel():
    mymodel = {}
    path = '/home/wzt/project/GeneFusion/dbGAP/results/'
    model_9 = path + 'hg_xgb_9.dat'
    model_10 = path + 'hg_xgb_10.dat'
    model_11 = path + 'hg_xgb_11.dat'
    with open(model_9, 'rb') as fin1, open(model_10, 'rb') as fin2, open(model_11, 'rb') as fin3:
        mymodel[9] = pickle.load(fin1)
        mymodel[10] = pickle.load(fin2)
        mymodel[11] = pickle.load(fin3)
    return mymodel

def do_hydro_vector(pep):
    hydro_vector = []
    for i in pep:
        hydro_vector.append(hydro_score[i.upper()])
    return hydro_vector

def getH(pep,mymodel):
    hydro_vector = do_hydro_vector(pep)
    if len(hydro_vector) == 8:
        H = .5
        return H
    else:
        H = round(mymodel[len(hydro_vector)].predict_proba(hydro_vector)[0][1],3)
        return H

def main(line):
    lines = line.strip().split('\t')
    index = lines[-1]
    hla = lines[0].replace('*','')
    sample = lines[2] + '_' + lines[7]
    mtpep = lines[1]
    mtpep_score, mtpep_aff, mtpep_rank = lines[3],lines[4],lines[5]
    wtpep = getwildpep(mtpep)
    wtpep_score, wtpep_aff, wtpep_rank = cal_netMHCpan(hla, wtpep,index)
    mtpep_tap, mtpep_cleavage, mtpep_comb = cal_netctlpan(hla, mtpep,index)
    wtpep_tap, wtpep_cleavage, wtpep_comb = cal_netctlpan(hla, wtpep,index)
    A = Neoantigen(mtpep, wtpep, mtpep_aff, wtpep_aff).getA()
    R = getR(mtpep,iedb_seq)
    H = getH(mtpep,mymodel)
    self_similar = cal_similarity(mut_pep=mtpep, wild_pep=wtpep)
    mismatch = cal_mismatch(mtpep, wtpep)
    out = '\t'.join((hla,sample,str(mismatch),mtpep,mtpep_score, mtpep_aff, mtpep_rank,mtpep_comb,
    wtpep, wtpep_score, wtpep_aff, wtpep_rank, wtpep_comb,str(H),str(self_similar),str(R),str(A)))
    out = out + '\n'
    return out
filein = '/home/wzt/project/GeneFusion/TCGA/netMHC.filter_wild.pep'
fileout = '/home/wzt/project/GeneFusion/TCGA/Fusion_Score.tsv'

with open(filein,'r') as fin, open(fileout,'w') as fout:
    fout.write('HLA\tSample\tmismatch\t')
    fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
    fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
    fout.write('Hydro_Model\tSelf_similar\tR\tA\n')
    iedb_seq = getiedbseq()
    mymodel = getmodel()
    fin.readline()
    lines = [line.strip() for line in fin]
    pool = Pool(processes=48)
    results = pool.map(main,lines)
    pool.close()
    pool.join()
    for result in results:
        fout.write(result)

