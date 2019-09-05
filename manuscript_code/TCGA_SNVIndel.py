import os,sys,re
import numpy as np
import subprocess
from multiprocessing import Process,Pool
import pickle
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqIO.FastaIO import SimpleFastaParser
import math
from math import log, exp
from getNatureA import Neoantigen
import glob
import pickle
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from collections import defaultdict
import tempfile

hydro_score = dict(R=-4.5,K=-3.9,N=-3.5,D=-3.5,Q=-3.5,E=-3.5,H=-3.2,P=-1.6,Y=-1.3,W=-0.9,S=-.8,
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



def getH(pep,mymodel):
    hydro_vector = do_hydro_vector(pep)
    if len(hydro_vector) == 8:
        H = .5
        return H
    else:
        H = round(mymodel[len(hydro_vector)].predict_proba(hydro_vector)[0][1],3)
        return H


def do_hydro_vector(pep):
    hydro_vector = []
    for i in pep:
        hydro_vector.append(hydro_score[i.upper()])
    return hydro_vector

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

def getiedbseq():
    iedb_seq = []
    with open('/home/wzt/project/GeneFusion/iedb.fasta', 'r') as fin:
        for t, seq in SimpleFastaParser(fin):
            iedb_seq.append(seq)
    return iedb_seq

def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x-ma),v))) + ma

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

def cal_netctlpan(hla,pep):
    pepfa = tempfile.NamedTemporaryFile(mode='w+',dir='.')
    netctl = tempfile.NamedTemporaryFile(mode='w+',dir='.')
    pepfa.write('>pep\n{}\n'.format(pep))
    pepfa.flush()
    cmd = 'python2  /home/zhouchi/software/netchop/predict.py --method netctlpan --length {}  --allele {}  ' \
          '--epitope_threshold -99.9  --threshold -99.9  --noplot  {} > {} '.format(len(pep),hla,pepfa.name,netctl.name)
    os.system(cmd)
    BOOL = True
    while BOOL:
        lines = subprocess.getoutput('grep {} {}'.format(pep,netctl.name)).split('\t')
        if len(lines) == 0:
            os.system(cmd)
        else:
            pepfa.close()
            netctl.close()
            return lines[3],lines[4],lines[5]

def isfloat(x):
    try:
        return float(x)
    except:
        return 0

def main(file):
    sample = file.split('.')[0]
    fileout = '{}.filter.mupexi'.format(sample)
    with open(file,'r') as fin,open(fileout,'w') as fout:
        fout.write('HLA\tSample\tmismatch\t')
        fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
        fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
        fout.write('Hydro_Model\tSelf_similar\tR\tA\n')
        for line in fin:
            if line.startswith('#') or line.startswith('HLA_allele'):
                continue
            lines = line.strip().split('\t')
            if isfloat(lines[11]) <= 2 and isfloat(lines[26]) * isfloat(lines[16]) >= 1 and lines[25] == 'No':
                hla = lines[0]
                mtpep,mtpep_score, mtpep_aff, mtpep_rank = lines[7],lines[12],lines[10],lines[11]
                wtpep,wtpep_score, wtpep_aff, wtpep_rank = lines[1], lines[6], lines[4], lines[5]
                mtpep_tap, mtpep_cleavage, mtpep_comb = cal_netctlpan(hla, mtpep)
                wtpep_tap, wtpep_cleavage, wtpep_comb = cal_netctlpan(hla, wtpep)
                A = Neoantigen(mtpep, wtpep, mtpep_aff, wtpep_aff).getA()
                R = getR(mtpep, iedb_seq)
                H = getH(mtpep, mymodel)
                self_similar = cal_similarity(mut_pep=mtpep, wild_pep=wtpep)
                mismatch = lines[17]
                out = '\t'.join((hla, sample, str(mismatch), mtpep, mtpep_score, mtpep_aff, mtpep_rank, mtpep_comb,
                wtpep, wtpep_score,wtpep_aff, wtpep_rank, wtpep_comb,str(H),str(self_similar),str(R),str(A)))
                fout.write('{}\n'.format(out))
def domain(file):
    try:
        main(file)
    except:
        print (file)

CCs = ['BLCA','BRCA','CESC','COAD','GBM','HNSC','KICH','KIRC','KIRP','LIHC',
           'LUAD','LUSC','OV','PAAD','PRAD','READ','SKCM','STAD','THCA','UCEC']
#CC = sys.argv[1]
for CC in CCs:
    wd = '/home/wzt/project/GeneFusion/SNVIndel/{}/MuPeXI'.format(CC)
    os.chdir(wd)
    mymodel = getmodel()
    iedb_seq = getiedbseq()
    temp = glob.glob('*mupexi')
    temp = [i for i in temp if 'filter' not in i]
    files = []
    for file in temp:                 ##把已经计算过得去除掉
        sample = file.split('.')[0]
        fileout = '{}.filter.mupexi'.format(sample)
        if not os.path.isfile(fileout):
            files.append(file)
    pool = Pool(processes=40)
    pool.map(domain,files)
    pool.close()
    pool.join()