import os,re,sys,glob
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool
import string
import pickle
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqIO.FastaIO import SimpleFastaParser
import math
from math import log, exp
from getNatureA import Neoantigen
import pickle
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from collections import defaultdict
import time


import warnings
warnings.filterwarnings('ignore')
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
    with open('temp.fa', 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'python2  /home/zhouchi/software/netchop/predict.py --method netctlpan --length {}  --allele {}  ' \
          '--epitope_threshold -99.9  --threshold -99.9  --noplot  temp.fa > temp.tsv'.format(len(pep),hla)
    subprocess.call(cmd,shell=True)
    while True:
        lines = subprocess.getoutput('grep {} temp.tsv'.format(pep)).strip().split()
        if len(lines) == 0:
            os.system(cmd)
        else:
            return lines[3],lines[4],lines[5]

def cal_netMHCpan(hla,pep):
    with open('temp.fa', 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'netMHCpan -a {}  -l {} -f temp.fa -s  -BA > temp.tsv'.format(hla,len(pep))
    subprocess.call(cmd,shell=True)
    cmd = 'grep {} temp.tsv'.format(pep)
    results = subprocess.getoutput(cmd).strip().split()
    return results[11],results[12],results[13]

def tranform(x):
    return 50000**(1-x)


def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x-ma),v))) + ma


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

RNA_dict = {'SRR3184279':'Pt1','SRR3184280':'Pt2','SRR3184281':'Pt4','SRR3184282':'Pt5','SRR3184283':'Pt6',
            'SRR3184284': 'Pt7','SRR3184285':'Pt8','SRR3184286':'Pt9','SRR3184287':'Pt10','SRR3184288':'Pt12',
            'SRR3184289': 'Pt13','SRR3184290':'Pt14','SRR3184291':'Pt15','SRR3184292':'Pt16','SRR3184293':'Pt19',
            'SRR3184294': 'Pt20','SRR3184295':'Pt22','SRR3184296':'Pt23','SRR3184297':'Pt25','SRR3184300':'Pt28',
            'SRR3184301':'Pt29','SRR3184302':'Pt31','SRR3184303':'Pt32','SRR3184304': 'Pt35','SRR3184305':'Pt37',
            'SRR3184306':'Pt38'}

SRRs = ['SRR3184279', 'SRR3184280', 'SRR3184281', 'SRR3184282', 'SRR3184283', 'SRR3184284', 'SRR3184285', 'SRR3184286',
        'SRR3184287', 'SRR3184288', 'SRR3184289', 'SRR3184290', 'SRR3184291', 'SRR3184292', 'SRR3184293', 'SRR3184294',
        'SRR3184295', 'SRR3184296', 'SRR3184297', 'SRR3184300', 'SRR3184301', 'SRR3184302', 'SRR3184303', 'SRR3184304',
        'SRR3184305', 'SRR3184306']



def starFusion(sample):
    des = '/home/wzt/project/GeneFusion/Data2/{}'.format(sample)
    os.chdir(des)
    cmd = 'docker run -d -v /home/wzt/project/GeneFusion/Data2/:/home/wzt/project/GeneFusion/Data2/ ' \
          '-v /home/wzt/database/STAR_Fusion/GRCh38/ctat:/home/wzt/database/STAR_Fusion/GRCh38/ctat ' \
          '--name {}  trinityctat/ctatfusion  /usr/local/src/STAR-Fusion/STAR-Fusion --left_fq  /home/wzt/project/GeneFusion/Data2/{}/{}.clean_1.fq ' \
          '--right_fq /home/wzt/project/GeneFusion/Data2/{}/{}.clean_2.fq  --genome_lib_dir /home/wzt/database/STAR_Fusion/GRCh38/ctat/  ' \
          '--FusionInspector validate --examine_coding_effect --output_dir  /home/wzt/project/GeneFusion/Data2/{}/STAR-Fusion ' \
          '--CPU 24'.format(sample,sample,sample,sample,sample,sample)
    subprocess.call(cmd,shell=True)
    cmd = 'docker wait {}'.format(sample)
    subprocess.call(cmd,shell=True)
    cmd = 'docker rm -f {}'.format(sample)
    subprocess.call(cmd,shell=True)
    # results = 'STAR-Fusion/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect'
    # with open(results,'r') as fin,open('star-fusion.txt','w') as fout:
    #     for index,line in enumerate(fin):
    #         lines = line.strip().split('\t')
    #         if index == 0:
    #             fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
    #                        .format(lines[0],lines[1], lines[2], lines[3], lines[5], lines[6],lines[8], lines[9], lines[19], lines[21],lines[23], lines[28]))
    #         else:
    #             if lines[10] == 'YES' and float(lines[19]) >= 0.1 and lines[28] != "." and int(lines[1]) >= 2 and int(lines[2]) >= 2:
    #                 fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
    #                            .format(lines[0],lines[1], lines[2], lines[3], lines[5], lines[6],lines[8], lines[9], lines[19], lines[21],lines[23], lines[28]))
    # os.system('echo "starFusion ok" >> pipeline.log')

def dofilter():
    for SRR in SRRs:
        os.chdir('/home/wzt/project/GeneFusion/Data2/{}'.format(SRR))
        indexs = [i + j for i in string.ascii_uppercase for j in string.ascii_uppercase]
        filein = 'STAR-Fusion/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect'
        fileout = '{}_fusion.filter.txt'.format(SRR)
        if os.path.isfile(filein):
            with open(filein,'r') as fin,open(fileout,'w') as fout:
                for index,line in enumerate(fin):
                    lines = line.strip().split('\t')
                    if line.startswith('#') or lines[28] == "." or lines[10] == 'NO' or float(lines[19]) <= 0.1:
                        continue
                    else:
                        pep = lines[-3].split('*')[0]
                        fout.write('>{}\n{}\n'.format(indexs[index],pep))

def gethla(SRR):
    myset = set()
    filein = '/home/wzt/project/GeneFusion/Data2/{}/hla_type.txt'.format(SRR)
    with open(filein,'r') as fin:
        for line in fin:
            line = line.strip()
            line = line[:-2] + ':' + line[-2:]
            myset.add(line)
    return ','.join(list(myset))



def netMHCpan():
    for SRR in SRRs:
        os.chdir('/home/wzt/project/GeneFusion/Data2/{}'.format(SRR))
        filein = '{}_fusion.filter.txt'.format(SRR)
        if os.path.isfile(filein) and os.path.getsize(filein) > 0:
            hlas = gethla(SRR)
            cmd = 'netMHCpan -a {} -f {} -s -BA > {}.netMHC.txt'.format(hlas,filein,SRR)
            subprocess.call(cmd,shell=True)

def getwildpep(mut_pep):
    path = '/home/zhouchi/software/MuPeXI/bin/pepmatch_db_x86_64'
    with open('temp.fa','w') as fout:
        fout.write('{}\n'.format(mut_pep))
    ref_seq = '/home/zhouchi/database/Annotation/references/reference_peptide_{}.txt'.format(len(mut_pep))
    cmd = '{}  -thr 10  temp.fa  {}  > temp.tsv'.format(path,ref_seq)
    subprocess.call(cmd,shell=True)
    cmd = "grep -v '#' temp.tsv"
    lines = subprocess.getoutput(cmd).strip().split()
    return lines[3],lines[5]

def parseMHC(SRR):
    os.chdir('/home/wzt/project/GeneFusion/Data2/{}'.format(SRR))
    filein = '{}.netMHC.txt'.format(SRR)
    fileout = 'fusion_score.tsv'
    if os.path.isfile(filein) and os.path.getsize(filein) > 0:
        with open(fileout,'w') as fout,open(filein,'r') as fin:
            fout.write('HLA\tSample\tmismatch\t')
            fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
            fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
            fout.write('Hydro_Model\tSelf_similar\tR\tA\n')
            for line in fin:
                if re.search('=\s+WB',line) or re.search('=\s+SB',line):
                    lines = line.strip().split()
                    mtpep,mtpep_score, mtpep_aff, mtpep_rank = lines[2],lines[11], lines[12], lines[13]
                    wtpep,mismatch = getwildpep(mtpep)
                    hla = lines[1].replace('*','')
                    if int(mismatch) == 0:
                        continue
                    sample = SRR + '_' + lines[10]
                    wtpep_score, wtpep_aff, wtpep_rank = cal_netMHCpan(hla, wtpep)
                    mtpep_tap, mtpep_cleavage, mtpep_comb = cal_netctlpan(hla, mtpep)
                    wtpep_tap, wtpep_cleavage, wtpep_comb = cal_netctlpan(hla, wtpep)
                    A = Neoantigen(mtpep,wtpep,mtpep_aff,wtpep_aff).getA()
                    R = getR(mtpep,iedb_seq)
                    H = getH(mtpep, mymodel)
                    self_similar = cal_similarity(mut_pep=mtpep, wild_pep=wtpep)
                    out = '\t'.join((hla, sample, str(mismatch), mtpep, mtpep_score, mtpep_aff, mtpep_rank, mtpep_comb,
                    wtpep,wtpep_score,wtpep_aff,wtpep_rank,wtpep_comb,str(H),str(self_similar),str(R),str(A)))
                    out = out + '\n'
                    fout.write(out)
#starFusion('SRR3184279')
#dofilter()
#netMHCpan()
iedb_seq = getiedbseq()
mymodel = getmodel()

pool = Pool(processes=30)
pool.map(parseMHC,SRRs)
pool.close()
pool.join()
