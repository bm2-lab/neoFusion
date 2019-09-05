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

def getWildpep(mut_pep):
    path = '/home/zhouchi/software/MuPeXI/bin/pepmatch_db_x86_64'
    with open('/tmp/tmp.fa','w') as fout:
        fout.write('{}\n'.format(mut_pep))
    ref_seq = '/home/zhouchi/database/Annotation/references/reference_peptide_{}.txt'.format(len(mut_pep))
    cmd = '{}  -thr 10  /tmp/tmp.fa  {}  > /tmp/tmp.tsv'.format(path,ref_seq)
    subprocess.call(cmd,shell=True)
    cmd = "grep -v '#' /tmp/tmp.tsv"
    line = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE)
    lines = line.stdout.readline().decode().strip().split()
    return lines[3],lines[5]

def cal(x):
    return 1 / (1 + np.exp(5 *(x-2)))

def f1(Rm,Rn,H,RA,mismatch,comb):
    score = Rm * (1 - Rn / 2 ** mismatch) * R * comb * H
    return  score

def f2(Rm,Rn,mismatch,comb):
    score = Rm * comb
    return  score



def getR(neo_seq,iedb_seq):
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

def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x-ma),v))) + ma

def getIEDBseq():
    iedb_seq = []
    with open('/home/wzt/project/GeneFusion/iedb.fasta', 'r') as fin:
        for t, seq in SimpleFastaParser(fin):
            iedb_seq.append(seq)
    return iedb_seq

def getModel():
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



def cal_netctlpan(hla,pep):
    with open('/tmp/tmp.fa', 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'python2  /home/zhouchi/software/netchop/predict.py --method netctlpan --length {}  --allele {}  ' \
          '--epitope_threshold -99.9  --threshold -99.9  --noplot  /tmp/tmp.fa > /tmp/tmp.tsv'.format(len(pep),hla)
    subprocess.call(cmd,shell=True)
    while True:
        line = subprocess.Popen('grep {} /tmp/tmp.tsv'.format(pep),shell=True,stdin=PIPE,stdout=PIPE)
        lines = line.stdout.readline().decode().strip().split()
        if len(lines) == 0:
            os.system(cmd)
        else:
            return lines[3],lines[4],lines[5]

def cal_netMHCpan(hla,pep):
    with open('/tmp/tmp.fa', 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'netMHCpan -a {}  -l {} -f /tmp/tmp.fa -s  -BA > /tmp/tmp.tsv'.format(hla,len(pep))
    subprocess.call(cmd,shell=True)
    cmd = 'grep {} /tmp/tmp.tsv'.format(pep)
    result = subprocess.Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
    results = result.stdout.readline().decode().strip().split()
    return results[11],results[12],results[13]


def gethla(SRR):
    myset = set()
    filein = '/home/wzt/project/GeneFusion/Cell_line/Data2/{}/hla_type.txt'.format(SRR)
    with open(filein,'r') as fin:
        for line in fin:
            line = line.strip()
            line = line[:-2] + ':' + line[-2:]
            myset.add(line)
    return ','.join(list(myset))


def netMHCpan(SRR):
    os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/{}'.format(SRR))
    filein = '{}_fusion.filter.txt'.format(SRR)
    if os.path.isfile(filein) and os.path.getsize(filein) > 0:
        hlas = gethla(SRR)
        cmd = 'netMHCpan -a {} -l 9,10,11 -f {} -s -BA > netMHC.txt'.format(hlas, filein)
        subprocess.call(cmd, shell=True)

def parsernetMHCpan(SRR):
    print('\n######## Running  neoAntigen Scoring  ############')
    os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/{}'.format(SRR))
    filein = 'netMHC.txt'
    fileout = '/tmp/{}NeoScore.txt'.format(SRR)
    ref = ''
    iedb_seq = getIEDBseq()
    mymodel  = getModel()
    with open('/home/wzt/database/STAR_Fusion/GRCh38/ctat/ref_annot_fusion.pep','r') as fin:
        for t,seq in SimpleFastaParser(fin):
            ref = seq
            break
    with open(filein,'r') as fin,open(fileout,'w') as fout:
        fout.write('HLA\tSample\tmismatch\t')
        fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
        fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
        fout.write('Hydro_Model\tR\tScore1\tScore2\n')
        for line in fin:
            if re.search('=\s+WB',line) or re.search('=\s+SB',line):
                lines = line.strip().split()
                mtpep = lines[2]
                if ref.find(mtpep) != -1 or len(mtpep) == 8:
                    continue
                sample = lines[10]
                mtpep, mtpep_score, mtpep_aff, mtpep_rank = lines[2], lines[11], lines[12], lines[13]
                wtpep, mismatch = getWildpep(mtpep)
                hla = lines[1].replace('*', '')
                wtpep_score, wtpep_aff, wtpep_rank = cal_netMHCpan(hla, wtpep)
                mtpep_tap, mtpep_cleavage, mtpep_comb = cal_netctlpan(hla, mtpep)
                wtpep_tap, wtpep_cleavage, wtpep_comb = cal_netctlpan(hla, wtpep)
                A = Neoantigen(mtpep, wtpep, mtpep_aff, wtpep_aff).getA()
                R = getR(mtpep, iedb_seq)
                H = getH(mtpep, mymodel)
                Rm = cal(float(mtpep_rank))
                Rn = cal(float(wtpep_rank))
                RA = float(A) * float(R)
                score1 = f1(Rm=Rm,Rn=Rn,H=float(H),RA=RA,mismatch=float(mismatch),comb=float(mtpep_comb))
                if score1 < 0 :
                    score1 = 0

                score2 = f2(Rm=Rm,Rn=Rn,mismatch=float(mismatch),comb=float(mtpep_comb))
                if score2 < 0 :
                    score2 = 0
                out = '\t'.join((hla, sample, str(mismatch), mtpep, mtpep_score, mtpep_aff, mtpep_rank, mtpep_comb,
                                 wtpep, wtpep_score, wtpep_aff, wtpep_rank, wtpep_comb, str(H), str(RA)))
                fout.write('{}\t{}\t{}\n'.format(out,str(score1),str(score2)))
    filein = '/tmp/{}NeoScore.txt'.format(SRR)
    fileout = 'NeoScore.txt'
    if os.path.isfile(filein) and os.path.getsize(filein) > 0:
        dat = pd.read_csv(filein,sep='\t',header=0)
        dat.sort_values(by='Score2',ascending=False,inplace=True)
        dat.to_csv(fileout,sep='\t',header=True,index=False)
    print('\n######## Complete Successfuly  ############')

 def comet(SRR):
     print (SRR)
     os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/{}'.format(SRR))
     mzML = glob.glob('*mzML')[0]
     cmd = 'comet -P{}  -DrefWithFusion.pep -N{} {}'.format(params,prefix,mzML)
     print (cmd)
     subprocess.call(cmd,shell=True)
     cmd = 'percolator {}.pin -r  {}_percolator.txt  --override --protein-enzyme no_enzyme'.format(prefix,prefix)
     subprocess.call(cmd,shell=True)
     filein = '{}_percolator.txt'.format(prefix)
     fileout = '{}_percolator.filter.txt'.format(prefix)
     with open(filein, 'r') as fin, open(fileout, 'w') as fout:
         fout.write(fin.readline())
         for line in fin:
             lines = line.strip().split('\t')
             if float(lines[2]) <= 0.01:
                 fout.write(line)


parsernetMHCpan('SRR925705')
parsernetMHCpan('SRR925708')


