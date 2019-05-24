import argparse
parser = argparse.ArgumentParser(description='It is used for FusionNeoantigen Analysis',formatter_class=argparse.RawDescriptionHelpFormatter,add_help=True)
subparsers = parser.add_subparsers(help='commands')

denovo = subparsers.add_parser('denovo',help='use STAR-Fuion to call Fusion transcript',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
Req = denovo.add_argument_group('Required')
Req.add_argument('--left',action='store',metavar='<left.fq>',required=True,help='fastq file')
Req.add_argument('--genome',action='store',metavar='<reference>',required=True,help='download from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB')
Req.add_argument('--hla',action='store',metavar='<hla_allel>',required=True,help='format HLA-A02:01,HLA-B02:01')

Opt = denovo.add_argument_group('Optional')
Opt.add_argument('--right',action='store',metavar='<rigth.fq>',help='Optional,but highly recommended')
Opt.add_argument('-o',action='store',metavar='<output_dir>',default='STAR-Fusion',help='default STAR-Fusion')
Opt.add_argument('-t',action='store',metavar='<threads>',type=int,default=4, help='default 4')
Opt.add_argument('-l',action='store',metavar='<length of pep>',type=str,default='9,10,11', help='default 9,10,11')
Opt.add_argument('-e',action='store',metavar='<FFPM>',type=float,default=0.1,help='minimum fusion fragments per million rna-seq frags')
Opt.add_argument('-p',action='store',metavar='<prefix>',type=str,default='', help='outfile prefix,defalut None')
Opt.add_argument('-f','--filter',action='store_true',default=False,help='filter fusion in non-cancer tissue and cell')
Opt.add_argument('-h','--help',action='help', help='show this help message and exit')

midway = subparsers.add_parser('midway',help='user provide the fusion transcipt',add_help=False,formatter_class=argparse.RawTextHelpFormatter)
Req = midway.add_argument_group('Required')
Req.add_argument('--fusion',action='store',metavar='<protein>',required=True,help='fusion protein in fasta format')
Req.add_argument('--hla',action='store',metavar='<hla_allel>',required=True,help='format HLA-A02:01,HLA-B02:01')

Opt = midway.add_argument_group('Optional')
Opt.add_argument('-e',action='store',metavar='<FFPM>',type=float,default=0.1,help='minimum fusion fragments per million rna-seq frags')
Opt.add_argument('-p',action='store',metavar='<prefix>',type=str,default='', help='outfile prefix,defalut None')
Opt.add_argument('-l',action='store',metavar='<length of pep>',type=str,default='9,10,11', help='default 9,10,11')
Opt.add_argument('-h','--help',action='help', help='show this help message and exit')
args = parser.parse_args()


import os,re,sys,glob,subprocess,pickle,warnings
import pandas as pd,numpy as np
from subprocess import PIPE
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool
from math import log, exp
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from collections import defaultdict
warnings.filterwarnings('ignore')
hydro_score = dict(R=-4.5,K=-3.9,N=-3.5,D=-3.5,Q=-3.5,E=-3.5,H=-3.2,P=-1.6, Y=-1.3,W=-0.9,S=-.8,
                   T=-0.7,G=-0.4,A=1.8,M=1.9,C=2.5,F=2.8,L=3.8,V=4.2,I=4.5)

def getWildpep(mut_pep):
    path = '/usr/local/src/pepmatch_db_x86_64'
    with open('/tmp/tmp.fa','w') as fout:
        fout.write('{}\n'.format(mut_pep))
    ref_seq = '/usr/local/data/reference_peptide_{}.txt'.format(len(mut_pep))
    cmd = '{}  -thr 5  /tmp/tmp.fa  {}  > /tmp/tmp.tsv'.format(path,ref_seq)
    subprocess.call(cmd,shell=True)
    cmd = "grep -v '#' /tmp/tmp.tsv"
    line = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE)
    lines = line.stdout.readline().decode().strip().split()
    return lines[3],lines[5]

def cal(x):
    return 1 / (1 + np.exp(5 *(x-2)))

def f1(Rm,Rn,H,R,mismatch,comb):
    score = Rm * (1 - Rn / 2 ** mismatch) * R * comb * H
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
    with open('/usr/local/data/IEDB.fasta', 'r') as fin:
        for t, seq in SimpleFastaParser(fin):
            iedb_seq.append(seq)
    return iedb_seq

def getModel():
    mymodel = {}
    model_9 = '/usr/local/data/hg_xgb_9.dat'
    model_10 = '/usr/local/data/hg_xgb_10.dat'
    model_11 = '/usr/local/data/hg_xgb_11.dat'
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
        H = .38
        return H
    else:
        H = round(mymodel[len(hydro_vector)].predict_proba(hydro_vector)[0][1],3)
        return H



def cal_netctlpan(hla,pep):
    with open('/tmp/tmp.fa', 'w') as fout:
        fout.write('>pep\n{}\n'.format(pep))
    cmd = 'python  /usr/local/src/netchop/predict.py --method netctlpan --length {}  --allele {}  ' \
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
    cmd = '/usr/local/src/netMHCpan-4.0/netMHCpan -a {}  -l {} -f /tmp/tmp.fa -s  -BA > /tmp/tmp.tsv'.format(hla,len(pep))
    subprocess.call(cmd,shell=True)
    cmd = 'grep {} /tmp/tmp.tsv'.format(pep)
    result = subprocess.Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
    results = result.stdout.readline().decode().strip().split()
    return results[11],results[12],results[13]



def STAR_Fusion(args):
    print ('######## Running  STAR_Fusion  ########')
    if args.right:
        cmd = '/usr/local/src/STAR-Fusion/STAR-Fusion --left_fq  {} --right_fq {}  --examine_coding_effect --verbose_level 0 ' \
              '--output_dir {} --CPU {} --genome_lib_dir {} '.format(args.left, args.right, args.o,args.t, args.genome)
    else:
        cmd = '/usr/local/src/STAR-Fusion/STAR-Fusion --left_fq  {}  --examine_coding_effect --verbose_level 0 ' \
              '--output_dir {} --CPU {} --genome_lib_dir {}'.format(args.left, args.o, args.t, args.genome)

    subprocess.call(cmd,shell=True)
    print ('######## STAR_Fusion Complete Successfully ########')

def filterFusion(args):
    filein = '{}/star-fusion.fusion_predictions.abridged.coding_effect.tsv'.format(args.o)
    fileout = '{}fusion_filter.txt'.format(args.p)
    index = 1
    if os.path.isfile(filein):
        with open(filein,'r') as fin,open(fileout,'w') as fout:
            for line in fin:
                lines = line.strip().split('\t')
                if line.startswith('#') or lines[-3] == '.' or lines[8] == 'NO_LDAS' or float(lines[9]) <= args.e:
                    continue
                if args.filter is True:
                    if isNormal(lines[0]):
                        continue
                pep = lines[-3].split('*')[0]
                fout.write('>fusion{}\n{}\n'.format(index,pep))
                index += 1

def isNormal(x):
    mylist = []
    with open('/usr/local/data/normal_fusion.tsv','r') as fin:
        for line in fin:
            mylist.append(line.strip())
    xs = x.split('--')
    if xs[0] + '_' + xs[1] in mylist or xs[1] + '_' + xs[0] in mylist:
        return True
    else:
        return False

def netMHCpan(args,filein):
    print('\n######## Running  netMHCpan  #############')
    if not os.path.isfile(filein) or  os.path.getsize(filein) ==0:
        print ('fusion_filter.txt is empty,no fusion exists')
        sys.exit(2)
    else:
        cmd = '/usr/local/src/netMHCpan-4.0/netMHCpan -a {} -f {} -l {} -s -BA > {}netMHCpan.txt'.format(args.hla,filein,args.l,args.p)
        subprocess.call(cmd,shell=True)

    cmd = 'grep Error {}netMHCpan.txt'.format(args.p)
    line = subprocess.Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
    if line.stdout.readline():
        print ('netMHCpan error,please check input for netMHCpan: HLA,protein,length')
        sys.exit(2)
    else:
        print('######## netMHCpan Complete Successfully ########')

def netMHCpanD(args):
    filein = '{}fusion_filter.txt'.format(args.p)
    netMHCpan(args,filein)


def netMHCpanM(args):
    filein = args.fusion
    netMHCpan(args,filein)



def parsernetMHCpan(args):
    print('\n######## Running  neoAntigen Scoring  ############')
    filein = '{}netMHCpan.txt'.format(args.p)
    fileout = '/tmp/{}neoScore.txt'.format(args.p)
    ref = ''
    iedb_seq = getIEDBseq()
    mymodel  = getModel()
    with open('/usr/local/data/ref_annot_fusion.pep','r') as fin:
        for t,seq in SimpleFastaParser(fin):
            ref = seq
            break
    with open(filein,'r') as fin,open(fileout,'w') as fout:
        fout.write('HLA\tSample\tmismatch\t')
        fout.write('MTpep\tMTpep_score\tMTpep_aff\tMTpep_rank\tMTpep_comb\t')
        fout.write('WTpep\tWTpep_score\tWTpep_aff\tWTpep_rank\tWTpep_comb\t')
        fout.write('Hydro_Model\tR\tScore\n')
        for line in fin:
            if re.search('=\s+WB',line) or re.search('=\s+SB',line):
                lines = line.strip().split()
                mtpep = lines[2]
                if ref.find(mtpep) != -1:
                    continue
                sample = lines[10]
                mtpep, mtpep_score, mtpep_aff, mtpep_rank = lines[2], lines[11], lines[12], lines[13]
                wtpep, mismatch = getWildpep(mtpep)
                hla = lines[1].replace('*', '')
                wtpep_score, wtpep_aff, wtpep_rank = cal_netMHCpan(hla, wtpep)
                mtpep_tap, mtpep_cleavage, mtpep_comb = cal_netctlpan(hla, mtpep)
                wtpep_tap, wtpep_cleavage, wtpep_comb = cal_netctlpan(hla, wtpep)
                R = getR(mtpep, iedb_seq)
                H = getH(mtpep, mymodel)
                Rm = cal(float(mtpep_rank))
                Rn = cal(float(wtpep_rank))
                score = f1(Rm=Rm,Rn=Rn,H=float(H),R=R,mismatch=float(mismatch),comb=float(mtpep_comb))
                if score < 0 :
                    score = 0
                out = '\t'.join((hla, sample, str(mismatch), mtpep, mtpep_score, mtpep_aff, mtpep_rank, mtpep_comb,
                                 wtpep, wtpep_score, wtpep_aff, wtpep_rank, wtpep_comb, str(H), str(R)))
                fout.write('{}\t{}\n'.format(out,str(score)))
    filein = '/tmp/{}neoScore.txt'.format(args.p)
    fileout = '{}neoScore.txt'.format(args.p)
    if os.path.isfile(filein) and os.path.getsize(filein) > 0:
        dat = pd.read_csv(filein,sep='\t',header=0)
        dat.sort_values(by='Score',ascending=False,inplace=True)
        dat.to_csv(fileout,sep='\t',header=True,index=False)
    print('\n######## Complete Successfuly  ############')

def isDenovo(args):
    try:
        if args.left:
            return True
    except:
        return False


if isDenovo(args):
    STAR_Fusion(args)
    filterFusion(args)
    netMHCpanD(args)
else:
    netMHCpanM(args)
parsernetMHCpan(args)
