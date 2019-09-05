import os,re,glob
import numpy as np
from collections import defaultdict
import subprocess
import pandas as pd

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

Pts = ['Pt1', 'Pt2', 'Pt4', 'Pt5', 'Pt6', 'Pt7', 'Pt9', 'Pt10', 'Pt12', 'Pt13', 'Pt14', 'Pt15', 'Pt16',
       'Pt19', 'Pt20', 'Pt22', 'Pt23', 'Pt25', 'Pt28', 'Pt29', 'Pt31', 'Pt32', 'Pt35', 'Pt37', 'Pt38']



kmers = [9,10,11]


def cal(x):
    return 1 / (1 + np.exp(5 *(x-2)))

def f1(Rm,Rn,self_similar,H,R,A,mismatch,comb):
    score = Rm * (1 - Rn / 2 ** mismatch) * R * comb * H
    return  score

def getfusionScore():
    mydict = {}
    for Pt in Pts:
        filein = '/home/wzt/project/GeneFusion/Data2/{}/fusion_score.tsv'.format(Pt)
        if not os.path.isfile(filein):
            continue
        with open(filein,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if len(lines[3]) not in kmers:
                    continue
                Rm = cal(float(lines[6]))
                Rn = cal(float(lines[11]))
                R = float(lines[15])
                A = float(lines[16])
                self_similar,H,mismatch, comb = map(float,[lines[14],lines[13],lines[2],lines[7]])
                score = f1(Rm=Rm,Rn=Rn,self_similar = self_similar,H=H,R=R,mismatch=mismatch,comb=comb,A=A)
                if score < 0:
                    score = 0
                mydict[Pt] = mydict.get(Pt,0) + score
    return mydict


def getSNVIndelScore():
    mydict = {}
    for Pt in Pts:
        filein = '/home/wzt/project/GeneFusion/Data2/{}/{}.filter.mupexi'.format(Pt,Pt)
        with open(filein,'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if len(lines[3]) not in kmers:
                    continue
                Rm = cal(float(lines[6]))
                Rn = cal(float(lines[11]))
                R =  float(lines[15])
                A =  float(lines[16])
                self_similar,H,mismatch, comb = map(float,[lines[14],lines[13],lines[2],lines[7]])
                score = f1(Rm=Rm,Rn=Rn,self_similar = self_similar,H=H,R=R,mismatch=mismatch,comb=comb,A = A)
                if score < 0:
                    score = 0
                mydict[Pt] = mydict.get(Pt,0) + score
    return mydict

def main():
    myPFS = {}
    fileout ='/home/wzt/project/GeneFusion/Data2/results/OS_91011_Score_CTL.tsv'
    fusion_dict = getfusionScore()
    SNVIndel_dict = getSNVIndelScore()
    filein = '/home/wzt/project/GeneFusion/Data2/Sample_State_1.tsv'
    with open(filein,'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            if lines[0] == 'Pt8':
                continue
            myPFS[lines[0]] = '\t'.join((lines[3],str(float(lines[1])/30),lines[2],lines[4],lines[5],lines[6]))  ### CTL
    with open(fileout,'w') as fout:
        fout.write('Patient\tgroup\tOS\tEvent\tCTL\tGender\tAge\tFusion\tSNVIndel\n')
        for Pat in Pts:
            Fusion_score = str(fusion_dict.get(Pat,0))
            SNV_score = str(SNVIndel_dict.get(Pat,0))
            fout.write('{}\t{}\t{}\t{}\n'.format(Pat,myPFS[Pat],Fusion_score,SNV_score))

main()
