import os,sys,subprocess,re
import glob
import string
from Bio.SeqIO.FastaIO import SimpleFastaParser

### call fusion and neoantigen

SRRs = ['SRR925706','SRR925718','SRR925715','SRR925723','SRR925708',
        'SRR925703','SRR925697','SRR925698','SRR925704','SRR925705',
        'SRR934643','SRR934642','SRR934635','SRR934634','SRR925736',
        'SRR925726']

def fqdump():
    for SRR in SRRs:
        os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/{}'.format(SRR))
        cmd = 'fastq-dump --split-3 /home/wzt/ncbi/public/sra/{}.sra &'.format(SRR)
        subprocess.call(cmd,shell=True)

def CleanData():
    for SRR in SRRs:
        cmd = 'java -jar /home/wzt/bin/trimmomatic-0.36.jar PE -threads 32 -phred33 {}_1.fastq  {}_2.fastq ' \
          '{}.clean_1.fq {}.rm_1.fq {}.clean_2.fq {}.rm_2.fq  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ' \
          'LEADING:20 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:20 '.format(SRR,SRR,SRR,SRR,SRR,SRR,SRR)
        subprocess.call(cmd,shell=True)


def starFusion(SRR):
    cmd = 'docker run -d -v /home/wzt/project/GeneFusion/Cell_line/Data2/:/home/wzt/project/GeneFusion/Cell_line/Data2/ ' \
          '-v /home/wzt/database/STAR_Fusion/GRCh38/ctat:/home/wzt/database/STAR_Fusion/GRCh38/ctat ' \
          '--name {}  trinityctat/ctatfusion  /usr/local/src/STAR-Fusion/STAR-Fusion --left_fq  /home/wzt/project/GeneFusion/Cell_line/Data2/{}/{}_1.fastq ' \
          '--right_fq /home/wzt/project/GeneFusion/Cell_line/Data2/{}/{}_2.fastq  --genome_lib_dir /home/wzt/database/STAR_Fusion/GRCh38/ctat/  ' \
          '--FusionInspector validate --examine_coding_effect --output_dir  /home/wzt/project/GeneFusion/Cell_line/Data2/{}/STAR-Fusion ' \
          '--CPU 30'.format(SRR,SRR,SRR,SRR,SRR,SRR)
    subprocess.call(cmd,shell=True)
    cmd = 'docker wait {}'.format(SRR)
    subprocess.call(cmd,shell=True)
    cmd = 'docker rm -f {}'.format(SRR)
    subprocess.call(cmd,shell=True)



def hla(SRR):
    os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/{}'.format(SRR))
    if not os.path.isfile('fished_1.fq'):
        cmd = 'razers3 --percent-identity 95 --max-hits  1  -tc 16 --distance-range 0 -o fished_1.bam ' \
              '/home/wzt/software/OptiType/data/hla_reference_rna.fasta ' \
              '{}_1.fastq'.format(SRR)
        subprocess.call(cmd,shell=True)
        cmd = 'razers3 --percent-identity 95 --max-hits  1  -tc 16 --distance-range 0 -o fished_2.bam ' \
              '/home/wzt/software/OptiType/data/hla_reference_rna.fasta ' \
              '{}_2.fastq'.format(SRR)
        subprocess.call(cmd,shell=True)
        cmd = 'samtools fastq fished_1.bam > fished_1.fq'
        subprocess.call(cmd,shell=True)
        cmd = 'samtools fastq fished_2.bam > fished_2.fq'
        subprocess.call(cmd,shell=True)

    if not os.path.isfile('hla_type.txt'):
        cmd = 'python2 /home/zhouchi/software/OptiType/OptiTypePipeline.py  -i fished_1.fq  fished_2.fq --rna  -o HLA'
        subprocess.call(cmd,shell=True)
        dirs = glob.glob('HLA/*')[0]
        hla = glob.glob('HLA/*/*tsv')[0]
        pdf = glob.glob('HLA/*/*pdf')[0]
        cmd = 'cp {} HLA/hla.txt'.format(hla)
        subprocess.call(cmd,shell=True)
        cmd = 'cp {} HLA/coverage.pdf'.format(pdf)
        subprocess.call(cmd,shell=True)
        cmd = 'rm -r {}'.format(dirs)
        subprocess.call(cmd,shell=True)
        with open('HLA/hla.txt','r') as fin:
            data = fin.read()
        lines = data.strip().split('\n')[1].split('\t')[1:-2]
        hla_type = []
        for line in lines:
            line_ = 'HLA-'+line.replace('*','').replace(':','')
            hla_type.append(line_)
        hla_type = list(set(hla_type))
        with open('hla_type.txt','w') as fout:
            for i in hla_type:
                fout.write('{}\n'.format(i))



def dofilter():
    for SRR in SRRs:
        myset = set()
        os.chdir('/home/wzt/project/GeneFusion/Cell_line/Data2/results/{}'.format(SRR))
        indexs = [i + j for i in string.ascii_uppercase for j in string.ascii_uppercase]
        index = 0
        filein1 = '../../{}/STAR-Fusion/FusionInspector-validate/finspector.fusion_predictions.final.abridged.FFPM.annotated.coding_effect'.format(SRR)
        filein2 = '../../{}/STAR-Fusion/star-fusion.fusion_predictions.abridged.coding_effect.tsv'.format(SRR)
        fileout = '{}_fusion.filter.txt'.format(SRR)
        with open(fileout, 'w') as fout:
            if os.path.isfile(filein1):
                with open(filein1, 'r') as fin:
                    for line in fin:
                        lines = line.strip().split('\t')
                        if line.startswith('#') or lines[28] == "." or lines[10] == 'NO' or float(lines[19]) <= 0.1:
                            continue
                        else:
                            pep = lines[-3].split('*')[0]
                            if pep in myset:
                                continue
                            fout.write('>Neo_{}\n{}\n'.format(indexs[index], pep))
                            index += 1
                            myset.add(pep)

            if os.path.isfile(filein2):
                with open(filein2,'r') as fin:
                    for line in fin:
                        lines = line.strip().split('\t')
                        if line.startswith('#') or lines[-3] == '.':
                            continue
                        else:
                            pep = lines[-3].split('*')[0]
                            if pep in myset:
                                continue
                            fout.write('>Neo_{}\n{}\n'.format(indexs[index], pep))
                            index += 1
                            myset.add(pep)


