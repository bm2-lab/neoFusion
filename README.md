# neoFusion: predict fusion neoantigen and assess immunogenecity  
  
Fusion, an important class of somatic mutation, is an ideal source of tumor-derived neoantigens for the capability of creating open reading frame. In this study, the first one-top pipeline for identifying and prioritizing fusion neoantigen are presented.            
     
## Dependencies   

#### Required software:  
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [OptiType](https://github.com/FRED-2/OptiType)
* [netchop](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netchop)
* [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion)
* [ pepmatch ]  module of MuPeXI, for convenience, we provide this binary tool.


#### Python2 package     
     xgboost
     biopython
     scikit-learn
     pandas
     numpy
     subprocess
     multiprocessing
     pickle   

#### reference  
     download from (https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/)    
     please refer to STAR-Fusion for more detail.    
     All the other needed materials can be found at data directory   

## Installation   
#### Install via Docker, highly recommended   
Docker image of neoFusion is available at https://hub.docker.com/r/weizhiting/neofusion/.   
if you have docker installed, you can pull the image like so:   
    
    docker pull weizhiting/neofusion

#### Install from source, not recommended   
To install  neoFusion, you just need download neoFusion.py and all the described dependencies.


## Usage    
if you install from docker:
    
    docker run -it weizhiting/neofusion
    python /usr/local/neoFusion.py -h   

if you install from source:    

    edit the software path in neoFusion to the right path in your environment.  

Two mode are provided, denovo and midway. In denovo, you should provide fastq files.    
  
    python neoFusion denovo --left 1.fq  --right 2.fq  --hla  HLA-A02:01  --genome STAR_Fusion/GRCh38/ctat/

In midway, you can provide fusion transcripts in fa format.    
    
    python neoFusion midway --fusion fusion.fa --hla HLA-A02:01     

## output 
The output file "neoscore.txt" contains all putative neoantigens information.

| Column | Description |
| - | -: |
| HLA | HLA type |
| mismatch | The base number different between mtpep and wtpep |
| MTpep | fusion derived neopeptide |
| MTpep_score | netMHCpan predicted score of mtpep |
| MTpep_aff | netMHCpan predicted affinity of mtpep |
| MTpep_rank | netMHCpan predicted bind rank(%) of mtpep |
| MTpep_comb | netCTLpan predicted comb score of mtpep, combined MHC score, cleavage score and TAP score |
| WTpep | pepmatch_db_x86_64 extracted normal peptide |
| WTpep_score |  netMHCpan predicted score of wtpep |
| WTpep_aff | netMHCpan predicted affinity of wtpep |
| WTpep_rank | netMHCpan predicted bind rank(%) of wtpep |
| WTpep_comb | netCTLpan predicted comb score of wtpep, combined MHC score, cleavage score and TAP score |
| Hydro_Model | peptide immunity measurement based on amino acid hydrophobicity |
| R | T cell recognition score calculated based on TCR cross-reactivity |
| Score | Immunogenicity score of neoantigens |

## Contact   

Zhiting Wei   
1632738@tongji.edu.cn

Biological and Medical Big data Mining Lab    
Tongji University    
Shanghai, China.

## neoFusion flowchart
![](workflow.png)


