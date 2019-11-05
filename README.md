# Bioinformática Instrumental - Módulo II
# Identificação de Variantes Genéticas

## Aula05: Chamada de Variantes e Controle de Qualidade
Seguiremos utilizando as boas práticas de chamada de variantes preconizadas pelo pacote de ferramentas [GATK](https://software.broadinstitute.org/gatk/) (_Genome Analysis Toolkit_).

Nesta etapa utilizaremos os dados "limpos", pré-processados (após ApplyBQSR), para iniciar a chamada de variantes propriamente dita. Identificaremos SNPs e INDELs de linhagens tumorais e germinativas para produzir um único conjunto de dados no formato [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (_variant call format_). Por fim, filtraremos variantes potencialmente artefatuais com base em escores de qualidade gerados pelo software de chamada de variantes e pelos demais programas até aqui utlizados.
O workflow completo abordado nesta aula prática:   
![calling](https://github.com/cmasotti/BioinfoInstrumental-Aula05/blob/master/pipeline_aula05.png) 

A seguir, executaremos o passo-a-passo da chamada de variantes e refinamento dos dados:

### PASSO 1: ACESSO AO SERVIDOR REMOTO BIOINFO
Para relembrar como fazer a conexão ao servidor remoto via Putty, reveja os passos [aqui](https://github.com/cmasotti/BioinfoInstrumental-Aula05/blob/master/Acesso_servidor_remoto.pdf).

### PASSO 2: ORGANIZAR ESTRUTURA DE DIRETÓRIOS NO SERVIDOR REMOTO
No prompt da linha de comando no servidor remoto, criar diretórios para pré-processamento.
```bash   
aluno30@ea046e981f34:~$ 
aluno30@ea046e981f34:~$ cd /mnt/curso/aluno30/          # diretório de trabalho   
aluno30@ea046e981f34:/mnt/curso/aluno30$ ls             # confira os arquivos e diretórios existentes   
aluno30@ea046e981f34:/mnt/curso/aluno30$ mkdir calling  # criar diretório principal de trabalho da aula5  
aluno30@ea046e981f34:/mnt/curso/aluno30$ cd calling 
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir hg38   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir gvcf   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir combineGVCFs   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir genotype   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir hardFilters   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir leftNormalization  
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ mkdir annotation   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ ls     # confira os diretórios criados   
annotation  combineGVCFs  genotype  gvcf  hardFilters  hg38  leftNormalization   
```   
Criar link simbólico para o genoma de referência na pasta hg38/:
```bash   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ cd hg38   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling/hg38$ ln -s /mnt/dados/aula4/hg38/* .  
aluno30@ea046e981f34:/mnt/curso/aluno30/calling/hg38$ ls      # confira os arquivos do hg38 salvos
hg38.dict  hg38.fa  hg38.fa.amb  hg38.fa.ann  hg38.fa.bwt  hg38.fa.fai  hg38.fa.pac  hg38.fa.sa  
```  

### PASSO 3: SALVAR E CONFERIR ARQUIVOS DE DADOS PRÉ-PROCESSADOS 
Na aula passada, realizamos o passo-a-passo da organização e da "limpeza" dos dados para hoje poder prosseguir com a chamada de variantes. No último passo do pré-processamento realizamos o Base Quality Score Recalibration (BQSR) em duas etapas: (1) Detectando sistematicamente os erros de aferição de qualidade da base com [BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php) e (2) aplicando os novos escores aos dados com [ApplyBQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php).   

```bash   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling$ cd gvcf      # salvar dados pré-processados no diretório do próximo passo   
aluno30@ea046e981f34:/mnt/curso/aluno30/calling/gvcf$ ln -s ../../preprocessing/bqsr/TCGA-BH-A1F0-11B_BRCA_bqsr.ba* .  
aluno30@ea046e981f34:/mnt/curso/aluno30/calling/gvcf$ ln -s ../../preprocessing/bqsr/TCGA-BH-A1F0-01A_BRCA_bqsr.ba* .  
aluno30@ea046e981f34:/mnt/curso/aluno30/calling/gvcf$ ls      # conferir os arquivos salvos
TCGA-BH-A1F0-01A_BRCA_bqsr.bai  TCGA-BH-A1F0-01A_BRCA_bqsr.bam  TCGA-BH-A1F0-11B_BRCA_bqsr.bai  TCGA-BH-A1F0-11B_BRCA_bqsr.bam ```   

### PASSO 4: IDENTIFICAÇÃO DE VARIANTES GENÉTICAS 
Identificação de SNPs e INDELs germinativas em um ou mais indivíduos em conjunto num único arquivo VCF.




























