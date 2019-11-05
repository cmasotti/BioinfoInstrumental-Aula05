# Bioinformática Instrumental - Módulo II
# Identificação de Variantes Genéticas

## Aula05: Chamada de Variantes e Controle de Qualidade
Seguiremos utilizando as boas práticas de chamada de variantes preconizadas pelo pacote de ferramentas [GATK](https://software.broadinstitute.org/gatk/) (_Genome Analysis Toolkit_).

Nesta etapa utilizaremos os dados "limpos", pré-processados (após ApplyBQSR), para iniciar a chamada de variantes propriamente dita. Identificaremos SNPs e INDELs de linhagens tumorais e germinativas para produzir um único conjunto de dados no formato [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (_variant call format_). Por fim, filtraremos variantes potencialmente artefatuais com base em escores de qualidade gerados pelo software de chamada de variantes e pelos demais programas até aqui utlizados.
O workflow completo abordado nesta aula prática:
!(calling)[]








