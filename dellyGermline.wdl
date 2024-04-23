version 1.0

struct inputSamples {
  File bam
  File bai
  File vcf
}

struct dellyResources {
    String modules
    String referenceGenome
    String dellyExclude
    String genomeMap
}

workflow dellyGermline {
  input {
    String reference
    String outputFileNamePrefix
    Array[inputSamples] inputSamples
  }

  parameter_meta {
    reference: "The genome reference build. for example: hg19, hg38"
    inputSamples: "Collection of BAM, BAI, and VCF files for >= 20 samples"
  }

  Map[String,dellyResources] resources = {
    "hg19": {
      "modules": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg19/p13 hg19-delly/1.0 hg19-mappability-map/1.0",
      "referenceGenome": "$HG19_ROOT/hg19_random.fa",
      "dellyExclude": "$HG19_DELLY_ROOT/human.hg19.excl.tsv",
      "genomeMap": "$HG19_MAPPABILITY_MAP_ROOT/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz"
    },
    "hg38": {
      "modules": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg38/p12 hg38-delly/1.0 hg38-mappability-map/1.0",
      "referenceGenome": "$HG38_ROOT/hg38_random.fa",
      "dellyExclude": "$HG38_DELLY_ROOT/human.hg38.excl.tsv",
      "genomeMap": "$HG38_MAPPABILITY_MAP_ROOT/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
    }
  }   


  scatter(s in inputSamples) {
    String inputBams = "~{s.bam}"
    String inputBamIndexes = "~{s.bai}"
    String inputVcfs = "~{s.vcf}"
  }


  #Merge SV sites for each sample
  call mergeSVSites {
    input:
      inputVcfs = inputVcfs,
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }


  scatter(sampleIndex in range(length(inputBams))) {
    String sampleName = basename(inputBams[sampleIndex], ".bam")

    #Genotype SVs for each sample
    call svGenotype {
      input:
        mergedBcf = mergeSVSites.mergedBcf,
        inputBam = inputBams[sampleIndex],
        inputBamIndex = inputBamIndexes[sampleIndex],
        outputFileNamePrefix = sampleName,
        modules = resources[reference].modules,
        referenceGenome = resources[reference].referenceGenome,
        dellyExclude = resources[reference].dellyExclude
    }
  }


  call mergeSamples as mergeSVSamples {
    input:
      genotypedBcfFiles = svGenotype.genoBcf,
      genotypedBcfIndexes = svGenotype.genoBcfIndex, 
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }


  call svFilter {
    input:
      mergedBcf = mergeSVSamples.mergedBcf,
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }


  scatter(sampleIndex in range(length(inputBams))) {
    #Call CNVs for each sample
    call callCNV {
      input:
        inputBam = inputBams[sampleIndex],
        inputBamIndex = inputBamIndexes[sampleIndex],
        germlineSVBcf = svFilter.germlineSVBcf,
        outputFileNamePrefix = basename(inputBams[sampleIndex], ".bam"),
        modules = resources[reference].modules,
        referenceGenome = resources[reference].referenceGenome,
        genomeMap = resources[reference].genomeMap      
    }
  }


  #Merge CNVs into unified site list
  call mergeCNVSites {
    input:
      cnvBcfs = callCNV.cnvBcf,
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }


  scatter(sampleIndex in range(length(inputBams))) {
    #Genotype CNVs for each sample
    call cnvGenotype {
      input:
        mergedCNVBcf = mergeCNVSites.mergedCNVBcf,
        inputBam = inputBams[sampleIndex],
        inputBamIndex = inputBamIndexes[sampleIndex],
        outputFileNamePrefix = basename(inputBams[sampleIndex], ".bam"),
        modules = resources[reference].modules,
        referenceGenome = resources[reference].referenceGenome,
        genomeMap = resources[reference].genomeMap
    }
  }


  call mergeSamples as mergeCNVSamples {
    input:
      genotypedBcfFiles = cnvGenotype.genoCNVBcf,
      genotypedBcfIndexes = cnvGenotype.genoCNVBcfIndex, 
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }


  #Filter for germline CNVs
  call cnvFilter {
    input:
      mergedGenoBcf = mergeCNVSamples.mergedBcf,
      modules = resources[reference].modules,
      outputFileNamePrefix = outputFileNamePrefix
  }



  meta {
    author: "Hannah Driver, Muna Mohamed"
    email: "hdriver@oicr.on.ca, mmohamed@oicr.on.ca"
    description: "Joint germline structural variant and copy number variant calling using Delly"
    dependencies: [
      {
        name: "delly/0.9.1",
        url: "https://github.com/dellytools/delly/releases/download/v0.9.1/delly_v0.9.1_linux_x86_64bit"
      },
      {
        name: "bcftools/1.9",
        url: "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2"
      },
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      }
    ]
    output_meta: {
      germlineSVs: "Filtered vcf file containing PASS structural variant calls",
      germlineCNVs: "Filtered vcf file containing copy number variant calls"
    }
  }


  output {
    File germlineSVs = svFilter.germlineSVVcf
    File germlineCNVs = cnvFilter.germlineCNVs
  }

}


task mergeSVSites{
  input {
    Array[File] inputVcfs
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8 
    Int timeout = 4
  }

  parameter_meta {
    inputVcfs: "Array of vcf files"
    modules: "Modules needed to run Delly"
    outputFileNamePrefix: "Prefix for output file naming"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly merge -o ~{outputFileNamePrefix}.merged.sites.bcf ~{sep = " " inputVcfs}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedBcf = "~{outputFileNamePrefix}.merged.sites.bcf"
  }
}


task svGenotype{
  input {
    File mergedBcf
    File inputBam
    File inputBamIndex
    String outputFileNamePrefix
    String modules
    String referenceGenome
    String dellyExclude
    Int jobMemory = 30 
    Int timeout = 24
  }

  parameter_meta {
    mergedBcf: "Merged bcf file from >=20 unrelated samples"
    inputBam: "Bam file to be genotyped"
    inputBamIndex: "Index file for bam that will be genotyped"
    outputFileNamePrefix: "Prefix for output file naming"
    modules: "Modules needed to run Delly"
    referenceGenome: "Path to fasta file with genomic assembly"
    dellyExclude: "List of regions to exclude (telomeres and centromeres)"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly call -g ~{referenceGenome} -v ~{mergedBcf} -o ~{outputFileNamePrefix}.geno.bcf -x ~{dellyExclude} ~{inputBam}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

    output {
    File genoBcf = "~{outputFileNamePrefix}.geno.bcf"
    File genoBcfIndex = "~{outputFileNamePrefix}.geno.bcf.csi"
  }
}


task mergeSamples{
  input { 
    Array[File] genotypedBcfFiles
    Array[File] genotypedBcfIndexes
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8 
    Int timeout = 4
  }

  parameter_meta {
    genotypedBcfFiles: "Array of bcf files that have gone through Delly's genotype filter"
    genotypedBcfIndexes: "Array of index files for bcfs that have one through Delly's genotype filter"
    modules: "Modules needed to run Delly"
    outputFileNamePrefix: "Prefix for output file naming"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    #Merge all genotyped samples to get a single vcf/bcf using bcftools merge
    bcftools merge -m id -O b -o ~{outputFileNamePrefix}.merged.geno.bcf ~{sep = " " genotypedBcfFiles}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedBcf = "~{outputFileNamePrefix}.merged.geno.bcf"
  }
}


task svFilter {
  input {
    File mergedBcf
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8 
    Int timeout = 4
  }

  parameter_meta {
    mergedBcf: "Merged bcf file from >=20 unrelated samples that have gone through Delly's genotype filter"
    modules: "Modules needed to run Delly"
    outputFileNamePrefix: "Prefix for output file naming"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    #Index
    bcftools index ~{mergedBcf}
    #Merge
    delly filter -f germline -o ~{outputFileNamePrefix}.germlineSV.bcf ~{mergedBcf}
    #Convert to vcf
    bcftools view -O v -o ~{outputFileNamePrefix}.germlineSV.vcf ~{outputFileNamePrefix}.germlineSV.bcf
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File germlineSVVcf = "~{outputFileNamePrefix}.germlineSV.vcf"
    File germlineSVBcf = "~{outputFileNamePrefix}.germlineSV.bcf"
  }
}


task callCNV {
  input {
    File inputBam
    File inputBamIndex
    File germlineSVBcf
    String outputFileNamePrefix
    String modules
    String referenceGenome
    String genomeMap
    Int jobMemory = 30 
    Int timeout = 24
  }

  parameter_meta {
    inputBam: "Bam file that will be analysed for copy number variation"
    inputBamIndex: "Index file for bam that will be analysed for copy number variation"
    germlineSVBcf: "Filtered bcf file containing structural variant calls"
    outputFileNamePrefix: "Prefix for output file naming"
    modules: "Modules needed to run Delly"
    referenceGenome: "Path to fasta file with genomic assembly"
    genomeMap: "Haplotype map file needed to run Delly"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly cnv -o ~{outputFileNamePrefix}.bcf -g ~{referenceGenome} -m ~{genomeMap} -l ~{germlineSVBcf} ~{inputBam}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File cnvBcf = "~{outputFileNamePrefix}.bcf"
  }
}


task mergeCNVSites {
  input {
    Array[File] cnvBcfs
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8 
    Int timeout = 4
  }

  parameter_meta {
    cnvBcfs: "Array of bcf files obtained by CNV calling for all samples with Delly"
    modules: "Modules needed to run Delly"
    outputFileNamePrefix: "Prefix for output file naming"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly merge -e -p -o ~{outputFileNamePrefix}.mergedSites.bcf -m 1000 -n 100000 ~{sep = " " cnvBcfs}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedCNVBcf = "~{outputFileNamePrefix}.mergedSites.bcf"
  }
}


task cnvGenotype {
  input {
    File mergedCNVBcf
    File inputBam
    File inputBamIndex
    String outputFileNamePrefix
    String modules
    String referenceGenome
    String genomeMap
    Int jobMemory = 30 
    Int timeout = 24
  }

  parameter_meta {
    mergedCNVBcf: "Merged bcf file obtained from the CNVs of >=20 unrelated samples"
    inputBam: "Bam file to be genotyped"
    inputBamIndex: "Index file for bam that will be genotyped"
    outputFileNamePrefix: "Prefix for output file naming"
    modules: "Modules needed to run Delly"
    referenceGenome: "Path to fasta file with genomic assembly"
    genomeMap: "Haplotype map file needed to run Delly"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly cnv -u -v ~{mergedCNVBcf} -g ~{referenceGenome} -m ~{genomeMap} -o ~{outputFileNamePrefix}.geno.bcf ~{inputBam}  
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File genoCNVBcf = "~{outputFileNamePrefix}.geno.bcf"
    File genoCNVBcfIndex = "~{outputFileNamePrefix}.geno.bcf.csi"
  }
}


task cnvFilter {
  input {
    File mergedGenoBcf
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8 
    Int timeout = 4  
  }

  parameter_meta {
    mergedGenoBcf: "Bcf file obtained by merging the CNV genotypes from >=20 unrelated samples"
    modules: "Modules needed to run Delly"
    outputFileNamePrefix: "Prefix for output file naming"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    #Index merged bcf
    bcftools index ~{mergedGenoBcf}
    #Merge
    delly classify -f germline -o ~{outputFileNamePrefix}.germlineCNV.bcf ~{mergedGenoBcf}
    #Convert to vcf
    bcftools view -O v -o ~{outputFileNamePrefix}.germlineCNV.vcf ~{outputFileNamePrefix}.germlineCNV.bcf
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File germlineCNVs = "~{outputFileNamePrefix}.germlineCNV.vcf"
  }
}
