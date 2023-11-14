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
}

workflow dellyGermline {
  input {
    String reference
    Array[inputSamples] inputSamples
  }

  parameter_meta {
    reference: "The genome reference build. for example: hg19, hg38"
    inputSamples: "Collection of BAM, BAI, and VCF files for >= 20 samples"
  }

  Map[String,dellyResources] resources = {
    "hg19": {
      "modules": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg19/p13 hg19-delly/1.0",
      "referenceGenome": "$HG19_ROOT/hg19_random.fa",
      "dellyExclude": "$HG19_DELLY_ROOT/human.hg19.excl.tsv"
    },
    "hg38": {
      "modules": "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg38/p12 hg38-delly/1.0",
      "referenceGenome": "$HG38_ROOT/hg38_random.fa",
      "dellyExclude": "$HG38_DELLY_ROOT/human.hg38.excl.tsv"
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
      modules = resources[reference].modules
  }


  scatter(sampleIndex in range(length(inputBams))) {
    String sampleName = basename(inputBams[sampleIndex], ".bam")

    # Genotype SVs for each sample
    call genotype {
      input:
        mergedVcf = mergeSVSites.mergedVcf,
        inputBam = inputBams[sampleIndex],
        inputBamIndex = inputBamIndexes[sampleIndex],
        outputFileNamePrefix = sampleName,
        modules = resources[reference].modules,
        referenceGenome = resources[reference].referenceGenome,
        dellyExclude = resources[reference].dellyExclude
    }
  }


  call mergeSamples {
    input:
      genotypedVcfFiles = genotype.genoVcf,
      genotypedVcfIndexes = genotype.genoVcfIndex, 
      modules = resources[reference].modules
  }


  call filter {
    input:
      mergedBcf = mergeSamples.mergedBcf,
      modules = resources[reference].modules
  }



  meta {
    author: "Hannah Driver"
    email: "hannah.driver@oicr.on.ca"
    description: "Joint germline structural variant calling using Delly"
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
      germlineVcf: "filtered vcf file containing PASS structural variant calls"
    }
  }


  output {
    File germlineVcf = filter.germlineVcf
  }

}


task mergeSVSites{
  input {
    Array[File] inputVcfs
    String modules
    Int jobMemory = 24 
    Int timeout = 24
  }

  parameter_meta {
    inputVcfs: "Array of vcf files"
    modules: "modules needed to run Delly"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly merge -o merged.sites.vcf ~{sep = " " inputVcfs}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "merged.sites.vcf"
  }
}

task genotype{
  input {
    File mergedVcf
    File inputBam
    File inputBamIndex
    String outputFileNamePrefix
    String modules
    String referenceGenome
    String dellyExclude
    Int jobMemory = 24 
    Int timeout = 24
  }

  parameter_meta {
    mergedVcf: "Merged vcf file from >=20 unrelated samples"
    inputBam: "bam file to be genotyped"
    inputBamIndex: "index file for bam that will be genotyped"
    outputFileNamePrefix: "sample ID, this is provided to maivs and cannot include reseerved characters [;,_\\s] "
    modules: "modules needed to run Delly"
    referenceGenome: "Path to fasta file with genomic assembly"
    dellyExclude: "List of regions to exclude (telomeres and centromeres)"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    delly call -g ~{referenceGenome} -v ~{mergedVcf} -o ~{outputFileNamePrefix}.geno.vcf -x ~{dellyExclude} ~{inputBam}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

    output {
    File genoVcf = "~{outputFileNamePrefix}.geno.vcf"
    File genoVcfIndex = "~{outputFileNamePrefix}.geno.vcf.csi"
  }
}


task mergeSamples{
  input { 
    Array[File] genotypedVcfFiles
    Array[File] genotypedVcfIndexes
    String modules
    Int jobMemory = 24 
    Int timeout = 24
  }

  parameter_meta {
    genotypedVcfFiles: "Array of vcf files that have gone through Delly's genotype filter"
    genotypedVcfIndexes: "Array of index files for vcfs that have one through Delly's genotype filter"
    modules: "modules needed to run Delly"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    #Merge all genotyped samples to get a single VCF/BCF using bcftools merge
    bcftools merge -m id -O b -o merged.geno.bcf ~{sep = " " genotypedVcfFiles}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "merged.geno.bcf"
  }
}


task filter{
  input {
    File mergedBcf
    String modules
    Int jobMemory = 24 
    Int timeout = 24
  }

  parameter_meta {
    mergedBcf: "Merged vcf file from >=20 unrelated samples that have gone through Delly's genotype filter"
    modules: "modules needed to run Delly"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -eu -o pipefail
    #Index
    bcftools index ~{mergedBcf}
    #Merge
    delly filter -f germline -o germline.vcf ~{mergedBcf}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File germlineVcf = "germline.vcf"
  }
}
