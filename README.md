# dellyGermline

Joint germline structural variant calling using Delly

## Overview

## Dependencies

* [delly 0.9.1](https://github.com/dellytools/delly/releases/download/v0.9.1/delly_v0.9.1_linux_x86_64bit)
* [bcftools 1.9](https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2)
* [tabix 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2)


## Usage

### Cromwell
```
java -jar cromwell.jar run dellyGermline.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`reference`|String|The genome reference build. for example: hg19, hg38
`inputSamples`|Array[inputSamples]|Collection of BAM, BAI, and BCF files for >= 20 samples


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`mergeSVSites.jobMemory`|Int|24|Memory allocated for this job
`mergeSVSites.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`genotype.jobMemory`|Int|24|Memory allocated for this job
`genotype.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`mergeSamples.jobMemory`|Int|24|Memory allocated for this job
`mergeSamples.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`filter.jobMemory`|Int|24|Memory allocated for this job
`filter.timeout`|Int|24|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description
---|---|---
`germlineBcf`|File|filtered bcf file containing PASS structural variant calls


## Commands
 This section lists command(s) run by the dellyGermline workflow
 
 * Running dellyGermline
 
 === Merge SV Sites Across Samples ===.
 
 ```
     set -eu -o pipefail
     delly merge -o MERGED_SITES_BCF INPUT_BCFS
   ```
 
 === Genotype Samples Across Merged SV Sites ===.
 
 ```
     set -eu -o pipefail
     delly call -g REFERENCE_GENOME -v MERGED_SITES_BCF -o GENOTYPED_BCF -x DELLY_EXCLUDE INPUT_BAMS
   ```
 
 === Merge Genotyped Samples into Single BCF ===.
 
 ```
     set -eu -o pipefail
     bcftools merge -m id -O b -o MERGED_GENOTYPED_BCF GENOTYPED_BCFS
   ```
 
 === Apply Germline Filter to Merged BCF ===.
 
 ```
     set -eu -o pipefail
     #Index
     bcftools index MERGED_GENOTYPED_BCF
     #Merge
     delly filter -f germline -o GERMLINE_BCF MERGED_GENOTYPED_BCF
   ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
