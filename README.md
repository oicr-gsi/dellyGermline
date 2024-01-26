# dellyGermline

Joint germline structural variant and copy number variant calling using Delly

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
`inputSamples`|Array[inputSamples]|Collection of BAM, BAI, and VCF files for >= 20 samples


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`mergeSVSites.jobMemory`|Int|12|Memory allocated for this job
`mergeSVSites.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`genotype.jobMemory`|Int|12|Memory allocated for this job
`genotype.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`mergeSVSamples.jobMemory`|Int|12|Memory allocated for this job
`mergeSVSamples.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`filter.jobMemory`|Int|12|Memory allocated for this job
`filter.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`callCNV.jobMemory`|Int|12|Memory allocated for this job
`callCNV.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`mergeCNVSites.jobMemory`|Int|12|Memory allocated for this job
`mergeCNVSites.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`cnvGenotype.jobMemory`|Int|12|Memory allocated for this job
`cnvGenotype.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`mergeCNVSamples.jobMemory`|Int|12|Memory allocated for this job
`mergeCNVSamples.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`cnvFilter.jobMemory`|Int|12|Memory allocated for this job
`cnvFilter.timeout`|Int|12|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description
---|---|---
`germlineSVVcf`|File|Filtered vcf file containing PASS structural variant calls
`germlineCNVs`|File|Filtered vcf file containing copy number variant calls


## Commands
 This section lists command(s) run by the dellyGermline workflow
 
 * Running dellyGermline
 
 === Merge SV Sites Across Samples ===.
 
 '''
     set -eu -o pipefail
     delly merge -o MERGED_SV_SITES_BCF INPUT_VCFS
 '''
 
 === Genotype Samples Across Merged SV Sites ===.
 
 '''
     set -eu -o pipefail
     delly call -g REFERENCE_GENOME -v MERGED_SV_SITES_BCF -o GENOTYPED_BCF -x DELLY_EXCLUDE INPUT_BAM
 '''
 
 === Merge Genotyped Samples into Single BCF ===.
 
 '''
     set -eu -o pipefail
     bcftools merge -m id -O b -o MERGED_GENOTYPED_BCF GENOTYPED_BCFS
 '''
 
 === Apply Germline Filter to Merged BCF with SV Calls and Convert to VCF ===.
 
 '''
     set -eu -o pipefail
     bcftools index MERGED_GENOTYPED_BCF
     delly filter -f germline -o GERMLINE_SV_BCF MERGED_GENOTYPED_BCF
     bcftools view -O v -o GERMLINE_SV_VCF GERMLINE_SV_BCF
 '''
 
 === Call CNVs for Samples While Refining Breakpoints Using the SV Calls ===.
 
 '''
     set -eu -o pipefail
     delly cnv -o CNV_BCF -g REFERENCE_GENOME -m MAPPABILITY_MAP -l GERMLINE_SV_BCF INPUT_BAM
 '''
 
 === Merge CNV Sites Across Samples ===.
 
 '''
     set -eu -o pipefail
     delly merge -e -p -o MERGED_SITES_CNV_BCF -m 1000 -n 100000 CNV_BCFS
 '''
 
 === Genotype Samples Across Merged CNV Sites ===.
 
 '''
     set -eu -o pipefail
     delly cnv -u -v MERGED_SITES_CNV_BCF -g REFERENCE_GENOME -m MAPPABILITY_MAP -o GENOTYPED_BCF INPUT_BAM
 '''
 
 === Merge Genotyped CNV Samples into Single BCF ===.
 
 '''
     set -eu -o pipefail
     bcftools merge -m id -O b -o MERGED_GENOTYPED_CNV_BCF GENOTYPED_BCFS
 '''
 
 === Apply Germline Filter to Merged BCF with CNV Calls and Convert to VCF ===.
 
 '''
     set -eu -o pipefail
     bcftools index MERGED_GENOTYPED_CNV_BCF
     delly classify -f germline -o GERMLINE_CNVS_BCF MERGED_GENOTYPED_CNV_BCF
     bcftools view -O v -o GERMLINE_CNVS_VCF GERMLINE_CNVS_BCF
 ''' ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
