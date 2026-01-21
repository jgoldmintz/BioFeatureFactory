# Example pipeline workflow

## Files
- `ABCB1_mutations.csv` – single-column CSV (`mutant`) listing 50 SNVs.

**Note**: The mutations will best align on GRCh38 reference genome which can be downloaded [here](
  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz).
## Usage
For the vast majority of the pipelines in this repository properly mapped mutations are required for biologically accurate predictions. Therefore running the `exon_aware_mapping.py` first is highly recommended.
#### 
Example usage:

```bash
   python3 utils/exon_aware_mapping.py \
     --mutations /path/to/mutations/ \
     --annotation /path/to/annotations.gtf \
     --reference /path/to/reference_genome.fa \
     --out-fasta /path/to/output_fastas/ \
     --out-chromosome-mapping /path/to/chromosome_mappings/ \
     --out-genomic-mapping /path/to/genomic_mappings/ \
     --out-transcript-mapping /path/to/transcript_mappings/ \
     --out-aa-mapping /path/to/aa_mappings/ \
     --orf /path/to/orf_fastas/ \
     --force-cds transcript_overrides.csv \
     --verbose
```

Adjust external paths as needed for your environment.
