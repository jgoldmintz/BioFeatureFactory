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
python3 BioFeatureFactory/dependencies/exon_aware_mapping.py \
    --mutations BioFeatureFactory/tests/data/ABCB1_mutations.csv \
    --annotation /path/to/annotation.gtf \
    --reference /path/to/reference_genome.fna \
    --out-fasta /tmp/output_fastas \
    --out-chromosome-mapping /tmp/chromosome_mappings \
    --out-transcript-mapping /tmp/transcript_mappings
```

Adjust external paths as needed for your environment.
