# Kraken2 Fungal Classification and Binning Workflow

This pipeline documents the process of building a **fungal-only Kraken2 database**, classifying reads, and assembling/binning fungal MAGs from metagenomic data. Manual fixes were required due to **broken NCBI FTP links** as of early 2024.

- **Kraken2 version:** 2.0.9-beta  
- **Kraken2 GitHub:** [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)

---

## NCBI FTP Path Bug
Kraken2’s `--download-library fungi` command fails with an error like:

```
rsync_from_ncbi.pl: unexpected FTP path (new server?) for https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/... 
```

This is due to a change in NCBI's directory structure. See open issue:  
https://github.com/DerrickWood/kraken2/issues/226

## Solution: Manual Download + Header Reformatting

### Step 1: Download and clean fungal genome list
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt -O clean_assembly_summary.txt
```

### Step 2: Fix missing header (replace with full 44-column format)
```bash
HEADER=$'assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tasm_submitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material\tasm_not_live_date\tassembly_type\tgroup\tgenome_size\tgenome_size_ungapped\tgc_percent\treplicon_count\tscaffold_count\tcontig_count\tannotation_provider\tannotation_name\tannotation_date\ttotal_gene_count\tprotein_coding_gene_count\tnon_coding_gene_count\tpubmed_id\ttranscriptome_assemblies\tannotation_pipeline\tannotation_method\tfeature_count\tassembly_notes\tassembly_type_notes'

grep -v '^#' clean_assembly_summary.txt > no_header.txt
(echo -e "$HEADER"; cat no_header.txt) > clean_assembly_summary_fixed.txt

# Confirm header fix
head -n 1 clean_assembly_summary_fixed.txt | awk -F '\t' '{print NF}'
```

### Step 3: Generate fungal genome links
```bash
awk -F '\t' '$12 == "Complete Genome" && $11 == "latest" && $5 == "reference genome" {print $20}' clean_assembly_summary_fixed.txt | \
    sed 's|ftp://|https://|' | \
    awk -F '/' '{print $0"/"$NF"_genomic.fna.gz"}' > fungal_links.txt

# Sanity check
head fungal_links.txt
wc -l fungal_links.txt
```

### Step 4: Download curated fungal genomes
```bash
wget -i fungal_links.txt
```

### Step 5: Test genome integrity
```bash
# If all files are okay, gunzip -t returns no output or errors
gunzip -t fungal_genomes/*.gz
```

### Step 6: Fix FASTA headers to include kraken:taxid
Kraken2 requires headers like:
```
>kraken_seq_<accession>_<count>|kraken:taxid|<taxid>
```

```bash
SUMMARY="clean_assembly_summary_fixed.txt"
GENOME_DIR="fungal_genomes"
OUT_DIR="kraken_ready_genomes"

mkdir -p "$OUT_DIR"

grep -v '^#' "$SUMMARY" | while IFS=$'\t' read -r accession _ _ _ _ taxid _ _ _ _ _ _ _ _ _ _ _ _ _ ftp_path _; do
    filename=$(basename "$ftp_path")_genomic.fna.gz
    input_file="$GENOME_DIR/$filename"
    out_file="$OUT_DIR/$filename"

    if [ -f "$out_file" ]; then echo "Already processed: $filename"; continue; fi
    if [ ! -f "$input_file" ]; then echo "Missing input: $input_file"; continue; fi

    echo "🔧 Fixing: $filename with taxid $taxid"

    zcat "$input_file" | \
        awk -v taxid="$taxid" -v asm="$accession" '
            /^>/ {header_count += 1; print ">kraken_seq_" asm "_" header_count "|kraken:taxid|" taxid; next}
            { print }
        ' | gzip > "$out_file"
done

echo "All new genomes rewritten with kraken-compatible headers."
```

---

## Build Kraken2 Fungal Database
```bash
#!/bin/bash
#SBATCH --job-name=kraken2_fungal_db_build
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=3-00:00:00
#SBATCH --mem=250gb
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=First.Last@colostate.edu
#SBATCH --nodelist=zenith

DB_DIR="/home/projects/Agribiome/SPUR_green_roofs/MetaG/Kraken/kraken_fungi_db"
READY_GENOMES="/home/projects/Agribiome/SPUR_green_roofs/MetaG/Kraken/kraken_ready_genomes"

for genome in "$READY_GENOMES"/*.fna; do
  kraken2-build --add-to-library "$genome" --db "$DB_DIR"
done

kraken2-build --build --db "$DB_DIR" --threads 20
```

---

## Classify Metagenomic Reads
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=1-00:00:00
#SBATCH --mem=100gb
#SBATCH --nodelist=zenith
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=First.last@colostate.edu

cd ../classification

DB_DIR="/home/projects/Agribiome/SPUR_green_roofs/MetaG/Kraken/kraken_fungi_db"
READDIR="/home/projects/Agribiome/SPUR_green_roofs/MetaG/trimmed_reads"

SAMPLES=("A" "B" "C")
SNUMS=("52" "51" "53")

for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"
  SNUM="${SNUMS[$i]}"
  echo "🔍 Classifying sample GHRZ_${SAMPLE}..."

  R1="${READDIR}/GHRZ_${SAMPLE}_S${SNUM}_R1_BBDuktrimmed.fastq.gz"
  R2="${READDIR}/GHRZ_${SAMPLE}_S${SNUM}_R2_BBDuktrimmed.fastq.gz"

  OUT_PREFIX="GHRZ_${SAMPLE}"
  REPORT="${OUT_PREFIX}_fungal_report.txt"
  OUTPUT="${OUT_PREFIX}_fungal_output.txt"
  CLASSIFIED="${OUT_PREFIX}_classified#"
  UNCLASSIFIED="${OUT_PREFIX}_unclassified#"

  kraken2 \
    --db "$DB_DIR" \
    --threads 20 \
    --paired "$R1" "$R2" \
    --report "$REPORT" \
    --output "$OUTPUT" \
    --classified-out "${CLASSIFIED}.fq" \
    --unclassified-out "${UNCLASSIFIED}.fq"

done

echo "✅ All Kraken2 fungal classifications complete!"
```

---

## ✅ Summary
This workflow provides a reliable workaround for Kraken2 fungal classification when default downloads fail. Manual curation + FASTA header correction restores compatibility. This setup is ideal for microbial ecology or soil metagenome studies targeting the fungal component.

Feel free to fork or cite this workflow for similar projects!
