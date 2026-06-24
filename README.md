# SARS-CoV-2 and Tuberculosis Bioinformatics Project

Проект разделен на organism-specific workflows:

- COVID / SARS-CoV-2: входные файлы в `data/raw/covid_data/`, результаты в `results/covid_data/`.
- Tuberculosis: входные файлы в `data/raw/tuberculosis_data/`, результаты в `results/tuberculosis_data/`.

## Структура

```text
data/
  raw/
    covid_data/
      sequence.fasta
      sequence.gb
      sequence1.fasta
      sequence1.gb
    tuberculosis_data/
      accessions_103.txt
      diploma_35_strains/
      ncbi_assemblies/
results/
  covid_data/
    genbank_table.xlsx
    nextclade.tsv
    aligned.fasta
    phylogenetics/
  tuberculosis_data/
    tb_metadata.xlsx
reports/
  covid_data/
    sars2_analysis_report.md
  tuberculosis_data/
```

## Commands

Общая справка:

```powershell
python main.py --help
```

SARS-CoV-2 pipeline по умолчанию:

```powershell
python main.py sars2
python scripts/generate_report.py
```

По умолчанию SARS-CoV-2 читает:

```text
data/raw/covid_data/sequence.gb
data/raw/covid_data/sequence.fasta
```

И пишет:

```text
results/covid_data/genbank_table.xlsx
results/covid_data/nextclade.tsv
results/covid_data/aligned.fasta
reports/covid_data/sars2_analysis_report.md
```

SARS-CoV-2 запуск для `sequence1`:

```powershell
python scripts/genbank_to_excel.py sars2 --input data/raw/covid_data/sequence1.gb --fasta data/raw/covid_data/sequence1.fasta --output results/covid_data/genbank_table_sequence1.xlsx --nextclade-output results/covid_data/nextclade_sequence1.tsv --nextclade-aligned-fasta results/covid_data/aligned_sequence1.fasta
```

Tuberculosis metadata workbook:

```powershell
python main.py tuberculosis
```

По умолчанию TB workflow читает скачанные NCBI assembly-файлы из:

```text
data/raw/tuberculosis_data/ncbi_assemblies/
```

И пишет:

```text
results/tuberculosis_data/tb_metadata.xlsx
```

Tuberculosis input preparation for Snippy/TB-Profiler/FastTree stages:

```powershell
python main.py tuberculosis prepare-inputs
```

This command prepares the next analysis layer:

```text
results/tuberculosis_data/tb_input_manifest.tsv
results/tuberculosis_data/tb_input_manifest.xlsx
results/tuberculosis_data/tb_input_validation.md
results/tuberculosis_data/inputs/fasta/
data/raw/tuberculosis_data/reference/NC_000962.3.fasta
data/raw/tuberculosis_data/reference/NC_000962.3.gb
```

`tb_input_manifest.tsv` is the central table for later Snippy/TB-Profiler/FastTree automation. It contains stable sample IDs, paths to FASTA/GenBank/assembly-report files, metadata, length/N/GC checks, and warnings.

Tuberculosis Snippy plan generation:

```powershell
python main.py tuberculosis snippy-plan
```

This command does not run Snippy. It creates WSL-ready command files from `tb_input_manifest.tsv`:

```text
results/tuberculosis_data/snippy/snippy_multi.tsv
results/tuberculosis_data/snippy/snippy_commands.sh
results/tuberculosis_data/snippy/README_snippy.md
```

Run the generated `snippy_commands.sh` later from WSL after activating the `snippy` conda environment.

After Snippy has finished, validate the outputs and build the FastTree Newick tree:

```powershell
python main.py tuberculosis snippy-summary
python main.py tuberculosis build-tree
python main.py tuberculosis mutation-summary
python main.py tuberculosis itol-export
python main.py tuberculosis report
```

These commands create:

```text
results/tuberculosis_data/snippy/snippy_summary.xlsx
results/tuberculosis_data/snippy/snippy_summary.md
results/tuberculosis_data/snippy/core.no_reference.aln
results/tuberculosis_data/phylogenetics/tb_fasttree.nwk
results/tuberculosis_data/phylogenetics/fasttree.log
results/tuberculosis_data/variants/tb_mutation_summary.xlsx
results/tuberculosis_data/variants/tb_mutation_summary.md
results/tuberculosis_data/itol/itol_country_strip.txt
results/tuberculosis_data/itol/itol_variant_count_bars.txt
results/tuberculosis_data/itol/itol_target_gene_heatmap.txt
reports/tuberculosis_data/tb_analysis_report.md
```

`build-tree` excludes Snippy's `Reference` sequence by default, so `tb_fasttree.nwk` represents the 103 analyzed TB assemblies. Use `--include-reference` if a reference-containing tree is needed.

The TB workflow is intentionally completed without TB-Profiler because raw FASTQ reads are not part of this training dataset. Drug-resistance genes are summarized as candidate mutation loci from Snippy annotations, not as formal TB-Profiler resistance calls.

В TB Excel создаются листы:

```text
TB_Metadata_counts
TB_Assembly_metadata
TB_Run_Metadata
```

## Helper Scripts

```powershell
python scripts/genbank_to_excel.py sars2
python scripts/genbank_to_excel.py tuberculosis
python scripts/generate_report.py
python scripts/create_alignment_short.py
python scripts/prepare_phylogenetic_network_inputs.py
python scripts/create_popart_inputs.py
python scripts/extract_fasta_accessions.py
python scripts/download_ncbi_assemblies.py --metadata-only
python main.py tuberculosis prepare-inputs
python main.py tuberculosis snippy-plan
python main.py tuberculosis snippy-summary
python main.py tuberculosis build-tree
python main.py tuberculosis mutation-summary
python main.py tuberculosis itol-export
python main.py tuberculosis report
```

Извлечение accession-номеров из COVID FASTA:

```powershell
python scripts/extract_fasta_accessions.py --input data/raw/covid_data/sequence.fasta --output data/raw/covid_data/accessions_from_fasta.txt --strip-version
```

## Docker

Docker Compose запускает SARS-CoV-2 workflow и COVID-отчет:

```powershell
docker compose up --build pipeline
```

Контейнер использует:

```text
Nextclade CLI: 3.21.2
Nextclade dataset: nextstrain/sars-cov-2/wuhan-hu-1/orfs
Nextclade dataset tag: 2026-04-21--09-39-50Z
```

## Notes

- SARS-CoV-2 workflow использует Nextclade и подходит только для COVID/SARS-CoV-2.
- Tuberculosis workflow сейчас делает metadata/counts workbook по NCBI Assembly GenBank-файлам; Nextclade для TB не запускается.
- TB analysis is intentionally assembly-based for the training project: it uses RefSeq assembly FASTA/GenBank files, not raw FASTQ reads.
- Raw reads are not downloaded by default because they can require many gigabytes of storage and substantially longer processing time. This is a documented methodological limitation, not a bug.
- The TB workflow now uses `tb_input_manifest.tsv`, Snippy `core.aln`/`core.tab`, FastTree Newick output, and assembly-based mutation summaries as its reproducible input contracts.
- Raw-данные и результаты игнорируются Git через `.gitignore`.
- Если Excel-файл открыт в Excel/LibreOffice, Windows может заблокировать перезапись. Закрой файл и повтори команду.
