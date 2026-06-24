# Биоинформатический проект SARS-CoV-2 и туберкулеза

Проект разделен на два organism-specific workflow:

- COVID / SARS-CoV-2: входные файлы находятся в `data/raw/covid_data/`, результаты пишутся в `results/covid_data/`.
- Tuberculosis / Mycobacterium tuberculosis: входные файлы находятся в `data/raw/tuberculosis_data/`, результаты пишутся в `results/tuberculosis_data/`.

Ссылка на Яндекс-диск с исходными сырыми данными:

https://disk.360.yandex.ru/d/Cxt39v_-ygHNlQ

## Команды

Общая справка:

```powershell
python main.py --help
```

## SARS-CoV-2:

Пример входных данных на рассматриваемых автором последовательностях.

Положите файлы из Яндекс-диска, относящиеся к интересующему вас исследованию в папку:

`data/raw/covid_data`

## SARS-CoV-2 workflow по умолчанию:

```powershell
python main.py sars2
python scripts/generate_report.py
```

По умолчанию SARS-CoV-2 workflow читает:

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

Запуск SARS-CoV-2 workflow для `sequenceX`. В качестве примера используется форма записи sequence1:

```powershell
python scripts/genbank_to_excel.py sars2 --input data/raw/covid_data/sequence1.gb --fasta data/raw/covid_data/sequence1.fasta --output results/covid_data/genbank_table_sequence1.xlsx --nextclade-output results/covid_data/nextclade_sequence1.tsv --nextclade-aligned-fasta results/covid_data/aligned_sequence1.fasta
```

## Workflow для туберкулеза

Пример входных данных на рассматриваемых автором последовательностях.

Из Яндекс-диска нужно будет взять accessions_103.txt, после его положить файл сюда:

`data/raw/tuberculosis_data/accessions_103.txt`

После этого запускаем скачивание из NCBI:

```powershell
python scripts/download_ncbi_assemblies.py --accessions data/raw/tuberculosis_data/accessions_103.txt --output-dir data/raw/tuberculosis_data/ncbi_assemblies
```

Создание базовой таблицы метаданных:

```powershell
python main.py tuberculosis
```

По умолчанию TB workflow читает скачанные NCBI Assembly файлы из:

```text
data/raw/tuberculosis_data/ncbi_assemblies/
```

И пишет:

```text
results/tuberculosis_data/tb_metadata.xlsx
```

Подготовка входных файлов для следующих этапов:

```powershell
python main.py tuberculosis prepare-inputs
```

Команда создает:

```text
results/tuberculosis_data/tb_input_manifest.tsv
results/tuberculosis_data/tb_input_manifest.xlsx
results/tuberculosis_data/tb_input_validation.md
results/tuberculosis_data/inputs/fasta/
data/raw/tuberculosis_data/reference/NC_000962.3.fasta
data/raw/tuberculosis_data/reference/NC_000962.3.gb
```

`tb_input_manifest.tsv` - центральная таблица для последующей автоматизации Snippy, FastTree, анализа мутаций и iTOL. В ней хранятся стабильные sample ID, пути к FASTA/GenBank/assembly-report файлам, метаданные, проверки длины/N/GC и предупреждения.

## Snippy

Создание плана запуска Snippy:

```powershell
python main.py tuberculosis snippy-plan
```

Эта команда не запускает Snippy. Она создает файлы команд для Docker/Linux/WSL на основе `tb_input_manifest.tsv`:

```text
results/tuberculosis_data/snippy/snippy_multi.tsv
results/tuberculosis_data/snippy/snippy_commands.sh
results/tuberculosis_data/snippy/README_snippy.md
```

В Docker этот этап запускается командой `docker compose run --rm tb-snippy`. В Linux/WSL можно запустить `snippy_commands.sh` напрямую после активации окружения, где доступен `snippy`.

После завершения Snippy выполняются следующие этапы:

```powershell
python main.py tuberculosis snippy-summary
python main.py tuberculosis build-tree
python main.py tuberculosis mutation-summary
python main.py tuberculosis itol-export
python main.py tuberculosis report
```

Эти команды создают:

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
results/tuberculosis_data/itol/itol_regions_strip.txt
results/tuberculosis_data/itol/itol_tb_lineage_strip.txt
results/tuberculosis_data/itol/itol_drug_resistance_symbols.txt
reports/tuberculosis_data/tb_analysis_report.md
```

`build-tree` по умолчанию исключает служебную последовательность Snippy `Reference`, поэтому `tb_fasttree.nwk` представляет 103 анализируемые TB-сборки. Если нужно дерево с референсом, используется флаг `--include-reference`.

TB workflow намеренно завершен без TB-Profiler, потому что в учебном наборе данных нет raw FASTQ reads. Гены лекарственной устойчивости суммируются как candidate mutation loci по Snippy-аннотациям, а не как формальные TB-Profiler resistance calls.

В TB Excel создаются листы:

```text
TB_Metadata_counts
TB_Assembly_metadata
TB_Run_Metadata
```

## Вспомогательные скрипты

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

Docker-образ содержит Python-зависимости проекта, Nextclade для SARS-CoV-2, а также Snippy 4.6.0 и FastTree 2.1.11 для TB workflow.

Сборка образа:

```powershell
docker compose build pipeline
```

Запуск SARS-CoV-2 workflow и COVID-отчета:

```powershell
docker compose up pipeline
```

Подготовка TB metadata, manifest, FASTA inputs и Snippy-плана:

```powershell
docker compose run --rm tb-prepare
```

Запуск Snippy для всех TB-сборок внутри Docker:

```powershell
docker compose run --rm tb-snippy
```

Постобработка уже готовых Snippy-результатов: Snippy summary, FastTree, mutation summary, iTOL и финальный TB-отчет:

```powershell
docker compose run --rm tb-post
```

Полный TB workflow одной командой, включая повторный Snippy-запуск для всех 103 сборок:

```powershell
docker compose run --rm tb-full
```

`tb-full` и `tb-snippy` могут выполняться долго, потому что Snippy анализирует 103 бактериальные сборки. Для повторной генерации отчетов без повторного Snippy-запуска используй `tb-post`.

Контейнер использует:

```text
Nextclade CLI: 3.21.2
Nextclade dataset: nextstrain/sars-cov-2/wuhan-hu-1/orfs
Nextclade dataset tag: 2026-04-21--09-39-50Z
Snippy: 4.6.0
FastTree: 2.1.11
```

## Примечания

- SARS-CoV-2 workflow использует Nextclade и подходит только для COVID/SARS-CoV-2.
- Tuberculosis workflow строится по NCBI Assembly GenBank/FASTA файлам; Nextclade для TB не запускается.
- TB-анализ намеренно assembly-based: он использует RefSeq assembly FASTA/GenBank файлы, а не raw FASTQ reads.
- Raw reads не скачиваются по умолчанию, потому что они могут занимать много гигабайт и резко увеличивают время обработки. Это документированное методическое ограничение, а не ошибка.
- TB workflow использует `tb_input_manifest.tsv`, Snippy `core.aln`/`core.tab`, FastTree Newick output и assembly-based mutation summaries как воспроизводимые входные контракты между этапами.
- Raw-данные и результаты игнорируются Git через `.gitignore`.
- Если Excel-файл открыт в Excel/LibreOffice, Windows может заблокировать перезапись. Закрой файл и повтори команду.
