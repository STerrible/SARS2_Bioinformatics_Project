# SARS-CoV-2 Bioinformatics Project

Возможности проекта: 

- читает GenBank-файл `data/raw/sequence.gb`;
- извлекает метаданные, CDS-аннотации и последовательности;
- запускает Nextclade CLI для FASTA-файла `data/raw/sequence.fasta`;
- сохраняет raw-результат Nextclade в `results/nextclade.tsv`;
- сохраняет выравненные Nextclade-последовательности в `results/aligned.fasta`;
- добавляет в итоговый Excel raw-лист Nextclade и удобные summary-листы.

В ближайшее время планируется добавление филогенетики и улучшенных автоматических отчетных файлов, которые создаются из Excel-файла.

## Структура проекта

Основная точка запуска: `main.py`.

Код разделен на общий пакет `bioseq_pipeline` и organism-specific профили:

- `src/bioseq_pipeline/core`: общие GenBank/FASTA/QC/Excel-функции.
- `src/bioseq_pipeline/tools`: адаптеры внешних инструментов, сейчас Nextclade.
- `src/bioseq_pipeline/organisms/sars2`: текущий SARS-CoV-2 workflow.
- `src/bioseq_pipeline/organisms/tuberculosis`: заготовка под будущий TB workflow.

CLI поддерживает явный выбор профиля: `python main.py sars2`. Для совместимости `python main.py` без профиля также запускает SARS-CoV-2 pipeline.

## Команды CLI

Общая справка по доступным профилям:

```powershell
python main.py --help
```

Запуск текущего SARS-CoV-2 pipeline:

```powershell
python main.py sars2
```

Справка по параметрам SARS-CoV-2 pipeline:

```powershell
python main.py sars2 --help
```

Совместимый старый запуск без явного профиля:

```powershell
python main.py
```

Алиасы SARS-CoV-2 pipeline:

```powershell
python main.py sars-cov-2
python main.py covid
```

Заготовка под будущий pipeline для туберкулеза:

```powershell
python main.py tuberculosis
python main.py tb
```

Вспомогательные скрипты:

```powershell
python scripts/genbank_to_excel.py sars2
python scripts/generate_report.py
python scripts/create_alignment_short.py
python scripts/prepare_phylogenetic_network_inputs.py
python scripts/create_popart_inputs.py
python scripts/extract_fasta_accessions.py
```

Извлечение accession-номеров из FASTA для последующего ручного скачивания данных из NCBI:

```powershell
python scripts/extract_fasta_accessions.py `
  --input data/raw/sequence.fasta `
  --output data/raw/accessions_from_fasta.txt `
  --strip-version
```

Если `--output` не указан, accession-номера печатаются в консоль.

## Требования

Рекомендуемый способ запуска - Docker. Он фиксирует версию Nextclade CLI и версию SARS-CoV-2 dataset внутри образа, поэтому результат не зависит от локальной установки `nextclade.exe`.

Зафиксированные версии Docker-сборки:

```text
Nextclade CLI: 3.21.2
Nextclade dataset: nextstrain/sars-cov-2/wuhan-hu-1/orfs
Nextclade dataset tag: 2026-04-21--09-39-50Z
```

В датасете Nextclade используется референс `MN908947 (Wuhan-Hu-1/2019)`. Его последовательность совпадает с `NC_045512.2` из локального GenBank-файла.

Для локального запуска без Docker нужен Python 3.12 или совместимая версия Python 3, Python-библиотеки из `requirements.txt` и установленный Nextclade CLI:

```powershell
pip install -r requirements.txt
```

Пути к локальному Nextclade можно передать через аргументы CLI или переменные окружения:

```powershell
$env:NEXTCLADE_EXE = "C:\Games\Nextclade\nextclade.exe"
$env:NEXTCLADE_DATASET = "C:\Games\Nextclade\sars-cov-2"
```

## Входные файлы

По умолчанию используются:

```text
data/raw/sequence.gb
data/raw/sequence.fasta
```

`sequence.gb` нужен для листов с метаданными, CDS и последовательностями.

`sequence.fasta` нужен для запуска Nextclade.

Файлы должны уже лежать в проекте. Автоматической загрузки из NCBI в pipeline нет. Возможно, добавлю в будущем.

## Запуск

Воспроизводимый запуск через Docker Compose:

```powershell
docker compose up --build pipeline
```

Команда собирает образ, скачивает зафиксированный dataset tag внутрь образа, запускает `python main.py sars2` и затем обновляет `reports/sars2_analysis_report.md`.

После успешного запуска создаются или обновляются:

```text
results/nextclade.tsv
results/aligned.fasta
results/genbank_table.xlsx
reports/sars2_analysis_report.md
```

Локальный запуск без Docker:

Из корня проекта:

```powershell
python main.py sars2
python scripts/generate_report.py
```

Если `results/genbank_table.xlsx` открыт в Excel, Windows может заблокировать запись. В этом случае закройте файл Excel и запустите команду снова.

## Параметры запуска

Можно переопределить входные и выходные пути:

```powershell
python main.py sars2 `
  --input data/raw/sequence.gb `
  --fasta data/raw/sequence.fasta `
  --output results/genbank_table.xlsx `
  --nextclade-output results/nextclade.tsv `
  --nextclade-aligned-fasta results/aligned.fasta `
  --nextclade-exe nextclade `
  --nextclade-dataset /opt/nextclade/datasets/sars-cov-2
```

Параметры:

- `--input`: GenBank-файл.
- `--fasta`: FASTA-файл для Nextclade.
- `--output`: итоговый Excel-файл.
- `--nextclade-output`: TSV-файл raw-результата Nextclade.
- `--nextclade-aligned-fasta`: FASTA-файл с выравненными Nextclade-последовательностями.
- `--nextclade-exe`: путь к `nextclade.exe` или имя команды из `PATH`, например `nextclade`.
- `--nextclade-dataset`: путь к локальному датасету Nextclade.

## Листы итогового Excel

`Metadata_counts`
: Метаданные по каждому GenBank record: accession, описание, организм, страна, дата сбора, host, isolate, длина последовательности и состав нуклеотидов.

`CDS_features`
: CDS-аннотации из GenBank: gene, product, protein_id, координаты, strand, location и признак наличия translation.

`Sequences`
: Последовательности из GenBank.

`QC_summary`
: Базовая сводка по локальному GenBank-файлу: число records, число CDS и число records с `N`.

`Nextclade_results`
: Raw-таблица Nextclade без сокращения колонок. Значения `N/A` сохраняются как текст, а не превращаются в пустые значения.

`Nextclade_QC`
: Удобная QC-таблица по каждому образцу на основе колонок Nextclade: accession, accession version, страна, дата сбора, общий QC-статус, score, coverage, missing data, mixed sites, private mutations, SNP clusters, frameshifts, stop codons, warnings и errors.

`Nextclade_Mutations`
: Нормализованный список мутаций из готовых колонок Nextclade. Для каждой строки добавлены accession, accession version, страна и дата сбора:

- `substitutions`
- `deletions`
- `insertions`
- `frameShifts`
- `aaSubstitutions`
- `aaDeletions`
- `aaInsertions`

Проект не пересчитывает эти мутации самостоятельно. Лист только раскладывает уже готовые значения Nextclade по строкам.

`Nextclade_Summary`
: Общая статистика по Nextclade: количество образцов, распределение `qc.overallStatus`, `clade` и `Nextclade_pango`.

`Nextclade_Gene_Summary`
: Сводка аминокислотных замен, делеций и инсерций по генам/белкам. Лист показывает, на какие гены приходится больше всего изменений по данным Nextclade.

`Nextclade_Top_Mutations`
: Самые частые нуклеотидные и аминокислотные мутации из готовых колонок Nextclade. Для каждой мутации указано число образцов и процент от общего числа образцов.

`Country_Summary`
: Сводка по странам: число образцов, распределение QC-статусов, clade и `Nextclade_pango`.

`Country_Mutations`
: Мутации по странам. Лист показывает, какие готовые Nextclade-мутации встречаются в образцах из каждой страны и у какой доли образцов страны они обнаружены.

`Amino_Acid_Changes`
: Частоты аминокислот, участвующих в аминокислотных заменах. Этот лист помогает проверять тезисы о наиболее часто затрагиваемых аминокислотах, например лейцине, треонине и гистидине.

`Amino_Acid_Changes_By_Gene`
: Частоты аминокислотных замен отдельно по каждому gene. Лист помогает увидеть, какие аминокислоты чаще затрагиваются в конкретных белках.

`aligned_fasta`
: Выравненные последовательности из FASTA-файла, созданного Nextclade. Лист содержит `seq_id`, `description` и `aligned_sequence`; это наглядное представление результата Nextclade, а не собственное pairwise alignment в коде проекта.

`Run_Metadata`
: Техническая информация о запуске: версия Nextclade, Docker image, requested/actual dataset tag, пути к входным и выходным файлам, путь к датасету, время запуска и размеры обработанных таблиц.

Все листы Excel автоматически получают ширину столбцов по содержимому.

## Проверка результата

Быстрая проверка без перезаписи основного Excel:

```powershell
python main.py sars2 --output C:\tmp\genbank_table_check.xlsx
```

Ожидаемый набор листов:

```text
Metadata_counts
CDS_features
Sequences
QC_summary
Nextclade_results
Nextclade_QC
Nextclade_Mutations
Nextclade_Summary
Nextclade_Gene_Summary
Nextclade_Top_Mutations
Country_Summary
Country_Mutations
Amino_Acid_Changes
Amino_Acid_Changes_By_Gene
aligned_fasta
Run_Metadata
```

## Markdown-отчет

После создания Excel можно сгенерировать Markdown-отчет с краткой сводкой, QC, географией выборки, clade/lineage и основными мутациями:

```powershell
python scripts/generate_report.py
```

По умолчанию скрипт читает:

```text
results/genbank_table.xlsx
```

И создает:

```text
reports/sars2_analysis_report.md
```

При необходимости пути можно переопределить:

```powershell
python scripts/generate_report.py `
  --input results/genbank_table.xlsx `
  --output reports/sars2_analysis_report.md
```

## Частые проблемы

`Permission denied: results\genbank_table.xlsx`
: Итоговый Excel открыт в Excel или другом приложении. Закройте файл и повторите запуск.

`Nextclade executable was not found`
: При Docker-запуске пересоберите образ через `docker compose up --build pipeline`. При локальном запуске проверьте путь к `nextclade.exe`, переменную `NEXTCLADE_EXE` или передайте правильный путь через `--nextclade-exe`.

`Nextclade dataset directory was not found`
: При Docker-запуске пересоберите образ, чтобы dataset был скачан внутрь контейнера. При локальном запуске проверьте путь к датасету, переменную `NEXTCLADE_DATASET` или передайте правильный путь через `--nextclade-dataset`.

`Input FASTA file was not found`
: Проверьте наличие `data/raw/sequence.fasta` или передайте путь через `--fasta`.

`No GenBank records found`
: Проверьте, что `data/raw/sequence.gb` существует и содержит записи в формате GenBank.

## Важные ограничения

- Pipeline работает с уже подготовленными локальными файлами.
- Pipeline не загружает данные из NCBI.
- Docker-сборка загружает только зафиксированный Nextclade dataset tag, а не новые GenBank-данные.
- Pipeline не выполняет собственный анализ мутаций.
- Pipeline не выполняет pairwise alignment в коде проекта.
- Pipeline не запускает Pangolin.
- Pipeline не строит филогенетические деревья.
- Основной источник мутаций, clade, lineage и QC - Nextclade.
