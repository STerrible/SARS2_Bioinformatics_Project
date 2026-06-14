# SARS-CoV-2 Bioinformatics Project

Возможности проекта: 

- читает GenBank-файл `data/raw/sequence.gb`;
- извлекает метаданные, CDS-аннотации и последовательности;
- запускает локальный Nextclade для FASTA-файла `data/raw/sequence.fasta`;
- сохраняет raw-результат Nextclade в `results/nextclade.tsv`;
- сохраняет выравненные Nextclade-последовательности в `results/aligned.fasta`;
- добавляет в итоговый Excel raw-лист Nextclade и удобные summary-листы.

В ближайшее время планируется добавление филогенетики и улучшенных автоматических отчетных файлов, которые создаются из Excel-файла.

## Структура проекта

Основная точка запуска: `main.py`.

## Требования

Нужен Python 3.12 или совместимая версия Python 3.

Python-библиотеки:

```powershell
pip install pandas openpyxl biopython
```

Также нужен установленный Nextclade. По умолчанию проект ожидает:

```text
C:\Games\Nextclade\nextclade.exe
C:\Games\Nextclade\sars-cov-2
```

(Если я не найду способ более общего оформления в коде, я изменю директории в ближайшее время)

В датасете Nextclade используется референс `MN908947 (Wuhan-Hu-1/2019)`. Его последовательность совпадает с `NC_045512.2` из локального GenBank-файла.

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

Из корня проекта:

```powershell
python main.py
```

После успешного запуска создаются или обновляются:

```text
results/nextclade.tsv
results/aligned.fasta
results/genbank_table.xlsx
```

Если `results/genbank_table.xlsx` открыт в Excel, Windows может заблокировать запись. В этом случае закройте файл Excel и запустите команду снова.

## Параметры запуска

Можно переопределить входные и выходные пути:

```powershell
python main.py `
  --input data/raw/sequence.gb `
  --fasta data/raw/sequence.fasta `
  --output results/genbank_table.xlsx `
  --nextclade-output results/nextclade.tsv `
  --nextclade-aligned-fasta results/aligned.fasta `
  --nextclade-exe C:\Games\Nextclade\nextclade.exe `
  --nextclade-dataset C:\Games\Nextclade\sars-cov-2
```

Параметры:

- `--input`: GenBank-файл.
- `--fasta`: FASTA-файл для Nextclade.
- `--output`: итоговый Excel-файл.
- `--nextclade-output`: TSV-файл raw-результата Nextclade.
- `--nextclade-aligned-fasta`: FASTA-файл с выравненными Nextclade-последовательностями.
- `--nextclade-exe`: путь к `nextclade.exe`.
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
: Техническая информация о запуске: версия Nextclade, пути к входным и выходным файлам, путь к датасету, время запуска и размеры обработанных таблиц.

Все листы Excel автоматически получают ширину столбцов по содержимому.

## Проверка результата

Быстрая проверка без перезаписи основного Excel:

```powershell
python main.py --output C:\tmp\genbank_table_check.xlsx
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

После создания Excel можно сгенерировать Markdown-отчет, который сопоставляет текущие результаты pipeline с контрольными выводами дипломной работы:

```powershell
python scripts/generate_report.py
```

По умолчанию скрипт читает:

```text
results/genbank_table.xlsx
```

И создает:

```text
reports/diploma_comparison.md
```

При необходимости пути можно переопределить:

```powershell
python scripts/generate_report.py `
  --input results/genbank_table.xlsx `
  --output reports/diploma_comparison.md
```

## Частые проблемы

`Permission denied: results\genbank_table.xlsx`
: Итоговый Excel открыт в Excel или другом приложении. Закройте файл и повторите запуск.

`Nextclade executable was not found`
: Проверьте путь к `nextclade.exe` или передайте правильный путь через `--nextclade-exe`.

`Nextclade dataset directory was not found`
: Проверьте путь к датасету или передайте правильный путь через `--nextclade-dataset`.

`Input FASTA file was not found`
: Проверьте наличие `data/raw/sequence.fasta` или передайте путь через `--fasta`.

`No GenBank records found`
: Проверьте, что `data/raw/sequence.gb` существует и содержит записи в формате GenBank.

## Важные ограничения

- Pipeline работает с уже подготовленными локальными файлами.
- Pipeline не загружает данные из NCBI.
- Pipeline не выполняет собственный анализ мутаций.
- Pipeline не выполняет pairwise alignment в коде проекта.
- Pipeline не запускает Pangolin.
- Pipeline не строит филогенетические деревья.
- Основной источник мутаций, clade, lineage и QC - Nextclade.
