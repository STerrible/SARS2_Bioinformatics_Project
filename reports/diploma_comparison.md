# Сравнение автоматизированного анализа с выводами дипломной работы

Отчет сгенерирован: `2026-06-15T02:00:01+05:00`.

## Назначение отчета

Этот отчет сопоставляет результаты текущего автоматизированного pipeline с контрольными выводами исходной дипломной работы. Сравнение носит ориентировочный характер: размер и состав текущей выборки могут отличаться от выборки диплома.

## Источники данных

- Итоговый Excel: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\genbank_table.xlsx`
- Время запуска pipeline по `Run_Metadata`: `2026-06-15T01:59:29+05:00`
- Версия Nextclade: `nextclade 3.21.2`
- Датасет Nextclade: `C:\Games\Nextclade\sars-cov-2`
- Мутации, clade, lineage и QC берутся только из листа `Nextclade_results`.

## Масштаб выборки

| Показатель | Диплом | Текущий pipeline |
| --- | --- | --- |
| Объем данных | 4000 нуклеотидных последовательностей; в разделе филогенетики указаны 8000 последовательностей | 35 образцов |
| Средняя длина генома | 29870 нуклеотидов | см. `Metadata_counts` |
| Процент сходства | 85,99% | не пересчитывается; проект не делает собственное множественное выравнивание |

## Сравнение ключевых выводов

| Тезис диплома | Текущий результат | Комментарий |
| --- | --- | --- |
| Наиболее вариабельный ген: `S` | `ORF1a` (29 наблюдений; тип: amino_acid_substitution) | Текущий результат рассчитан по готовым аминокислотным изменениям Nextclade и зависит от состава выборки. |
| Частые аминокислоты: лейцин, треонин, гистидин | см. таблицу ниже | Сравнение выполнено по листу `Amino_Acid_Changes`. |
| Предковый штамм: `NC_045512 (China, 2019/12)` | не определяется | Pipeline не строит филогенетическое дерево и не выявляет предковый штамм. |
| Наиболее вариабельные страны: Египет, Нидерланды, США | см. таблицу ниже | Текущая выборка меньше и может не содержать эти страны. |

## Топ генов по аминокислотным изменениям

| gene | mutation_type | mutation_observations | sample_count |
| --- | --- | --- | --- |
| ORF1a | amino_acid_substitution | 29 | 20 |
| ORF8 | amino_acid_substitution | 12 | 11 |
| ORF1b | amino_acid_substitution | 11 | 7 |
| S | amino_acid_substitution | 10 | 10 |
| ORF3a | amino_acid_substitution | 8 | 7 |
| ORF1a | amino_acid_deletion | 5 | 1 |
| N | amino_acid_substitution | 3 | 3 |
| E | amino_acid_substitution | 1 | 1 |
| ORF1b | amino_acid_deletion | 1 | 1 |
| ORF1b | amino_acid_insertion | 1 | 1 |
| S | amino_acid_deletion | 1 | 1 |

## Топ мутаций по текущей выборке

| mutation_type | gene | mutation | raw_value | sample_count | mutation_observations | percent_of_samples |
| --- | --- | --- | --- | --- | --- | --- |
| amino_acid_substitution | ORF8 | L84S | ORF8:L84S | 11 | 11 | 31.43 |
| nucleotide_substitution |  | C8782T | C8782T | 11 | 11 | 31.43 |
| nucleotide_substitution |  | T28144C | T28144C | 11 | 11 | 31.43 |
| amino_acid_substitution | ORF1a | L3606F | ORF1a:L3606F | 7 | 7 | 20.0 |
| amino_acid_substitution | ORF3a | G251V | ORF3a:G251V | 6 | 6 | 17.14 |
| nucleotide_substitution |  | G11083T | G11083T | 6 | 6 | 17.14 |
| nucleotide_substitution |  | G26144T | G26144T | 6 | 6 | 17.14 |
| nucleotide_substitution |  | C18060T | C18060T | 4 | 4 | 11.43 |
| amino_acid_substitution | ORF1b | P1427L | ORF1b:P1427L | 3 | 3 | 8.57 |
| amino_acid_substitution | ORF1b | Y1464C | ORF1b:Y1464C | 3 | 3 | 8.57 |
| nucleotide_substitution |  | A17858G | A17858G | 3 | 3 | 8.57 |
| nucleotide_substitution |  | C17373T | C17373T | 3 | 3 | 8.57 |
| nucleotide_substitution |  | C17747T | C17747T | 3 | 3 | 8.57 |
| amino_acid_substitution | N | P344S | N:P344S | 2 | 2 | 5.71 |
| amino_acid_substitution | ORF1a | R3323C | ORF1a:R3323C | 2 | 2 | 5.71 |

## Проверка аминокислот из выводов диплома

| amino_acid_name_ru | present_in_current_data | mutation_observations | sample_count |
| --- | --- | --- | --- |
| лейцин | да | 21 | 18 |
| треонин | да | 4 | 4 |
| гистидин | да | 2 | 2 |

## Частоты аминокислот в текущей выборке

| role | amino_acid | amino_acid_name | amino_acid_name_ru | mutation_observations | sample_count |
| --- | --- | --- | --- | --- | --- |
| reference_amino_acid | L | leucine | лейцин | 21 | 18 |
| reference_amino_acid | P | proline | пролин | 11 | 11 |
| reference_amino_acid | G | glycine | глицин | 8 | 7 |
| reference_amino_acid | S | serine | серин | 5 | 4 |
| reference_amino_acid | T | threonine | треонин | 4 | 4 |
| reference_amino_acid | Y | tyrosine | тирозин | 4 | 4 |
| reference_amino_acid | A | alanine | аланин | 3 | 3 |
| reference_amino_acid | I | isoleucine | изолейцин | 3 | 3 |
| reference_amino_acid | R | arginine | аргинин | 3 | 3 |
| reference_amino_acid | F | phenylalanine | фенилаланин | 3 | 2 |
| reference_amino_acid | D | aspartic acid | аспарагиновая кислота | 2 | 2 |
| reference_amino_acid | H | histidine | гистидин | 2 | 2 |
| reference_amino_acid | M | methionine | метионин | 2 | 2 |
| reference_amino_acid | E | glutamic acid | глутаминовая кислота | 1 | 1 |
| reference_amino_acid | V | valine | валин | 1 | 1 |

## Частоты аминокислот по gene

| gene | role | amino_acid | amino_acid_name | amino_acid_name_ru | mutation_observations | sample_count |
| --- | --- | --- | --- | --- | --- | --- |
| E | reference_amino_acid | L | leucine | лейцин | 1 | 1 |
| E | sample_amino_acid | H | histidine | гистидин | 1 | 1 |
| N | reference_amino_acid | P | proline | пролин | 3 | 3 |
| N | sample_amino_acid | S | serine | серин | 3 | 3 |
| ORF1a | reference_amino_acid | L | leucine | лейцин | 8 | 8 |
| ORF1a | reference_amino_acid | P | proline | пролин | 4 | 4 |
| ORF1a | reference_amino_acid | I | isoleucine | изолейцин | 3 | 3 |
| ORF1a | reference_amino_acid | G | glycine | глицин | 2 | 2 |
| ORF1a | reference_amino_acid | M | methionine | метионин | 2 | 2 |
| ORF1a | reference_amino_acid | R | arginine | аргинин | 2 | 2 |
| ORF1a | reference_amino_acid | S | serine | серин | 2 | 2 |
| ORF1a | reference_amino_acid | T | threonine | треонин | 2 | 2 |
| ORF1a | reference_amino_acid | A | alanine | аланин | 1 | 1 |
| ORF1a | reference_amino_acid | D | aspartic acid | аспарагиновая кислота | 1 | 1 |
| ORF1a | reference_amino_acid | E | glutamic acid | глутаминовая кислота | 1 | 1 |
| ORF1a | reference_amino_acid | F | phenylalanine | фенилаланин | 1 | 1 |
| ORF1a | sample_amino_acid | F | phenylalanine | фенилаланин | 8 | 8 |
| ORF1a | sample_amino_acid | I | isoleucine | изолейцин | 5 | 5 |
| ORF1a | sample_amino_acid | L | leucine | лейцин | 3 | 3 |
| ORF1a | sample_amino_acid | S | serine | серин | 3 | 3 |

## Страны из выводов диплома в текущей выборке

| diploma_country | present_in_current_data | current_sample_count |
| --- | --- | --- |
| Египет | нет | 0 |
| Нидерланды | нет | 0 |
| США | да | 12 |

## Страны с наибольшим числом наблюдений мутаций в текущей выборке

| country | mutation_observations | unique_mutations | country_sample_count |
| --- | --- | --- | --- |
| USA | 116 | 78 | 12 |
| China | 40 | 34 | 9 |
| India | 22 | 22 | 2 |
| South Korea | 15 | 15 | 1 |
| Sweden | 11 | 11 | 1 |
| Australia | 6 | 6 | 1 |
| Brazil | 6 | 6 | 1 |
| Viet Nam | 6 | 4 | 2 |
| Japan | 5 | 5 | 1 |
| Finland | 3 | 3 | 1 |
| Italy | 3 | 3 | 1 |
| Taiwan | 3 | 3 | 2 |
| Nepal | 1 | 1 | 1 |

## QC и распределение clade/lineage

| section | value | count | percent |
| --- | --- | --- | --- |
| samples | total | 35 | 100.0 |
| qc.overallStatus | good | 34 | 97.14 |
| qc.overallStatus | mediocre | 1 | 2.86 |
| clade | 19A | 23 | 65.71 |
| clade | 19B | 11 | 31.43 |
| clade | 20A | 1 | 2.86 |
| Nextclade_pango | B | 23 | 65.71 |
| Nextclade_pango | A | 7 | 20.0 |
| Nextclade_pango | A.1 | 3 | 8.57 |
| Nextclade_pango | B.1 | 1 | 2.86 |
| Nextclade_pango | A.3 | 1 | 2.86 |

## Методические ограничения

- Текущая выборка существенно меньше выборки диплома, поэтому совпадения и расхождения нельзя трактовать как окончательный биологический вывод.
- Pipeline не выполняет собственное множественное выравнивание, pairwise alignment или расчет p-distance.
- Pipeline не загружает данные из NCBI и анализирует только локально подготовленные `sequence.gb` и `sequence.fasta`.
- Мутации и QC берутся из Nextclade; листы отчета только агрегируют уже готовые результаты.
- Филогенетика, определение предкового штамма и Pangolin не входят в текущий pipeline.

## Автоматический вывод

Текущий pipeline успешно автоматизирует табличную часть анализа для 35 образцов: извлекает GenBank-метаданные, запускает Nextclade, сохраняет raw-результаты и строит сводные листы. По сравнению с дипломной работой наиболее сильная сторона проекта - воспроизводимость и скорость обновления отчета. Главное ограничение - меньший размер выборки и отсутствие филогенетического блока.
