# Краткий отчет по анализу SARS-CoV-2

Сгенерировано: `2026-06-15T02:12:22+05:00`.
Источник: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\genbank_table.xlsx`.

## Выборка

- Образцов в текущем анализе: **35**.
- В дипломной работе для сравнения фигурируют 4000 нуклеотидных последовательностей и 8000 последовательностей в филогенетическом блоке.
- Средняя длина генома в текущей выборке: **29869.80 н.**
- Средняя длина генома в выводах диплома: **29870 н.**

Самый короткий геном:

| accession_version | country | collection_date | length |
| --- | --- | --- | --- |
| MT163720.1 | USA | 01-Mar-2020 | 29732 |

Самый длинный геном:

| accession_version | country | collection_date | length |
| --- | --- | --- | --- |
| MT121215.1 | China | 02-Feb-2020 | 29945 |

## Средний состав A/G/C/T

| base | mean_count | mean_percent |
| --- | --- | --- |
| A | 8932.0 | 29.9 |
| G | 5859.91 | 19.62 |
| C | 5488.0 | 18.37 |
| T | 9589.86 | 32.11 |

## Clade в текущей выборке

| clade | count | percent |
| --- | --- | --- |
| 19A | 23 | 65.71 |
| 19B | 11 | 31.43 |
| 20A | 1 | 2.86 |

## Страны в текущей выборке

| country | sample_count |
| --- | --- |
| USA | 12 |
| China | 9 |
| India | 2 |
| Taiwan | 2 |
| Viet Nam | 2 |
| Australia | 1 |
| South Korea | 1 |
| Italy | 1 |
| Nepal | 1 |
| Sweden | 1 |
| Brazil | 1 |
| Japan | 1 |
| Finland | 1 |

## Топ генов по числу изменений

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

## Топ аминокислотных замен

| gene | mutation | sample_count | percent_of_samples |
| --- | --- | --- | --- |
| ORF8 | L84S | 11 | 31.43 |
| ORF1a | L3606F | 7 | 20.0 |
| ORF3a | G251V | 6 | 17.14 |
| ORF1b | P1427L | 3 | 8.57 |
| ORF1b | Y1464C | 3 | 8.57 |
| N | P344S | 2 | 5.71 |
| ORF1a | R3323C | 2 | 5.71 |
| ORF1a | T1840I | 2 | 5.71 |
| E | L37H | 1 | 2.86 |
| N | P46S | 1 | 2.86 |

## Проверка D614G

Обнаружена: 1 образцов; страны: USA.

## Короткое сравнение с дипломом

- В дипломе наиболее вариабельным указан ген **S**.
- В текущей выборке по листу `Nextclade_Gene_Summary` на первом месте: **ORF1a**.
- Это не противоречие само по себе: текущая выборка намного меньше и имеет другой состав.
- Сильная сторона текущего проекта: отчет обновляется автоматически после пересборки Excel.
