# Входные файлы для филогенетической сети

Эти файлы подготовлены для внешних инструментов построения филогенетических сетей.

- `alignment_short.fasta`: выравненные последовательности с короткими именами образцов.
- `network_metadata.tsv`: полная таблица метаданных для окраски и аннотации.
- `network_traits.tsv`: компактная таблица признаков для быстрого импорта.
- `alignment_short_mapping.tsv`: соответствие коротких имен исходным FASTA-описаниям.

Рекомендуемый порядок работы:

1. Используйте `alignment_short.fasta` как файл последовательностей.
2. Используйте `network_traits.tsv` или `network_metadata.tsv` для окраски по стране, месяцу, clade, lineage или D614G.
3. Постройте сеть в PopART, R или другом внешнем инструменте.

Воспроизводимость:

- Лучше пересоздавать эти файлы через Docker workflow, описанный в основном README проекта.
- Docker-образ фиксирует версию Nextclade CLI и tag SARS-CoV-2 dataset.

Сгенерированные пути:

- alignment: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\alignment_short.fasta`
- metadata: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\network_metadata.tsv`
- traits: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\network_traits.tsv`
- mapping: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\alignment_short_mapping.tsv`
