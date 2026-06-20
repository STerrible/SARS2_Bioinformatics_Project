# Phylogenetic Network Inputs

These files are prepared for external phylogenetic network tools.

- `alignment_short.fasta`: aligned sequences with short sample names.
- `network_metadata.tsv`: full metadata table for coloring and annotation.
- `network_traits.tsv`: compact trait table for quick import into network tools.
- `alignment_short_mapping.tsv`: mapping from short names to original FASTA descriptions.

Recommended workflow:

1. Use `alignment_short.fasta` as the sequence input.
2. Use `network_traits.tsv` or `network_metadata.tsv` for coloring by country, month, clade, lineage, or D614G.
3. Build the network in PopART, R, or another external tool.

Generated paths:

- alignment: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\alignment_short.fasta`
- metadata: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\network_metadata.tsv`
- traits: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\network_traits.tsv`
- mapping: `C:\Users\STerrible\PycharmProjects\SARS2_Bioinformatics_Project\results\phylogenetics\alignment_short_mapping.tsv`
