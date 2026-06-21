# PopART Inputs

These NEXUS files are prepared for PopART haplotype network analysis.

Recommended PopART workflow:

1. Open one of the `popart_*.nex` files in PopART.
2. Choose a network method, usually Median Joining or TCS.
3. Use the trait legend for coloring by country, month, clade, lineage, or D614G.
4. Export the final figure from PopART as SVG/PNG/PDF.

Reproducibility:

- Prefer regenerating these files through the Docker workflow described in the project README.
- The Docker image freezes the Nextclade CLI version and the SARS-CoV-2 dataset tag used upstream.

Generated files:

- `popart_country.nex`: coloring by `country`.
- `popart_month.nex`: coloring by `month`.
- `popart_clade.nex`: coloring by `clade`.
- `popart_lineage.nex`: coloring by `lineage`.
- `popart_d614g.nex`: coloring by `d614g`.

For each NEXUS file, a matching `*_trait_labels.tsv` file maps PopART-safe
labels such as `South_Korea` back to the original metadata values.

Note: with 35 genomes the network will be much smaller than the large diploma figures.
To make a similar dense network, the dataset needs many more genomes across months and countries.
