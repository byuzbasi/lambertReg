# Data, attribution and license

## QSAR aquatic toxicity

Davide Ballabio, Matteo Cassotti, Viviana Consonni, and Roberto Todeschini,
*QSAR aquatic toxicity* (2014), UCI Machine Learning Repository.
[Dataset DOI](https://doi.org/10.24432/C5SG7H) ·
[Provider](https://archive.ics.uci.edu/dataset/505/qsar%2Baquatic%2Btoxicity) ·
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

`qsar/data/qsar/` holds the unmodified archive and its original CSV member.
The CSV has no header and uses a semicolon delimiter. Columns 1–8 are
TPSA(Tot), SAacc, H-050, MLOGP, RDCHI, GATS1p, nN and C-040. Column 9 is
the quantitative LC50 response on the source's negative log scale.
There are 546 rows. The input manifest records the SHA-256 values.

## MAGIC maize

Matteo Dell’Acqua (2015), *MAGIC Maize RIL genotyping*, Figshare, version 1.
[Dataset DOI](https://doi.org/10.6084/m9.figshare.1437449.v1) ·
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
`magic-maize/data/maize/MMlines.geno.gz` is the unmodified genotype file.

Matteo Dell’Acqua et al. (2015), *Genetic properties of the MAGIC maize
population: a new platform for high definition QTL mapping in Zea mays*,
*Genome Biology* 16, 167. [Article and supplementary data](https://doi.org/10.1186/s13059-015-0716-z).
`Table_S3.xlsx` is the marker map and `Table_S6.xlsx` is the source phenotype
table. The publisher applies [CC0](https://creativecommons.org/publicdomain/zero/1.0/)
to data made available with the article unless otherwise stated. The article
is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

`marker_inventory.csv` is a derived annotation in genotype-column order.
It records marker positions and allele coding derived from the genotype
file and marker map. The runner aligns the `PH` response (plant height in
centimetres) by the source line identifier. The genotype file has 529 lines
and 54,234 marker columns before training-specific filtering. Input manifests
record the bytes and SHA-256 values of these files.

The third-party data retain their source terms. The package and study code
are under the repository's GPL-3 license.
