Explanation of the data tables:

  * table_of_crossovers_in_trisomies.csv - Contains the crossovers that were identified in trisomies via LD-PGTA.
  * table_of_crossovers_in_maternal_genomes_of_disomies.csv - Contains the crossovers that were identified in the maternal genome of disomies via LD-CHASE.
  * table_of_crossovers_in_paternal_genomes_of_disomies.csv - Contains the crossovers that were identified in the paternal genome of disomies via LD-CHASE.

The columns in all the tables are:
  (1) chr_id - Chromosome number.
  (2) embryo_id - A unique identifier that was given to each embryo.
  (3) crossover_position - The position of the crossover, where the chromosome length is normalized to 1.
  (4) z_score - a measurement of the confidence in calling the crossover, which is explained in the subsection "Identifying meiotic crossovers" of the article.

Technical details about calling the crossovers:

* Only cases where the genomic windows covered at least 50% of the chromosome were included in the analysis.
* The minimal z_score for calling a crossover was set to 1.96.
* A crossover was called by examining at least 15 genomic windows on both sides of the crossover position.
* The minimal distance between called crossovers is 2% of the chromosome length.
