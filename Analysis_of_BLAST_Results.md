Analysis of BLAST Results
================
Allison Bogisich <asbogisich@dons.usfca.edu>
October 16, 2017

Introduction
============

Only within the past decade or so has the technical capability to sequence complex biological samples been able to take off to such an astounding extent. Recent developments in phylogenetic community analysis (Lozupune and Knight,2005) with high-throughput pyrosequencing methods (Hamandy et al.,2008) have established the criteria for not only sequencing microbiome samples, but the very foundations of bacterial forensics. Fierer and his team have conducted a series of studies that demonstrate the strengths and weaknesses of the field (Fierer et al., 2010). They found that skin-associated bacterial communities are surprisingly diverse, and that they are independently stable enough to assist traditional forensic evaluations.

Conventional methods of obtaining forensic results from human DNA require suffient amounts of blood, tissue, semen, or saliva on an object. However, it is often difficult to obtain a large (and uncontaminated) enough sample to sequence. In order to boost standard results, recovering bacterial DNA from touched surfaces may be far easier, espeically for identifying objects where clear fingerprints can't be obtained (e.g. fabrics, smudged surfaces, highly textured surfaces) (Ibid.) Given how abundant bacterial cells are on the surface of skin and on shed epidermal cells (Fredricks, 2001) and how highly personalized bacterial communities are implicates that more research is needed. In this post-analysis of Fierer's research, I compare the accuracy and reliability of the resulting sequence matches from male versus female subjects.

Methods
=======

This post-experimental review is based on .csv output files from a BLAST search conducted from trimmed and quality checked fasta files from the Fierer team study. The BLAST search matched his swabbed and sequenced samples against the GenBank database from NCBI Sequence Read Archive study number ERP022657. The above process and the resulting files can be found in a Github repository 2017-usfca-cs-640/ABogisich-week06-homework-QC-and-BLAST. Using bioinformatic applications in R and R Studio, the sample data were reorganized and joined into easily readable histograms to compare the reliability of microbial sequencing data from fanua present on female and male subjects.

Sample origin and sequencing
----------------------------

Fierer and his team swabbed for skin-associated bacteria from different epidermal regions for each of their three studies. For the keyboard study, three participants were swabbed on the ventral surface of the distal joint of each fingertip. All individuals were healthy 20-35 year olds that had not taken antibiotics at least six months prior to swabbing. In the "storage" study, autoclaved cotton-tipped swabs moistened with sterile solution were again used, this time to sample the right axilary (armpit) skin surface sixteen times, from two healthy adults. Lastly, in the computer mouse study nine healthy adults were recruited (four female and five male, all 20-35 years of age) who worked in the same University of Colorado building. Using the swabbing technique outlined above, the entirety of the exposed surface of their dominant hand's palm (used to control the computer mice) were awabbed. Palm surfaces were sampled at midday and the particpants had been following typical hand hygiene practices prior to sampling. Swabs were stored at -80 degrees before DNA extraction. The microbial communities on these participants were compared to a previously compiled database from 270 other hands sampled by Fierer and collaborators (Fierer et al., 2008; Costello et al., 2009). The 270 hands in the database were from left and right palm surfaces belonging to an equal proportion of both healthy male and female volunteers, aged 18-40 years old.

Post DNA extraction from the samples, sequences were processed and analyzed using pyrosequencing procedures akin to Fierer's work in 2008. Sequences were removed if shorter than 200 bp or larger than 300 in length, had a quality score lower than 25, had ambiguous characters or uncorrectable barcodes, or if it did not contain the primer sequence. Remaining sequences were assigned to samples via examination of the 12-nt barcode and then clustered into operational taxonomic units (OTUs). Representative sequences were selected for each OTU based on the longest length sequnce with the greast number of hits to other sequences within the OTU.

Computational
-------------

The program R version 3.4.2 in conjunction with R Studio version 1.1.383 interface for Windows desktop was used. In addition, several packages listed below from the R library which were downloaded in order to utilize more specific exploratory data analysis functions. The blast data set and its associated metadata were joined together for use within a single dataframe. First, output files from BLAST sequence search results in .csv format were then read into proper dataframes using subsequent data qualifiers to parse query sequence identifiers into appropriately labeled columns/reading frames. The empty matrix for the joined data was then created, with a function that would loop and append together all the results for each of the BLAST results. Then the tab delimited metadata was read in and joined with the BLAST output using dpyler, forming one large data table for manipulation. Finally, in order to produce visible data based histograms, the diplyr piping system was used to select a subset of rows matching certain criterion and then pull out columns.

Results
=======

``` r
# Be sure to install these packages before running this script
# They can be installed either with the intall.packages() function
# or with the 'Packages' pane in RStudio

# load packages
library("dplyr")
library("tidyr")
library("knitr")
```

``` r
# Output format from BLAST is as detailed on:
# https://www.ncbi.nlm.nih.gov/books/NBK279675/
# In this case, we used: '10 sscinames std'
# 10 means csv format
# sscinames means unique Subject Scientific Name(s), separated by a ';'
# std means the standard set of result columns, which are:
# 'qseqid sseqid pident length mismatch
# gapopen qstart qend sstart send evalue bitscore',


# this function takes as input a quoted path to a BLAST result file
# and produces as output a dataframe with proper column headers
# and the 'qseqid' column split into sample and seq number
read_blast_output <- function(filename) {
  data_in <- read.csv(filename,
                      header = FALSE, # files don't have column names in them
                      col.names = c("sscinames", # unique Subject Sci Name(s)
                                    "qseqid",    # Query Seq-id
                                    "sseqid",    # Subject Seq-id
                                    "pident",    # Percntge of identical matches
                                    "length",    # Alignment length
                                    "mismatch",  # Number of mismatches
                                    "gapopen",   # Number of gap openings
                                    "qstart",    # Start of alignment in query
                                    "qend",      # End of alignment in query
                                    "sstart",    # Start of alignment in subj
                                    "send",      # End of alignment in subject
                                    "evalue",    # Expect value
                                    "bitscore"))  # Bit score

  # Next we want to split the query sequence ID into
  # Sample and Number components so we can group by sample
  # They originally look like "ERR1942280.1"
  # and we want to split that into two columns: "ERR1942280" and "1"
  # we can use the separate() function from the tidyr library to do this
  # Note that we have to double escape the period for this to work
  # the syntax is
  # separate(column_to_separate,
  # c("New_column_name_1", "New_column_name_2"),
  # "seperator")
  data_in <- data_in %>%
    separate(qseqid, c("sample_name", "sample_number"), "\\.")
}
```

``` r
# this makes a vector of all the BLAST output file names, including
# the name(s) of the directories they are in
files_to_read_in <- list.files(path = "output/blast",
                               full.names = TRUE)

# We need to create an empty matrix with the right number of columns
# so that we can rbind() each dataset on to it
joined_blast_data <- matrix(nrow = 0,
                            ncol = 14)

# now we loop over each of the files in the list and append them
# to the bottom of the 'joined_blast_data' object
# we do this with the rbind() function and the function we
# made earlier to read in the files, read_blast_output()
for (filename in files_to_read_in) {
  joined_blast_data <- rbind(joined_blast_data,
                             read_blast_output(filename))
}
```

``` r
# Next we want to read in the metadata file so we can add that in too
# This is not a csv file, so we have to use a slightly different syntax
# here the `sep = "\t"` tells the function that the data are tab-delimited
# and the `stringsAsFactors = FALSE` tells it not to assume that things are
# categorical variables
metadata_in <- read.table(paste0("data/metadata/",
                                 "fierer_forensic_hand_mouse_SraRunTable.txt"),
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)

# Finally we use the left_join() function from dplyr to merge or 'join' the
# combined data and metadata into one big table, so it's easier to work with
# in R the `by = c("Run_s" = "sample_name")` syntax tells R which columns
# to match up when joining the datasets together
joined_blast_data_metadata <- metadata_in %>%
  left_join(joined_blast_data,
            by = c("Run_s" = "sample_name"))
```

``` r
# Here we're using the dplyr piping syntax to select a subset of rows matching a
# criteria we specify (using the filter) function, and then pull out a column
# from the data to make a histogram. We don't need to tell the hist() function
# which data to use, because that's piped in, but we do have to give the
# hist() function the title and axis label we'd like to use for the figure
joined_blast_data_metadata %>%
  filter(env_material_s == "sebum") %>%
  pull(pident) %>%
  hist(main = "Percent Identity",
       xlab = "Percent")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github-ascii_identifiers/histogram%201-1.png)

Fig. 1 Percent Identity of Microbial DNA from Human Sebum

Percent identity refers to a quantitative measurement of the similarity between two sequences (DNA, amino acid or otherwise). In figure 1, we see the range of similarity between all microbial DNA sequences detected in human sebum samples. The vast majority of the microbial DNA collected matched completely with another sequence collected from skin sebum secretions, occurring at a very high frequency of at least 2,500. The next largest percent similarity occurred at about 87% similarity with other sequences, and no sequences had percent identities that fell below 80%.

``` r
joined_blast_data_metadata %>%
  filter(sex_s == "female") %>%
  pull(mismatch) %>%
  hist(main = "Female Microbial Survey",
       xlab = "Number of Mismatches")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github-ascii_identifiers/histogram%202-1.png) Fig. 2 Frequency of Mismatched Microbial Sequences from Female Samples

The number of mismatched sequences from female microbiome samples is skewed to the right/ a positive skew. The average sample had five or fewer mismatches in the sequence, with the tail trailing to the right. Very few samples (less than one-hundred) had more than twenty-five mismatches, while about six-hundred of the samples had between five and twenty-five mismatches.

``` r
joined_blast_data_metadata %>%
  filter(sex_s == "male") %>%
  pull(mismatch) %>%
  hist(main = "Male Microbial Survey",
       xlab = "Number of Mismatches")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github-ascii_identifiers/histogram%203-1.png) Fig. 3 Frequency of Mismatched Microbial Sequences from Male Samples

The distribution of mismatched sequences from male microbiome samples is generally skewed to the right, however it is not as strong a distribution as that of the female samples. Most frequently (about 1500 occurrences), male samples had between fifteen and twenty mismatches in the sequence, with the tail trailing to the right. Very few samples (about one-hundred) had more than twenty mismatches, while the majority of the samples (about two thousand) had between zero and twenty mismatches.

``` r
joined_blast_data_metadata %>%
  filter(sex_s == "female") %>%
  pull(bitscore) %>%
  hist(main = "Bit scores for Female Sequences",
       xlab = "Bit score Value")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github-ascii_identifiers/histogram%204-1.png) Fig. 4 Bit score Distribution for Sequenced Female Microbiome Samples

The bit score (S) formula takes into account gaps in the sequence; the higher the score, the better the alignment for the sequence/highly significant the match is. High bit scores can occur when the alignment length is longer, even when there are more mismatches than for other hits. Here in figure 4 we see that the distribution is strongly skewed left/ negatively skewed for sequences derived from female skin samples. With the tail trailing left, the majority of the sequenced microbial DNA had highly significant bit scores. The mode score falls between 400-425.

``` r
joined_blast_data_metadata %>%
  filter(sex_s == "male") %>%
  pull(bitscore) %>%
  hist(main = "Bit scores for Male Sequences",
       xlab = "Bit score Value")
```

![](Analysis_of_BLAST_Results_files/figure-markdown_github-ascii_identifiers/histogram%205-1.png) Fig. 5 Bit score Distribution for Sequenced Male Microbiome Samples

For microbial sequences derived from male skin samples, no strong pattern of distribution emerges.The most common bit score falls between 150-200, with a smaller, secondary clustering of scores falling between 325-400. The largest frequency for males is higher (~1800), and occurs at a lower bit score (150-200) versus the female bitscore mode of 400-425, at a lower frequency (~1000).

Discussion
==========

While homologies can be inferred from "excess" similarity, and "excess" similarity is based on statistical estimates, most bioinformatics researchers tend to describe similarity in terms of 'percent identity'. Although the common rule of thumb is that any two sequences are homologous if they are greater than 30% identical over their entire respective lengths, the 30% criterion misses many relatively detectable homologs (Pearson, 2014). In Figure 1, we can see that the percent identities of the sequenced microbial DNA from the human skin samples all were above eighty percent and most were at one-hundred percent. This implies that essentially all of the sample sequences either had an exactly identical homolog or were very near to it. This suggests that across all male and female skin samples, there were very few sequences that were entirely unique.

Several statistics allow comparison of hits across different searches returned from BLAST, and the ones examined here are mismatch and bit scores. In figures 2 and 3, the number of mismatches over the length of the alignment gives a rough idea of how closely sequences matched for female and male samples respectively. Based on the distributions, more mismatches at greater frequencies occurred in DNA matching of the male microbial surveys than in the female microbial surveys. This would suggest that the species making up the microbial communities of males are less documented/more novel than those species making up the female participants epidermal metagenomes. That not only are the microbiomes of invidiuals significantly different from each other, but that the microbiomes of males may be significantly different from those of female hominids. Further testing of this would need to be conducted in other environmental settings: different states, rural areas, etc. The bit score gender comparisons in figures 4 and 5, further portray biases in the accuracy and reliability of microbial sequence BLAST matches. While female bitscores in figure 4 maintained a highly negatively skewed distribution- with most sequences having very high alignment scores- males bitscores in figure 5 were markedly less reliable and robust. Since bit scores provide a constant statistical indicator for searching different databases of arious sizes or for searching the same database as it grows larger, bit scores are highly useful.

These apparent genedered differences have greater implications for the use of bacterial forensics, as biases in the ability to properly identify male and female microbiomes could potentially result in imbalanced or false convictions along the gender divide. Even if the overall reliablity of the statisitical analyses is strong for both genders, *if linking touched objects to one gender is even slightly easier and more reliable than for the opposite gender- that could have both bioethical and legal ramifications*. Further studies will still need to be conducted with greater and more equivalently gendered sample sizes to help improve the reliability and accuracy of microbial biome matching. Error could have possibly been compounded in what men versus women deem as "typical hygeine practices" prior to sampling, despite the fact that at least on palm surfaces, bacterial communites can recover within hours of washing (Fiere et al., 2008). The swabbing techniques used to pick up microbial communities in sebum mayb be slightly less effective on male versus female skin, which could further skew computational results. Especially needed are additional studies to assess how well the method works when pairing skin sebum to objects of different textures, or objects that come into contact with multiple skin locations on any given individual. As the reliablity of this method currently stands, submitting bacterial matches as forensic evidence could potentially become an institutialized form of judicial sexism.
