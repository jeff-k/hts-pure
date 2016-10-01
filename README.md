# hts-pure
Module for manipulating next generation sequencing data, currently supporting
bam and bam index files.

## Examples
Filter reads with a mapping quality less than or equal to 25:

  import Bio.Alignment.Bam

  main = do
    bam <- openBam "bamfile.bam" Nothing
    rs <- filter (\r -> 25 > mapq r) <$> alignments bam
    print rs

Using an index to seek within a bam file:

  import Bio.Alignment.Bam
  import Bio.Alignment.BamIndex

  main = do
    index <- openIndex "bamfile.bam.bai"
    bam <- openBam "bamfile.bam.bai" (Just index)
    rs <- pileup bam $ Pos 5 20000 30000
    print rs

## Issues
  * Bamfile handle is being closed after `pileup`, `alignments`, or `header`
  * To include all reads that overlap the requested interval, subtract the
    max read length from the begining of the position interval
