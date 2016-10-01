# hts-pure
Module for manipulating next generation sequencing data, currently supporting
bam and bam index files.

## Examples
Filter reads with a mapping quality less than or equal to 25:

```haskell
import Bio.Alignment.Bam

main = do
  bam <- openBam "bamfile.bam" Nothing
  filter (\r -> mapq r <= 25) <$> alignments bam >>= print
```

Using an index to seek within a bam file:

```haskell
import Bio.Alignment.Bam
import Bio.Alignment.BamIndex

main = do
  index <- openIndex "bamfile.bam.bai"
  bam <- openBam "bamfile.bam.bai" (Just index)
  pileup bam (Pos 5 20000 30000) >>= print
```

## Issues
  * Bamfile handle is being closed after `pileup`, `alignments`, or `header`
  * To include all reads that overlap the requested interval, subtract the
    max read length from the begining of the position interval
