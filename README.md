# hts-pure
Module for manipulating next generation sequencing data, currently supporting
Bam and Bam index files.

Alignments are instances of Binary and Show to support reading from Bam and
writing to Sam format. The namespace anticipates Cram and additional binary
filetypes will be supported.

This implementation intends to match the SAMv1 specification available at
https://github.com/samtools/hts-specs

## Examples
Filter reads with a mapping quality less than 25:

```haskell
import Bio.Alignment.Bam

main = do
  bam <- openBam "bamfile.bam" Nothing
  filter (\r -> mapq r < 25) <$> alignments bam >>= print
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
  * A conduit or pipes implementation would be better than lazy IO
  * To include all reads that overlap an interval, subtract the max read length
    from the begining of the position interval. TODO: parse CIGAR strings 
  * Native 0-based positions are assumed. This is an ontological can of worms.
  * `put` function of Binary instance for alignments is not implemented.
  * Indexing does not "jump gaps"; it only seeks to the first alignment in the
    requested interval and reads until the last covering read is hit
