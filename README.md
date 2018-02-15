# hts-pure
Module for manipulating high-throughput sequencing data with streaming data
structures, targeting the samtools family of file types.

This implementation intends to match the SAMv1 specification available at
https://github.com/samtools/hts-specs

## Examples
Filter reads with a mapping quality less than 25:

```haskell
import Conduit
import Bio.Alignment.Bam

main = do
  bam <- openBam "bamfile.bam"
  runConduitRes
        $ readBam bam .| filterC (\r -> mapq r < 25) .| printAlignment
```

Using an index to seek within a bam file:

```haskell
import Conduit
import Bio.Alignment.Bam
import Bio.Alignment.BamIndex

main = do
  bam <- openBam "bamfile.bam"
  runConduitRes $ pileup bam ('Chr4', 40000, 50000) .| printAlignment
```

## Notes
  * **This is an experimental project that doesn't work**
  * To include all reads that overlap an interval, subtract the max read length
    from the begining of the position interval. TODO: parse CIGAR strings 
  * Native 0-based positions are assumed. This is an ontological can of worms.
  * `put` function of Binary instance for alignments is not implemented.
  * Indexing does not "jump gaps"; it only seeks to the first alignment in the
    requested interval and reads until the last covering read is hit
