module Bio.Data.Cigar (Cigar, readCig, cigLen) where

import Data.Bits
import Data.Word

data Cigar = Cigar { opLen :: Int, op :: Char }

instance Show Cigar where
    show s = (show . opLen $ s) ++ [op s]

readCig :: Word32 -> Cigar 
readCig s = Cigar opLen $ t !! op
    where opLen = fromIntegral $ s `shiftR` 4
          op = fromIntegral $ 7 .&. s
          t = "MIDNSHP=X"

-- Sequence length should be sum of M/I/S/=/X ops
cigLen :: [Cigar] -> Int
cigLen = undefined

-- Return the offset within the query sequence of the target position with
-- respect to a reference position. (compensate for indels and hard clipping).
--readOffset :: Int -> Int -> [Cigar] -> a -> Int
--readOffset _ _ [] a = a
--readOffset start target (c:cs) a -> readOffset start target (x + a)
--  where
--    tof = target - start
--    x = case (op c) of
--      'M' -> opLen c
--      'D' -> opLen c
--      'I' -> (-1) * $ opLen c

