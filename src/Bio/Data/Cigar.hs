module Bio.Data.Cigar (Cigar, readcig) where

import Data.Bits
import Data.Word

data Cigar = Cigar {op_len::Int, op::Char}

instance Show Cigar where
    show s = (show . op_len $ s) ++ [(op s)]

readcig :: Word32 -> Cigar 
readcig s = Cigar op_len (t!!op)
    where op_len = fromIntegral $ s `shiftR` 4
          op = fromIntegral $ 7 .&. s
          t = "MIDNSHP=X"
