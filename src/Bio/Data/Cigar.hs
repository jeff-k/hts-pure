module Bio.Data.Cigar (Cigar, readcig) where

import Data.Bits
import Data.Word

data Cigar = Cigar { opLen :: Int, op :: Char }

instance Show Cigar where
    show s = (show . opLen $ s) ++ [op s]

readcig :: Word32 -> Cigar 
readcig s = Cigar opLen (t!!op)
    where opLen = fromIntegral $ s `shiftR` 4
          op = fromIntegral $ 7 .&. s
          t = "MIDNSHP=X"
