module Bio.Data.Cigar (Cigar, readcig) where

data Cigar = Cigar {op_len::Int, op::Char}

instance Show Cigar where
    show s = (show . op_len $ s) ++ [(op s)]

readcig :: Word32 -> Cigar 
readcig s =
    let op_len = fromIntegral $ s `shiftR` 4
        op = fromIntegral $ 7 .&. s
        t = "MIDNSHP=X" in
    Cigar op_len (t!!op)
