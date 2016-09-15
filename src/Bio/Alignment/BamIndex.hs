module Bio.Alignment.BamIndex (getIndex,getOffset,Index) where

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as Bchar

import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits

import Control.Applicative
import Control.Monad

import qualified Data.Map as M
import Bio.Data.Location

data Chunk = Chunk {beg::Word64, end::Word64}
data Bin = Bin {m::Int, b_chunks::[Chunk]}

data Index = Index {contigs::[ContigIndex], n_no_coor::Int}

data ContigIndex = ContigIndex {bins::[Bin], ioffsets::[Int]}

instance Show Chunk where
    show c = (show (beg c)) ++ "-" ++ (show (end c))

instance Show Index where
    show i = "Contigs:\t" ++ show (length (contigs i)) ++ "\n\t"
                          ++ show (contigs i) ++ "\n\t"

instance Show Bin where
    show b = show (m b) ++ (concat (map (\x -> "\t" ++ (show x) ++ "\n")
                                   (b_chunks b)))

instance Show ContigIndex where
    show c = show (bins c)

reg2bin :: Word64 -> Word64 -> Word64
reg2bin b e
    | log_shift 14 b e = calc_bin 14 15 b
    | log_shift 17 b e = calc_bin 17 12 b
    | log_shift 20 b e = calc_bin 20 9 b
    | log_shift 23 b e = calc_bin 23 6 b
    | log_shift 26 b e = calc_bin 26 3 b
    | otherwise = 0
    where calc_bin l c b = (((1 `shiftL` c) - 1) `div` 7) + (b `shiftR` l)
          log_shift l b e = (b `shiftR` l) == (e `shiftR` l)

bin2reg :: Word64 -> Word64 -> Word64
bin2reg b e =
    undefined
 
contigIndex :: Get ContigIndex
contigIndex = do
    n_bin <- fromIntegral <$> getWord32le
    bins <- replicateM n_bin bin
    n_intv <- fromIntegral <$> getWord32le
    offsets <- replicateM n_intv getWord64le
    return $ ContigIndex bins (fromIntegral <$> offsets)

chunk :: Get Chunk
chunk = do
    beg <- getWord64le
    end <- getWord64le
    return $ Chunk beg end

bin :: Get Bin
bin = do
    bin_id <- fromIntegral <$> getWord32le
    n_chunks <- fromIntegral <$> getWord32le
    b_chunks <- replicateM n_chunks chunk 
    return $ Bin bin_id b_chunks

getIndex :: Get Index
getIndex = do
    getByteString 4
    n_ref <- fromIntegral <$> getWord32le
    cs <- replicateM n_ref contigIndex
    empty <- isEmpty
    if empty
        then return $ Index cs 0
        else do n_no_coor <- fmap fromIntegral getWord32le
                return $ Index cs n_no_coor

voffs :: [Chunk] -> [(Word64, Word64)]
voffs i = map f i where
    f v = ((beg v) `shiftR` 16, (beg v) .&. 65535)

getOffset :: Index -> Pos -> [(Word64, Word64)]
getOffset i c =
    case (interval c) of
        Just (beg, end) ->  voffs (b_chunks (b!!bin_index)) where
            bin_index = (fromIntegral (reg2bin ((fromIntegral beg)::Word64) ((fromIntegral end)::Word64)))::Int
            b = bins ((contigs i)!!ref)
            offsets = ioffsets ((contigs i)!!ref)
            ref = 0 -- get ref
        Nothing -> [(0, 0)]
