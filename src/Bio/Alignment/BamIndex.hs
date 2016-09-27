module Bio.Alignment.BamIndex (openIndex,offsets,Index,contigs) where

import System.IO

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as Bchar

import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits
import Data.List

import Control.Applicative
import Control.Monad

import Bio.Data.Location

data Chunk = Chunk {beg::Word64, end::Word64}
data Bin = Bin {m::Int, b_chunks::[Chunk]}

data Index = Index { contigs :: [ContigIndex],
                     n_no_coor :: Int,
                     offsets :: Int -> Integer -> Integer -> [(Integer, Integer)] }

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
    let offsets ref beg end = voffs (b_chunks ((bins ((cs)!!ref))!!(bin_index beg end))) 

    if empty
      then return $ Index cs 0 offsets
      else do
        n_no_coor <- fmap fromIntegral getWord32le
        return $ Index cs n_no_coor offsets
    where
      bin_index beg end = (fromIntegral (reg2bin ((fromIntegral beg)::Word64)
                                                 ((fromIntegral end)::Word64)))::Int
--          offsets = ioffsets ((contigs i)!!r)


voffs :: [Chunk] -> [(Integer, Integer)]
voffs i = map f i where
    f v = ((fromIntegral ((beg v) `shiftR` 16)), (fromIntegral ((beg v) .&. 65535)))

openIndex :: String -> IO Index
openIndex path = do
  h <- openFile path ReadMode
  hSetBinaryMode h True
  s <- L.hGetContents h
  let 
    i = runGet getIndex
  return $ (runGet getIndex s)
