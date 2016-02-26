module Bam (bamfile,header) where

import System.IO

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as Bchar

import Codec.Compression.Zlib.Raw

import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits

import qualified Data.Map as M

import Control.Applicative
import Control.Monad

data Bgzf = Bgzf {id1::Int, cdata::B.ByteString, isize::Int} deriving Show

data Header = Header {text::String, refs::[Contig]}

data Cigar = Cigar {op_len::Int, op::Char}

data Alignment = Alignment {pos::Int, refID::Int, r::Int, read_name::String, seq_a::String, cigar::[Cigar]}

data Contig = Contig {name::String} deriving Show

instance Show Alignment where
    show a = (read_name a) ++ "\t" ++ (show . pos $ a) ++ "\t" ++ (show . refID $ a) ++ "\t" ++ (show . r $ a) ++ "\t" ++ (seq_a a) ++ "\t" ++ (show . cigar $ a) ++ "\n"

instance Show Cigar where
    show s = (show . op_len $ s) ++ [(op s)]

instance Show Header where
    show s = text s ++ show (length (refs s)) ++ "\n\t" ++ show (refs s) 

data Bamfile = Bamfile {header::Header, alignments::[Alignment]}

instance Show Bamfile where
    show b =  show (header b) ++ "\nt" ++ show (alignments b)

blocks :: Get [Bgzf]
blocks = do
    empty <- isEmpty
    if empty
        then return []
        else do
            block <- getBgzf 
            return (block:[])
--            rest <- blocks
--            return (block:rest)

readb :: Word8 -> String
readb s =
    let l = (fromIntegral $ 15 .&. s)
        b = fromIntegral $ s `shiftR` 4
        t = "=ACMGRSVTWYHKDBN" in
            [t!!b, t!!l]

readcig :: Word32 -> Cigar 
readcig s =
    let op_len = fromIntegral $ s `shiftR` 4
        op = fromIntegral $ 7 .&. s
        t = "MIDNSHP=X" in
    Cigar op_len (t!!op)

getAlignments :: Get [Alignment]
getAlignments = do
    empty <- isEmpty
    if empty
        then return []
        else do
            l <- fromIntegral <$> getWord32le
            refID <- fromIntegral <$> getWord32le
            pos <- fromIntegral <$> getWord32le
            l_read_name <- fromIntegral <$> getWord8
            mapq <- fromIntegral <$> getWord8 
            bin <- getWord16le
            n_cigar_op <- fromIntegral <$> getWord16le
            flag <- getWord16le
            l_seq <- fromIntegral <$> getWord32le
            next_refID <- getWord32le
            next_pos <- getWord32le
            tlen <- fromIntegral <$> getWord32le
            read_name <- getByteString l_read_name
            cigar_ops <- replicateM n_cigar_op getWord32le
            seq <- getByteString (div (l_seq + 1) 2)
            qual <- getByteString l_seq 
            tags <- getByteString (l - 32 - l_read_name - (n_cigar_op * 4) - l_seq - (div (l_seq + 1) 2))
--            rest <- getAlignments
            return ((Alignment refID pos mapq (Bchar.unpack read_name) (concat (map readb (B.unpack seq))) (map readcig cigar_ops)):[])

ref :: Get Contig
ref = do
    l_name <- fromIntegral <$> getWord32le
    name <- getByteString l_name 
    l_ref <- fromIntegral <$> getWord32le
    return $ Contig (Bchar.unpack name)

getHeader :: Get Header
getHeader = do
    getByteString 4
    l_header <- fromIntegral <$> getWord32le
    t <- getByteString l_header
    n_ref <- fromIntegral <$> getWord32le
    refs <- replicateM n_ref ref
    return $ Header (Bchar.unpack t) refs

getBamfile :: Get Bamfile
getBamfile = do
    h <- getHeader
    as <- getAlignments
    return $ Bamfile h as

getBgzf :: Get Bgzf
getBgzf = do
    id1 <- fromIntegral <$> getWord8
    id2 <- getWord8
    cm <- getWord8
    flg <- getWord8
    mtime <- fromIntegral <$> getWord32le
    xfl <- fromIntegral <$> getWord8
    os <- fromIntegral <$> getWord8
    xlen <- fromIntegral <$> getWord16le
    si1 <- getWord8
    si2 <- getWord8
    slen <- getWord16le
    bsize <- fromIntegral <$> getWord16le
    flags <- getByteString (xlen - 6)
    cdata <- getByteString (bsize - xlen - 19)
    crc32 <- getWord32le
    isize <- fromIntegral <$> getWord32le
    return $ Bgzf id1 cdata isize

dParam :: Bgzf -> DecompressParams
dParam block =
    DecompressParams (decompressWindowBits d) (isize block) Nothing
    where
        d = defaultDecompressParams

bamfile :: Handle -> IO Bamfile
bamfile h = do
    bs <- runGet blocks <$> L.hGetContents h
    return $ runGet getBamfile $ L.concat ((map (\x -> (decompressWith defaultDecompressParams (L.fromStrict . cdata $ x)))) bs)

