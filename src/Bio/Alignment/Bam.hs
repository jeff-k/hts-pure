module Bio.Alignment.Bam (openBam,pileup,alignments,header) where

import System.IO

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as Bchar

import Codec.Compression.Zlib.Raw

import Data.Binary
import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits

import System.IO (isEOF)

import Bio.Alignment.BamIndex

import Bio.Data.Cigar
import Bio.Data.Location

import Control.Applicative
import Control.Monad

data Bgzf = Bgzf {
    id1     :: Int,
    cdata   :: L.ByteString,
    isize   :: Int,
    bsize   :: Int
}

instance Binary Bgzf where
  put = undefined
  get = do
    id1 <- fromIntegral <$> getWord8
    id2 <- getWord8
    cm <- getWord8
    flg <- getWord8
    mtime <- fromIntegral <$> getWord32le
    xfl <- fromIntegral <$> getWord8
    os <- fromIntegral <$> getWord8
    xlen <- fromIntegral <$> getWord16le
    si1 <- fromIntegral <$> getWord8
    si2 <- getWord8
    slen <- fromIntegral <$> getWord16le
    bsize <- fromIntegral <$> getWord16le
    _ <- getByteString (xlen - 6)
    cdata <- getLazyByteString (fromIntegral (bsize - xlen - 19))
    crc32 <- getWord32le
    isize <- fromIntegral <$> getWord32le
    
    return $ Bgzf id1 cdata isize bsize


instance Show Bgzf where
  show bz = (show . id1 $ bz) ++ "\t" ++ (show . isize $ bz) ++ "\t" ++ (show . bsize $ bz) ++ "\n"

data Header = Header {text::Maybe String, refs::[Contig]}

data Alignment = Alignment {
    pos         :: Int,
    refID       :: Int,
    r           :: Int,
    read_name   :: String,
    seq_a       :: String,
    cigar       :: [Cigar]
}

data Contig = Contig {name::String}

instance Show Alignment where
    show a = (read_name a) ++ "\t" ++ (show . pos $ a) ++ "\t" ++
             (show . refID $ a) ++ "\t" ++ (show . r $ a) ++ "\t" ++
             (seq_a a) ++ "\t" ++ (show . cigar $ a) ++ "\n"

instance Binary Alignment where
  put = undefined
  get = do
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
    return $ Alignment refID
                       pos
                       mapq
                       (Bchar.unpack read_name)
                       (concat (map readb (B.unpack seq)))
                       (map readcig cigar_ops)

instance Show Contig where
    show c = name c

instance Show Header where
    show s = case (text s) of
      Nothing -> concat $ map (\x -> (show x) ++ "\n") (refs s) 
      Just t -> t ++ "\n" ++ (concat $ map (\x -> (show x) ++ "\n") (refs s))

data Bamfile = Bamfile {header::Header, handle::Handle}

instance Show Bamfile where
    show b =  show (header b)

readb :: Word8 -> String
readb s =
    let l = (fromIntegral $ 15 .&. s)
        b = fromIntegral $ s `shiftR` 4
        t = "=ACMGRSVTWYHKDBN" in
            [t!!b, t!!l]

getRef :: Get Contig
getRef = do
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
    refs <- replicateM n_ref getRef
    case l_header of
      0 -> return $ Header Nothing refs
      otherwise -> return $ Header (Just (Bchar.unpack t)) refs


deZ :: Bgzf -> L.ByteString
deZ block = decompressWith params (cdata block)
  where params = defaultDecompressParams
  --DecompressParams (isize block) 2**16 Nothing True

getBlocks :: Get [L.ByteString]
getBlocks = do
  empty <- isEmpty
  if empty
    then return []
    else do
      b <- deZ <$> (get :: Get Bgzf)
      bs <- getBlocks
      return (b:bs)

getReads :: Get [Alignment]
getReads = do
  empty <- isEmpty
  if empty
    then return []
    else do
      r <- get :: Get Alignment
      rs <- getReads
      return (r:rs)

openBam :: String -> IO Bamfile
openBam path = do
  h <- openFile path ReadMode
  hSetBinaryMode h True
  s <- L.hGetContents h
  return $ Bamfile (runGet getHeader s) h

alignments :: Bamfile -> IO [Alignment]
alignments b = do
  s <- L.hGetContents (handle b)
  return $ [runGet (get :: Get Alignment) $ L.concat $ runGet getBlocks s]

pileup :: Bamfile -> Index -> Pos -> IO [Alignment]
pileup b i pos =
  let (c, u) = (offsets i pos)!!0 in
    do
      hSeek (handle b) AbsoluteSeek c
      as <- alignments b
      return as
--      s <- L.hGetContents (handle b)
--      return $ [runGet (get :: Get Alignment) $ L.concat $ runGet getBlocks s]
