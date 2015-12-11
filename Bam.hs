import System.IO

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as Bchar

import System.Environment (getArgs)
import Codec.Compression.Zlib.Raw

import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits

import Control.Applicative
import Control.Monad

data Bgzf = Bgzf {id1::Int, cdata::B.ByteString, isize::Int} deriving Show
data Chunk = Chunk {beg::Int, end::Int} deriving Show
data Bin = Bin {m::Int, chunks::[Chunk]}

data Index = Index {contigs::[ContigIndex], n_no_coor::Int}
data Header = Header {text::String, refs::[Contig]}

data Cigar = Cigar {op_len::Int, op::Char}

data Alignment = Alignment {pos::Int, refID::Int, r::Int, read_name::String, seq_a::String, cigar::[Cigar]}

data Contig = Contig {name::String} deriving Show
data ContigIndex = ContigIndex {bins::[Bin], ioffsets::[Int]}

instance Show Alignment where
    show a = (read_name a) ++ "\t" ++ (show . pos $ a) ++ "\t" ++ (show . refID $ a) ++ "\t" ++ (show . r $ a) ++ "\t" ++ (seq_a a) ++ "\t" ++ (show . cigar $ a) ++ "\n"

instance Show Cigar where
    show s = (show . op_len $ s) ++ [(op s)]

instance Show Header where
    show s = text s ++ show (length (refs s)) ++ "\n\t" ++ show (refs s) 

instance Show Index where
    show i = "Contigs:\t" ++ show (length (contigs i)) ++ "\n\t" ++ show (contigs i) ++ "\n\t"

data Bamfile = Bamfile {header::Header, alignments::[Alignment]}

instance Show Bamfile where
    show b =  show (header b) ++ "\nt" ++ show (alignments b)

instance Show Bin where
    show b = show (m b)

instance Show ContigIndex where
    show c = show (bins c)

blocks :: Get [Bgzf]
blocks = do
    empty <- isEmpty
    if empty
        then return []
        else do
            block <- getBgzf 
            rest <- blocks
            return (block:rest)

readb :: Word8 -> String
readb s =
    let l = (fromIntegral $ 15 .&. s) in
        let b = fromIntegral $ s `shiftR` 4 in
            let t = "=ACMGRSVTWYHKDBN" in
                [t!!b, t!!l]

readcig :: Word32 -> Cigar 
readcig s =
    let op_len = fromIntegral $ s `shiftR` 4 in
        let op = fromIntegral $ 7 .&. s in
            let t = "MIDNSHP=X" in
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
            rest <- getAlignments
            return ((Alignment refID pos mapq (Bchar.unpack read_name) (concat (map readb (B.unpack seq))) (map readcig cigar_ops)):rest)

ref :: Get Contig
ref = do
    l_name <- fromIntegral <$> getWord32le
    name <- getByteString l_name 
    l_ref <- fromIntegral <$> getWord32le
    return $ Contig (Bchar.unpack name)

contigIndex :: Get ContigIndex
contigIndex = do
    n_bin <- fromIntegral <$> getWord32le
    bins <- replicateM n_bin bin
    n_intv <- fromIntegral <$> getWord32le
    offsets <- replicateM n_intv (getWord32le >> getWord32le)
    return $ ContigIndex bins (fromIntegral <$> offsets)

getHeader :: Get Header
getHeader = do
    getByteString 4
    l_header <- fromIntegral <$> getWord32le
    t <- getByteString l_header
    n_ref <- fromIntegral <$> getWord32le
    refs <- replicateM n_ref ref
    return $ Header (Bchar.unpack t) refs

chunk :: Get Chunk
chunk = do
    beg <- fromIntegral <$> getWord64le
    end <- fromIntegral <$> getWord64le
    return $ Chunk beg end

bin :: Get Bin
bin = do
    bin_id <- fromIntegral <$> getWord32le
    n_chunks <- fromIntegral <$> getWord32le
    chunks <- replicateM n_chunks chunk 
    return $ Bin bin_id chunks

getBamfile :: Get (Bamfile, String)
getBamfile = do
    h <- getHeader
    as <- getAlignments
    return $ (Bamfile h as, "Done")

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

getBlocks :: L.ByteString -> [B.ByteString]
getBlocks i = go decoder i
    where
        decoder = runGetIncremental getBgzf
        go :: Decoder Bgzf -> L.ByteString -> [B.ByteString]
        go (Done r _c block) input =
            (L.toStrict $ decompressWith (dParam block) (L.fromStrict . cdata $ block)) : go decoder (L.fromStrict r)
        go (Partial k) input =
            go (k $ Just (L.toStrict input)) input
        go (Fail rem _c msg) _input =
            error msg

main :: IO ()
main = do
    path <- getArgs
    i <- runGet getIndex <$> L.readFile (path!!0 ++ ".bai")

    h <- openFile (path!!0) ReadMode
    bs <- runGet blocks <$> L.hGetContents h
   
    x <- getBlocks <$> (L.hGetContents h) 
    print $ length x 
    let (y, _) = runGet getBamfile $ L.concat ((map (\x -> (decompressWith defaultDecompressParams (L.fromStrict . cdata $ x)))) bs) in
        mapM_ putStrLn (map show (alignments y))
