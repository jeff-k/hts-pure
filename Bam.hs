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
data Chunk = Chunk {beg::Word64, end::Word64} deriving Show
data Bin = Bin {m::Int, b_chunks::[Chunk]}

data Index = Index {contigs::[ContigIndex], n_no_coor::Int}
data Header = Header {text::String, refs::[Contig]}

data Rtree = Rleaf Int | Rtree [Rtree]

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
    show b = show (m b) ++ show (b_chunks b)

instance Show ContigIndex where
    show c = show (bins c)

reg2bin :: Word64 -> Word64 -> Word64
reg2bin b e
    | (b `shiftR` 14) == (e `shiftR` 14) = ((((1 `shiftL` 15) - 1) `div` 7) + (b `shiftR` 14))
    | True = 0

--reg2bins :: Int -> Int -> [Bin]
--reg2bins _ _ = [Bin 3 3]

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
--            rest <- getAlignments
            return ((Alignment refID pos mapq (Bchar.unpack read_name) (concat (map readb (B.unpack seq))) (map readcig cigar_ops)):[])

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
    offsets <- replicateM n_intv getWord64le
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
    beg <- getWord64le
    end <- getWord64le
    return $ Chunk beg end

bin :: Get Bin
bin = do
    bin_id <- fromIntegral <$> getWord32le
    n_chunks <- fromIntegral <$> getWord32le
    b_chunks <- replicateM n_chunks chunk 
    return $ Bin bin_id b_chunks

getBamfile :: Get Bamfile
getBamfile = do
    h <- getHeader
    as <- getAlignments
    return $ Bamfile h as

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
            go (k . chead $ input) (ctail input)
        go (Fail rem _c msg) _input =
            error msg

chead :: L.ByteString -> Maybe B.ByteString
chead i =
    case L.toChunks i of
        bs:_ -> Just bs
        _ -> Nothing

ctail :: L.ByteString -> L.ByteString
ctail i =
    case L.toChunks i of
        _:tail -> (L.fromStrict (B.concat tail))
        _ -> L.empty 

voff :: Index -> (Int, Word64, Word64)
voff i =
    (x, ((beg v) `shiftR` 16), ((beg v) .&. 65535))
    where
        v = ((b_chunks $ (bins ((contigs i)!!0))!!b)!!0)
        x = m $ (bins ((contigs i)!!0))!!b
        b = 1

main :: IO ()
main = do
    path <- getArgs
    (m, vo, bo) <- voff <$> runGet getIndex <$> L.readFile (path!!0 ++ ".bai")

    h <- openFile (path!!0) ReadMode

    hSeek h AbsoluteSeek (fromIntegral vo)
    bs <- runGet blocks <$> L.hGetContents h       
    let x = (map (\x -> (decompressWith defaultDecompressParams (L.fromStrict . cdata $ x))) bs)!!0 in
        print $ runGet getAlignments $ L.drop (fromIntegral bo) x
