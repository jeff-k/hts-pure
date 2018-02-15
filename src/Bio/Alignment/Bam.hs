module Bio.Alignment.Bam (openBam,closeBam,alignments,pileup,Bamfile,
                          header,refs,mapq) where

import System.IO

import qualified Data.ByteString as BS
import Data.ByteString.Char8

import Data.ByteString (ByteString)
import qualified Data.ByteString.Lazy as L

import Codec.Compression.Zlib.Raw

import Conduit
import Data.Conduit.Cereal

import Data.Serialize
import Data.Serialize.Get
--
--import Data.Binary
--import Data.Binary.Get
--import Data.Word
--import Data.Int
--import Data.Bits

import Data.List (findIndex)

import Bio.Alignment.BamIndex

import Bio.Data.Cigar
import Bio.Data.Location

import Control.Applicative
import Control.Monad

import GHC.IO.Handle (hDuplicate)

data Bgzf = Bgzf {
    id1     :: Int,
    cdata   :: ByteString,
    isize   :: Int,
    bsize   :: Int
}

instance Serialize Bgzf where
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
    cdata <- getByteString (fromIntegral (bsize - xlen - 19))
    crc32 <- getWord32le
    isize <- fromIntegral <$> getWord32le
    
    return $ Bgzf id1 cdata isize bsize

instance Show Bgzf where
  show bz = (show . id1 $ bz) ++ "\t" ++ (show . isize $ bz) ++
            "\t" ++ (show . bsize $ bz) ++ "\n"

data Header = Header { text :: Maybe Text, refs :: [Text] }

instance Show Header where
    show s =
      case text s of
        Nothing -> printRefs $ refs s
        Just t -> "@HD\t" ++ (unpack t) ++ "\n" ++ (printRefs $ refs s)
      where
        printRefs rs =
          concatMap (\ x -> "@SQ\tSN:" ++ (unpack x) ++ "\n") rs

data Alignment = Alignment {
    pos         :: Int,
    refID       :: Int,
    mapq        :: Int,
    read_name   :: Text,
    seq_a       :: B.ByteString,
    cigar       :: [Cigar],
    flag        :: Int
}

instance Show Alignment where
    show a = read_name a ++ "\t" ++ (show . flag $ a) ++ "\t" ++
             (show . pos $ a) ++ "\t" ++ (show . refID $ a) ++ "\t" ++
             (show . mapq $ a) ++ "\t" ++ seq_a a ++ "\t" ++
             (show . cigar $ a) ++ "\n"

instance Serialize Alignment where
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
    read_name <- E.decodeUtf8 <$> getByteString l_read_name
    cigar_ops <- replicateM n_cigar_op getWord32le
    seq <- getByteString (div (l_seq + 1) 2)
    qual <- getByteString l_seq 
    tags <- getByteString (l - 32 - l_read_name - (n_cigar_op * 4) - l_seq -
                           div (l_seq + 1) 2)
    return $ Alignment pos
                       refID
                       mapq
                       read_name)
                       seq -- (concatMap readb (B.unpack seq))
                       []  -- (map readCig cigar_ops)
                       0

data Bam = Bam { header      :: Header,
                 index       :: Maybe Index,
                 pileup      :: Pos -> ConduitM () Alignment IO (),
                 readBam     :: ConduitM () Alignment IO ()
               }

--readb :: Word8 -> String
readb s =
    let l = (fromIntegral $ 15 .&. s)
        b = fromIntegral $ s `shiftR` 4
        t = "=ACMGRSVTWYHKDBN" in
            [t!!b, t!!l]

getRef :: Get Text
getRef = do
    l_name <- fromIntegral <$> getWord32le
    name <- E.decodeUtf8 <$> getByteString l_name 
    l_ref <- fromIntegral <$> getWord32le
    return name

instance Serialize Header where
  put = undefined
  get = do
    getByteString 4
    l_header <- fromIntegral <$> getWord32le
    t <- E.decodeUtf8 <$> getByteString l_header
    n_ref <- fromIntegral <$> getWord32le
    refs <- replicateM n_ref getRef
    case l_header of
      0 -> return $ Header Nothing refs
      _ -> return $ Header (Just t) refs

deZ :: Bgzf -> ByteString
deZ block = L.toStrict $ decompressWith params (L.fromStrict (cdata block))
  where params = defaultDecompressParams
  --DecompressParams (isize block) 2**16 Nothing True

-- convert a stream of bgzf blocks into a stream of bytes 
bgzfConduit :: MonadThrow m => ConduitM ByteString Bgzf m ()
bgzfConduit = conduitGet2 get

bamConduit :: MonadThrow m => ConduitM ByteString Alignment m ()
bamConduit = conduitGet2 get

openBam :: FilePath -> IO Bam
openBam path = do
  h <- openFile path ReadMode
  hSetBinaryMode h True

  -- get first block to build header
  header <- --

  -- attempt to load index
  index <- --

  let
    pileup p = case index of
      Just i -> 
        hSeek h AbsoluteSeek $ fst $ head (uncurry (offsets i (ref p))
                                           (interval p))
      Nothing -> sourceHandleBS h .| bgzfConduit .| mapC deZ .| bamConduit

    readBam = sourceHandleBS h .| bgzfConduit .| mapC deZ .| bamConduit

  return $ Bamfile header index pileup readBam
