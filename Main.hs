import System.IO
import System.Directory

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString as B

import Codec.Compression.Zlib.Raw
import Data.Binary.Get
import Data.Word
import Data.Int
import Data.Bits

import System.Environment (getArgs)

import Control.Applicative
import Control.Monad

import Coord
import Index
import Bam

--dump :: IO (path, Maybe Coord)
--dump path c =
--    case c of
--        Just c -> 

main = do 
    [path, coords] <- getArgs
    index <- runGet getIndex <$> L.readFile (path ++ ".bai")
--    h <- openFile path ReadMode  
--    bamf <- bamfile h
    print $ intervals index
--    print $ header bamf
