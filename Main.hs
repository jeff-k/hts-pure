import System.IO
import System.Directory

import qualified Data.ByteString.Lazy as L

import Data.Binary.Get

import System.Environment (getArgs)

import Control.Applicative

import Index
import Bam
import HTS

main = do 
    args <- getArgs
    case args of
        [path] -> parseHeader path
        [path, coord] -> dumpRef path $ parseCoord coord

parseHeader path = do
    h <- openFile path ReadMode  
    bamf <- bamfile h
    print $ header bamf

dumpRef path coord = do
    index <- runGet getIndex <$> L.readFile (path ++ ".bai")
    print $ getOffset index coord
    h <- openFile path ReadMode
    stream <- bamSeek h index coord
    print stream
