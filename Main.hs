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
        [path, ref] -> dumpRef path $ Coord ref Nothing
        [path, ref, pos] -> dumpRef path $ Coord ref (Just (read pos, read pos))

parseHeader path = do
    h <- openFile path ReadMode  
    bamf <- bamfile h
    print $ header bamf

dumpRef path coord = do
    index <- runGet getIndex <$> L.readFile (path ++ ".bai")
    h <- openFile path ReadMode
    stream <- bamSeek h index coord
    print $ stream
