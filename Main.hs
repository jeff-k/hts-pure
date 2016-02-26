import System.IO
import System.Directory

import qualified Data.ByteString.Lazy as L

import Data.Binary.Get

import System.Environment (getArgs)

import Control.Applicative

import qualified Data.Map as M
import HTS
import Index
import Bam

main = do 
    args <- getArgs
    case args of
        [path] -> parseHeader path
        [path, ref] -> dumpRef path ref

parseHeader path = do
    h <- openFile path ReadMode  
    bamf <- bamfile h
--    print $ M.lookup coords (intervals index)
    print $ header bamf

dumpRef path ref = do
    index <- runGet getIndex <$> L.readFile (path ++ ".bai")
    putStrLn "asdf"
