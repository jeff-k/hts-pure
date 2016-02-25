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
    [path, coords] <- getArgs
    index <- runGet getIndex <$> L.readFile (path ++ ".bai")
--    h <- openFile path ReadMode  
--    bamf <- bamfile h
    print $ M.lookup coords (intervals index)
--    print $ header bamf
