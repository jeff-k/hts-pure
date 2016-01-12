import System.IO
import System.Directory

import System.Environment (getArgs)

import Control.Applicative
import Control.Monad

import Index
import Bam

main :: IO ()
main = do

    path <- getArgs
    index <- runGet getIndex <$> L.readFile (path!!0 ++ ".bai")

    h <- openFile (path!!0) ReadMode  
    bs <- runGet blocks <$> L.hGetContents h

    let y = runGet getBamfile $ L.concat  ((map (\x -> (decompressWith defaultDecompressParams (L.fromStrict . cdata $ x)))) bs) in
        mapM_ putStrLn (map show (alignments y))
