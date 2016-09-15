module Bio.Data.Location (Pos, parsePos, coverage, interval) where

import Data.Map as M

data Pos = Pos {ref::String, interval::Maybe (Int, Int)}

instance Show Pos where
    show c = case (interval c) of
        Just (beg,end) -> rangeStr (beg,end)
        Nothing -> ref c
        where
            rangeStr (b, e)
                | b == e = ref c ++ ":" ++ (show b)
                | otherwise = ref c ++ ":" ++ (show b) ++ "-" ++ (show e) 

coverage :: [Pos] -> [Pos]
coverage = undefined

parsePos :: String -> Pos
parsePos = undefined
--parsePos s = Pos s (Nothing)

