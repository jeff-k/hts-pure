module HTS (Coord(Coord), parseCoord, interval) where

import Data.Map as M

data Coord = Coord {ref::String, interval::Maybe (Int, Int)}

instance Show Coord where
    show c = case (interval c) of
        Just (beg,end) -> rangeStr (beg,end)
        Nothing -> ref c
        where
            rangeStr (b, e)
                | b == e = ref c ++ ":" ++ (show b)
                | otherwise = ref c ++ ":" ++ (show b) ++ "-" ++ (show e) 

parseCoord :: String -> Coord
parseCoord s = Coord s (Just (10,10000000))

