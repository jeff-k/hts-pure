module HTS (Coord, parseCoord, ITree) where

import Data.Map as M

data Coord = Coord {ref::String, interval::Maybe (Integer, Integer)}

instance Show Coord where
    show c = case (interval c) of
        Just (beg,end) -> rangeStr (beg,end)
        Nothing -> ref c
        where
            rangeStr (b, e)
                | b == e = ref c ++ ":" ++ (show b)
                | otherwise = ref c ++ ":" ++ (show b) ++ "-" ++ (show e) 

parseCoord :: String -> Coord
parseCoord s = Coord s Nothing

data ITree = Leaf (Integer, Integer) | Node ITree ITree

instance Show ITree where
    show t = "<interval tree>"
