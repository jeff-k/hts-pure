module Bio.Data.Location (Pos, parsePos, coverage, interval, t) where

data Pos = Pos {ref::String, interval::Maybe (Int, Int)}

instance Show Pos where
    show c = case (interval c) of
        Just (beg,end) -> rangeStr (beg,end)
        Nothing -> ref c
        where
            rangeStr (b, e)
                | b == e = ref c ++ ":" ++ (show b)
                | otherwise = ref c ++ ":" ++ (show b) ++ "-" ++ (show e) 

t :: Pos
t = Pos "chr1" (Just (50000000, 100000000))

coverage :: [Pos] -> [Pos]
coverage = undefined

parsePos :: String -> Pos
parsePos = undefined
