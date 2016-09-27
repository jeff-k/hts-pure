module Bio.Data.Location (Pos,coverage,interval,ref) where

data Pos = Pos {ref::String, interval::(Integer, Integer)}

instance Show Pos where
  show c = case (interval c) of
      (0,0) -> ref c
      (beg,end) -> rangeStr (beg,end)
      where
          rangeStr (b, e)
              | b == e = ref c ++ ":" ++ (show b)
              | otherwise = ref c ++ ":" ++ (show b) ++ "-" ++ (show e) 

--instance Read Pos where
--  read p = undefined

coverage :: [Pos] -> [Pos]
coverage = undefined

