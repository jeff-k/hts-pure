module Bio.Data.Location (mkPos,Pos,coverage,interval,ref) where

data Pos = Pos {ref::Int, interval::(Integer, Integer)}

instance Show Pos where
  show c = case (interval c) of
      (0,0) -> show $ ref c
      (beg,end) -> rangeStr (beg,end)
      where
          rangeStr (b, e)
              | b == e = show (ref c) ++ ":" ++ (show b)
              | otherwise = show (ref c) ++ ":" ++ (show b) ++ "-" ++ (show e) 

coverage :: [Pos] -> [Pos]
coverage = undefined

mkPos :: Int -> Integer -> Integer -> Pos
mkPos r b e = Pos r (b,e)
