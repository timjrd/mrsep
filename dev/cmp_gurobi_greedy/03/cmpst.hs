import Data.List

main = do
  (m,n) <- cmp
  putStrLn $ show m ++ " / " ++ show n ++ " = " ++ show (fromIntegral m / fromIntegral n)

strains = lines <$> readFile "cmp_4_8.txt"

cmp = do
  pairs <- transpose <$> strains
  return
    ( length $ filter same pairs
    , length pairs )
  where
    same (x:xs) = all (==x) xs
    
