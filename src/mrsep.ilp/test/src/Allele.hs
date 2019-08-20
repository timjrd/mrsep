module Allele where

import Control.Monad
import Control.Monad.Random

import System.Random.Shuffle

import Data.List
import qualified Data.Set    as S
import qualified Data.Vector as V

data Base = A | T | C | G
  deriving (Eq,Ord,Show)

type Allele = V.Vector Base

randomBase :: MonadRandom m => m Base
randomBase = do
  r <- getRandomR (1,4::Int)
  return $ case r of
    1 -> A
    2 -> T
    3 -> C
    _ -> G

randomAlleles :: MonadRandom m => Int -> Int -> Int -> m [Allele]
randomAlleles m n mutations =
  V.replicateM n randomBase >>= f S.empty
  where
    f alleles allele
      | S.size alleles >= m = shuffleM $ S.toList alleles
      | otherwise = do
          indexes <- replicateM mutations $ getRandomR (0,n-1)
          bases   <- replicateM mutations randomBase
          let newAllele = allele V.// (zip indexes bases)
          f (S.insert newAllele alleles) newAllele

variablePositions :: [Allele] -> [Allele]
variablePositions [] = []
variablePositions alleles@(first:_) =
  let const = foldl' (\bin x -> V.zipWith (&&) bin $ V.zipWith (==) first x)
              (V.replicate (V.length first) True) alleles
  in map (V.ifilter $ \i _ -> not $ const V.! i) alleles

showAlleles :: [Allele] -> String
showAlleles = unlines . map (concatMap show . V.toList)
