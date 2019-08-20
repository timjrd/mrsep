module Mapping where

import Control.Monad.Random

import System.Random.Shuffle

import Data.List
import qualified Data.Vector as V

import Allele
import qualified Read as R

data Mapping = Mapping Int Int Int Int R.Read
  deriving Show

readMappings :: MonadRandom m => [Allele] -> [R.Read] -> m [Mapping]
readMappings alleles reads = shuffleM
  $ concat
  $ zipWith m' [1..]
  $ filter (not . null)
  $ map r reads
  where
    m' i = map $ \(Mapping _ a b c d) -> (Mapping i a b c d)
    r read = zipWith m [1..] $ nub $ concatMap a alleles
      where
        n = length read
        m q k = Mapping undefined q (k+1) (k+n) read
        a allele = elemIndices read slices
          where
            slices = map slice [0 .. V.length allele - n]
            slice k = V.toList $ V.slice k n allele

showMappings :: [Mapping] -> String
showMappings = unlines . map f
  where
    f (Mapping i q fromK toK read) = unwords $
      show i : show q : concat (zipWith (\b k -> [show b, show k]) read [fromK..])
