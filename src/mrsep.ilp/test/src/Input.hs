module Input where

import Control.Monad.Random

import System.Random.Shuffle

import qualified Data.Vector as V

import Allele
import qualified Read as R
import Mapping

data Param = Param
  { known           :: Int
  , unknown         :: Int
  , alleleLength    :: Int
  , mutations       :: Int
  , knownInSample   :: Int
  , unknownInSample :: Int
  , coverage        :: Int
  , minReadLength   :: Int }

type Input = ([Allele],[Mapping])

mkInput :: MonadRandom m => Param -> m Input
mkInput p = do
  alleles <- variablePositions <$>
    randomAlleles (known p + unknown p) (alleleLength p) (mutations p)
  let (knownAlleles,unknownAlleles) = splitAt (known p) alleles
  ks <- replicateM (knownInSample   p) $ getRandomR (0, known   p - 1)
  us <- replicateM (unknownInSample p) $ getRandomR (0, unknown p - 1)
  sample <- shuffleM $
    zipWith (!!) (repeat knownAlleles  ) ks ++
    zipWith (!!) (repeat unknownAlleles) us
  reads <- R.reads (coverage p) (minReadLength p) sample
  mappings <- readMappings knownAlleles reads
  return (knownAlleles,mappings)

showInput :: Input -> String
showInput (as,ms) = showAlleles as ++ "\n" ++ showMappings ms

someInput :: Input
someInput =
  ([ v [A,T,C,G,A,T]
   , v [T,A,C,G,A,A]
   , v [T,A,G,C,T,A] ]

  ,[ m 1 1 3 5 [C,G,A]
   , m 2 1 1 2 [A,T]
   , m 2 2 5 6 [A,T]
   , m 3 1 3 6 [C,G,A,T]
   , m 4 1 1 1 [A]
   , m 4 2 5 5 [A]
   , m 4 3 2 2 [A]
   , m 4 4 6 6 [A]
   , m 5 1 1 2 [T,A]
   , m 5 2 5 6 [T,A]
   , m 6 1 1 2 [T,A]
   , m 6 2 5 6 [T,A]
   , m 7 1 1 2 [T,A]
   , m 7 2 5 6 [T,A]
   , m 8 1 3 4 [G,C] ])
  where
    v = V.fromList
    m = Mapping
    
