module Read where

import Prelude hiding (Read)

import Control.Monad.Random

import System.Random.Shuffle

import qualified Data.Vector as V

import Allele

type Read = [Base]

reads :: MonadRandom m => Int -> Int -> [Allele] -> m [Read]
reads coverage minLength sample =
  mapM (replicateM coverage . reads1 minLength) sample
  >>= shuffleM . concat . concat

reads1 :: MonadRandom m => Int -> Allele -> m [Read]
reads1 minLength bases
  | V.null bases = return []
  | otherwise = do
      let maxLength = V.length bases
      n <- getRandomR (minLength, maxLength)
      if maxLength - n < minLength
        then return [V.toList bases]
        else do let (read,remainder) = V.splitAt n bases
                nextReads <- reads1 minLength remainder
                return $ V.toList read : nextReads
