module Solve where

import Data.Maybe
import Data.List
import qualified Data.IntMap as M
import qualified Data.Vector as V

import Allele
import Read
import Mapping
import Input
import Run

type SInput = (V.Vector Allele, M.IntMap Mapping)
type Validity = M.IntMap Int
data Validity' = Invalid | MaybeValid Validity | Valid Validity

prepareInput :: Input -> SInput
prepareInput (alleles,mappings) =
  (V.fromList alleles, foldr f M.empty mappings)
  where f m@(Mapping i q _ _ _) ms = M.insert (pair i q) m ms

solve :: SInput -> [Solution]
solve input@(alleles,mappings) =
  let
    solutions = f (1,1,1) M.empty
    optimal = maximum $ map length solutions
  in
    sort
    $ map sort
    $ filter ((==optimal) . length) solutions
  where
    f assignment@(i,q,j) x
      | pair i 1 `M.notMember` mappings = []
      | pair i q `M.notMember` mappings = f (i+1, 1, 1) x
      | j > V.length alleles = f (i, q+1, 1) x
      | otherwise = case validity input assignment x of
          Invalid      -> next
          Valid      y -> [assignment] : add y ++ next
          MaybeValid y -> add y ++ next
      where
        next = f (i, q, j+1) x
        add y = map (assignment:) (f (i+1, 1, 1) y)

validity :: SInput -> Assignment -> Validity -> Validity'
validity (alleles,mappings) (i,q,j) x =
  let
    (Mapping _ _ fromK toK read) = mappings M.! (pair i q)
    match = read == V.toList (V.slice (fromK-1) (toK-fromK+1) (alleles V.! (j-1)))
    y = foldr (\k -> M.insertWith (+) (pair j k) 1) x [fromK..toK]
    maxJ = V.length alleles
    maxK = V.length $ V.head alleles
    coverages = [ [M.findWithDefault 0 (pair j' k) y | k <- [1..maxK]]
                | j' <- [1..maxJ] ]
    uniform = and $ map (\(c:cs) -> all (==c) cs) coverages
  in
    if      not match then Invalid
    else if uniform   then Valid      y
    else                   MaybeValid y

validity' :: SInput -> Assignment -> Validity' -> Validity'
validity' a b  Invalid       = Invalid
validity' a b (MaybeValid c) = validity a b c
validity' a b (Valid      c) = validity a b c

isValid :: SInput -> Solution -> Bool
isValid input solution =
  case foldr (validity' input) (MaybeValid M.empty) solution of
    Valid _ -> True
    _       -> False

pair :: Int -> Int -> Int
pair k1 k2 = (k1 + k2) * (k1 + k2 + 1) `div` 2 + k2
