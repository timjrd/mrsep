module Main where

import Input
import Test

param = Param
  { known           = 4
  , unknown         = 2
  , alleleLength    = 100
  , mutations       = 1
  , knownInSample   = 2
  , unknownInSample = 1
  , coverage        = 2
  , minReadLength   = 2 }

main = tests 1000 param >>= putStr . showAnomaly
