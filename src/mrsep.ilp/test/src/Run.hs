module Run where

import System.Environment
import System.Process

import Data.List

import Input

type Assignment = (Int,Int,Int)
type Solution = [Assignment]

run :: Input -> IO Solution
run input = do
  bin <- getEnv "SOLVER"
  out <- readCreateProcess (proc bin ["0","-"]) {std_err=CreatePipe}
    $ showInput input
  return
    $ sort
    $ map ((\[a,b,c]->(a,b,c)) . map read . words)
    $ tail
    $ lines out

showSolution :: Solution -> String
showSolution = unlines . map (\(a,b,c)-> unwords $ map show [a,b,c])
