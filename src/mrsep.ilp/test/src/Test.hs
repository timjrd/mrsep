module Test where

import Control.Exception
import Control.Concurrent.Async

import Input
import Run
import Solve

data Anomaly = Anomaly Input Solution [Solution]
  deriving Show

instance Exception Anomaly

test :: Param -> IO (Maybe Anomaly)
test param = do
  input <- mkInput param
  let sinput = prepareInput input
      solutions = solve sinput
  solution <- run input
  return $
    if isValid sinput solution
    && solution `elem` solutions
    then Nothing
    else Just $ Anomaly input solution solutions
  
tests :: Int -> Param -> IO (Maybe Anomaly)
tests n param = do
  a <- try
    $ replicateConcurrently_ n
    $ test param >>= maybe (return ()) throw
  return $ case a of
    Left  r  -> Just r
    Right () -> Nothing

showAnomaly :: Maybe Anomaly -> String
showAnomaly Nothing = "No anomalies were detected.\n"
showAnomaly (Just (Anomaly input solution solutions)) = concat
  [( paragraph 45 $ unwords $
     ["At least the following anomaly were detected:", valid, optimal] )
  , "\n\nINPUT:\n\n", showInput input
  , "\nEXTERNAL SOLUTION:\n\n", showSolution solution
  , concatMap showInternal $ zip [1..] solutions ]
  where
    sinput = prepareInput input
    valid = if isValid sinput solution
      then "The external solution is valid."
      else "The external solution is NOT valid."
    optimal = if null solutions
      then "The internal solver did not found any solutions, this is an ERROR."
      else if solution `elem` solutions
      then "According to the internal solver, the external solution is optimal."
      else "According to the internal solver, the external solution is NOT optimal."
    showInternal (i,s) =
      "\nINTERNAL SOLUTION #" ++ show i ++ ":\n\n" ++ showSolution s

paragraph :: Int -> String -> String
paragraph n = f 0 . words
  where
    f l (w:ws)
      | l + length w + 1 > n = "\n" ++ f 0 (w:ws)
      | otherwise = w ++ " " ++ f (l + length w + 1) ws
    f _ [] = []

