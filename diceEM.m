 (* This defines the accessors for the parameter set of the dice model. It is stored as a list
    {type1Prob, type2Prob, faceProbs1, faceProbs2}. The accessors make the code more readable. *)
 type1Prob[parameters_] := parameters[[1]];
 type2Prob[parameters_] := parameters[[2]];
 faceProbs1[parameters_] := parameters[[3]];
 faceProbs2[parameters_] := parameters[[4]];
 
 (* diceEM just sets up for the actual EM algorithm by counting the frequencies of the faces
    in each trial and providing initial random or arbitrary values to diceEMIterator.*)
 (* A trial is the result of drawing a die at random from the bag and rolling it n times. *)
 diceEM[trials_, maxIterations_, accuracy_, randomSeed_:314, debug_:False]:=
	Module[{numFaces, binCountsList, initialFaceProbs1, initialFaceProbs2},
		(*SeedRandom[] initializes the random number generator so you can get deterministic
		  results for testing purposes.*)
		SeedRandom[randomSeed];
		(* For strictly defined EM, which uses maximum likelihood estimations, any die faces
		   that do not occur anywhere in the input can be treated as though they don't exist.*)
		numFaces = Max[trials];
		(*initialize binCountsList here, by mapping the builtin function BinCounts over the trials.
		  Each entry in the resulting binCountsList should be a list of the face frequencies from one trial.
		  E.g., if there are four faces on the dice, one entry might be {3, 0, 1, 1}, indicating 
		  that on this trial face 1 was rolled 3 times, face 2 was not rolled at all, etc.  *)
	    (* Initialize initialFaceProbs1 and initialFaceProbs2 here. Use RandomReal and make it so
	       so that all entries are within a factor of two of one another. The built in functions
	       Normalize and Total may also be useful here. *)
	    binCountsList = Map[BinCounts[#, {Range[numFaces + 1]}] &, trials];
	    (*binCountsList sums up all the values in each sublist that are the integer values of Range[numFaces]. In binCount, numFaces+1
	    is used to make sure that the last bound is catched by the algorithm and thus the last face is counted*)
	    initialFaceProbs1 = Normalize[RandomReal[{1/numFaces,2/numFaces},numFaces], Total];
	    initialFaceProbs2 = Normalize[RandomReal[{1/numFaces,2/numFaces},numFaces], Total];
	    (*Used to randomly initialize the face probabilities, so that they all add to 1 and aren't very different from one another*)
	diceEMIterator[binCountsList, 
		           numFaces, 
		           {0.45, 0.55, initialFaceProbs1, initialFaceProbs2}, 
		           maxIterations, 
		           accuracy]]
		           
(* diceEMIterator implements the outer loop of the EM algorithm.
   It calls updateProbs on each iteration. *)		           
diceEMIterator[binCountsList_, numFaces_, initParams_, maxIterations_,
   accuracy_, debug_: False] := 
 Module[{oldParamEstimates, 
   newParamEstimates},
   (*Initialize the local variables.*)
   
  oldParamEstimates = initParams;
  
  Do[newParamEstimates = 
    updateProbs[binCountsList, oldParamEstimates]; 
   If[Total[Flatten[Abs[newParamEstimates - oldParamEstimates]]] < 
     accuracy, Break[], oldParamEstimates = newParamEstimates], 
   maxIterations];
   (*This Do loop continues to update newParamEstimates, for a total of maxIterations times. It stops before maxIterations
   if the total change in parameters is less than the assigned accuracy. Flatten is used to eliminate all sublists within the parameter
   estimates so I can nicely subtract them from one another*)
   
  
  If[type1Prob[newParamEstimates] <= type2Prob[newParamEstimates], 
   newParamEstimates, {type2Prob[newParamEstimates], 
    type1Prob[newParamEstimates], faceProbs2[newParamEstimates], 
    faceProbs1[newParamEstimates]}]
    
]
   
updateProbs[binCountsList_, oldParamEstimates_, debug_:False] :=
	Module[{posteriors,
		    (* type1Count and type2Count are the expected number of times a type1 or type2
		       die was drawn.*) 
		    type1Count, type2Count,
		    (* faceCounts1 is the expected number of times each face was rolled on a die 
		       of type 1.Likewise for faceCounts2.*) 
		    faceCounts1, faceCounts2},
		(*Create list of posterior probabilities of a Type1 die having been rolled on each draw 
		   by calling your dicePosteriors, which you should paste in to this file. *)
		posteriors = Map[dicePosterior[#,type1Prob[oldParamEstimates],type2Prob[oldParamEstimates],faceProbs1[oldParamEstimates],faceProbs2[oldParamEstimates]] &, binCountsList];
		(* Now use the posteriors to calculate EXPECTED number of times each die type was drawn. *) 
		type1Count = Total[posteriors];
		type2Count = Total[1-posteriors];
		(*These type 1 and type 2 counts are determined by finding the total posterior chances of picking dieTypeX. Since
		dicePosterior returns the probability of type1, you need to find the 1-posterior to get the probabilities for dieType2*)
		(* Now use the posteriors to calculate EXPECTED number of times each face was rolled
		   on each die typep. *) 
		faceCounts1 = Total[(binCountsList*posteriors)[[All,Range[Length[binCountsList[[1]]]]]]];
		faceCounts2 = Total[(binCountsList*(1-posteriors))[[All,Range[Length[binCountsList[[1]]]]]]];
		(*The expected face counts for each die are based on the multiplying respective posterior probabilities by the observed bins,
		and then adding the total count for each element within the sublists of the larger list to get the face totals for a given die *)
		{type1Count/(type1Count+type2Count),type2Count/(type1Count+type2Count),faceCounts1/(Total[faceCounts1]),faceCounts2/(Total[faceCounts2])}
		(*each respective faceProbability adds to one, which is why each element is divided by the total of that list. Additionally, type1Count
		and type2Count add to 1 which is why they're divided by the sum*)
		(* Finally, use these counts to compute maximum likelihood estimates for the parameters and 
		   return these estimates in a list: {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2} *)
	] 

 diceSample[numType1_, numType2_, type1_, type2_, draws_, rollsPerDraw_] :=
 Module[{ die1 = EmpiricalDistribution[type1 -> Range[4]], 
  die2 = EmpiricalDistribution[type2 -> Range[4]],
  dieList = 
   RandomChoice[{numType1/(numType1 + numType2), 
      numType2/(numType1 + numType2)} -> {0, 1}, draws] }, 
 Map[RandomVariate[If[# == 0, die1, die2], rollsPerDraw] &, dieList]
 ]
 
dicePosterior[binCounts_, type1Prior_, type2Prior_, faceProbs1_, faceProbs2_] :=
Module[{changeZeroes, getDieChance}, 
 (*initialize the two functions I define within parent function*)
 changeZeroes[faceProbs_] := 
  MapThread[If[#1 == 0, If[#2 == 0, 1, -1], #1] &, {faceProbs, binCounts}]; 
  (*If faceProbs contains 0, Return 1 if the outcome based on binCounts and faceProbs is possible and -1 if not*)
 getDieChance[faceProbs_, typePrior_] := 
  If[Min[faceProbs] < 0, 0, 
   typePrior*Apply[Times, Power[faceProbs, binCounts]]]; 
   (*Identify if it's an impossible situation if the minumum of faceProbs is -1 (-1<0); if impossible then return 0*)
 getDieChance[changeZeroes[faceProbs1],  
   type1Prior]/(getDieChance[changeZeroes[faceProbs1], type1Prior] + 
    getDieChance[changeZeroes[faceProbs2], type2Prior])]
    (*Here I call all my functions, using exhaustive conditionalizing to find bayes solution*)

myRound[x_, n_] :=
  N[IntegerPart[Round[x,10^-n]*10^n] / 10^n];
