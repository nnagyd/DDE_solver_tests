(*Definition of RK4*)

crkamat = {{1/2}, {0, 1/2}, {0, 0, 1}};
crkbvec = {1/6, 1/3, 1/3, 1/6};
crkcvec = {1/2, 1/2, 1};
ClassicalRungeKuttaCoefficients[4, p_] := N[{crkamat, crkbvec, crkcvec}, p];

(*Compile*)
logistic1 = Compile[
  {{parameterList, _Real, 1}},
  Table[
   sol = NDSolve[
     {
      x'[t] == x[t]*(parameterList[[i]] - x[t - 1]),
      x[t /; t <= 0] == 1.5 - Cos[t]
      },
     x[t],
     {t, 0, 10},
     Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 4,
       "Coefficients" -> ClassicalRungeKuttaCoefficients},
     StartingStepSize -> 0.001,
     MaxSteps -> 10100
     ];
   First@Evaluate[x[t] /. sol] /. t -> 10,
   {i, 1, Length[parameterList]}],
  CompilationTarget -> "C",
  "RuntimeOptions" -> "Speed"
  ]


(*Run and Benchmark*)
nrOfParameters = 32;
parameterList = N@Range[0, 4, 4/(nrOfParameters - 1)];
runtime = Timing[endvals = logistic1[parameterList]];
runtime[[1]]
Print[runtime[[1]]]

(*Plot Endvals*)
data = Transpose[{parameterList, endvals}];
plt = Show[
  ListPlot[
   data[[All, {1, 2}]],
   PlotStyle -> Red,
   Frame->True,
   FrameLabel -> {"p", "x(10)"},
   PlotLabel -> Style["logistic1 végértékek - NDSolve[]", Black, 14]
   ]
  ]
