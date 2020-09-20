(*Definition of RK4*)
crkamat = {{1/2}, {0, 1/2}, {0, 0, 1}};
crkbvec = {1/6, 1/3, 1/3, 1/6};
crkcvec = {1/2, 1/2, 1};
ClassicalRungeKuttaCoefficients[4, p_] := N[{crkamat, crkbvec, crkcvec}, p];

(*Compile*)
lorenz = Compile[
  {{parameterList, _Real, 1}},
  Table[
   sol = NDSolve[
     {
      x'[t] == 10 (y[t - 13/28] - x[t]),
      y'[t] == parameterList[[i]]*x[t] - y[t] - z[t]*x[t],
      z'[t] == x[t]*y[t] - 8/3  z[t],
      x[t /; t <= 0] == -8,
      y[t /; t <= 0] == -8 + Sin[2 \[Pi] t],
      z[t /; t <= 0] == -8
      },
     x[t],
     {t, 0, 10},
     Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 4, "Coefficients" -> ClassicalRungeKuttaCoefficients},
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
parameterList = N@Range[0, 50, 50/(nrOfParameters - 1)];
runtime = Timing[endvals = lorenz[parameterList]];
runtime[[1]]
Print[runtime[[1]]]

(*Save Endvals*)
data = Transpose[{parameterList, endvals}];

(*Plot Endvals*)
plt = Show[
  ListPlot[
   data[[All, {1, 2}]],
   PlotStyle -> Red,
   Frame->True,
   FrameLabel -> {"p", "x(10)"},
   PlotLabel -> Style["lorenz végértékek - NDSolve[]", Black, 14]
   ]
  ]
