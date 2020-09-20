(*Compile*)
lorenz = Compile[
  {{parameterList, _Real, 1}},
  Table[
   counter = 0;
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
     Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
     StartingStepSize -> 0.01,
     MaxSteps -> 10100,
     PrecisionGoal -> 8,
     AccuracyGoal -> 8,
     StepMonitor :> counter++
     ];
   {First@Evaluate[x[t] /. sol] /. t -> 10, counter},
   {i, 1, Length[parameterList]}],
  CompilationTarget -> "C",
  "RuntimeOptions" -> "Speed"
  ]

(*Run and Benchmark*)
nrOfParameters = 128;
parameterList = N@Range[0, 50, 50/(nrOfParameters - 1)];
runtime = Timing[endvals = lorenz[parameterList]];
runtime[[1]]
Print[runtime[[1]]]

(*Plot Endvals*)
data = Transpose[{parameterList, endvals[[All, 1]], endvals[[All, 2]]}];
plt = ListPlot[
  data[[All, {1, 2}]],
   PlotStyle -> Red,
   FrameLabel ->  {"p", "x(10)"},
   PlotLabel -> Style["lorenz végértékek - NDSolve[]", Black, 14]
   ]
