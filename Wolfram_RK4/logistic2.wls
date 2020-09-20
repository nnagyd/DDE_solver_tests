(*Definition of RK4*)

crkamat = {{1/2}, {0, 1/2}, {0, 0, 1}};
crkbvec = {1/6, 1/3, 1/3, 1/6};
crkcvec = {1/2, 1/2, 1};
ClassicalRungeKuttaCoefficients[4, p_] := N[{crkamat, crkbvec, crkcvec}, p];

(*Compile*)
nonContInitMedium[t_] := Piecewise[{
   {Cos[4 Pi *t], -2 <= t <= -1.5},
   {t^2, -1.5 < t <= -Sqrt[2]},
   {Exp[t], -Sqrt[2] < t <= -1.101},
   {0, -1.101 < t <= -0.5},
   {t + 0.5, -0.5 < t }
   }]
logistic2 = Compile[
 {{parameterList, _Real, 1}},
 Table[
  sol = NDSolve[
    {
     x'[t] == x[t]*x[t - 2]*(parameterList[[i]] - x[t - 1]),
     x[t /; t <= 0] == nonContInitMedium[t]
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
parameterList = N@Range[0, 2, 2/(nrOfParameters - 1)];
runtime = Timing[endvals = logistic2[parameterList]];
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
   PlotLabel -> Style["logistic2 végértékek - NDSolve[]", Black, 14]
   ]
  ]
