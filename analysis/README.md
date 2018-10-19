# Instructions

Plot pre, b only, and s+b fits:
```
root -l
.L make_plots.C++
make_plots()
```

Compare pre and post fit nuisances:
```
python ../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py ../model/fitDiagnostics.root --all
```
