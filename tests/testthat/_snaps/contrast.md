# print.survatr_result emits a stable banner

    Code
      print(res)
    Output
      <survatr_result>
        Type:        risk_difference
        Reference:   a0
        CI method:   none
        Time grid:   [1, 5] (5 unique times)
        Estimates:   10 rows
        Contrasts:   5 rows
      
         contrast  time    estimate    se ci_lower ci_upper
           <char> <int>       <num> <num>    <num>    <num>
      1: a1 vs a0     1 0.004163314    NA       NA       NA
      2: a1 vs a0     2 0.007510430    NA       NA       NA
      3: a1 vs a0     3 0.010161375    NA       NA       NA
      4: a1 vs a0     4 0.012220484    NA       NA       NA
      5: a1 vs a0     5 0.013778324    NA       NA       NA

