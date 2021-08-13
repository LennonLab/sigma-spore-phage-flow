Sporulation assay with IPTG-induced sigma factors
================

**Only consider counts with at least 100 events**

Quantities based on lower event counts are designated a value of 1
cell/mL, to prevent problems with ratio and logs.

# Replication outliers

I have N=3 for each flask at each time point. There are some
measurements which are obviously way off. To get rid of these I will
choose from each triplicate the 2 points that are in best agreement and
remove the 3rd point. There is a function to idenify such a point in
pacakge ‘outliers’ &gt; outlier {outliers} Finds value with largest
difference between it and sample mean, which can be an outlier.

> logical: if set to TRUE, gives vector of logical values, and possible
> outlier position is marked by TRUE

I will apply the outlier filtering on the number of total cells
(veg+spore)

# Overview of results

Concentrations of cell types:
![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
Remove blanks

Change in respnse to induction, log10(IPTG/noIPTG):

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony', 'pop'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

\#stats on log ratio of % spores

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

``` r
summary(aov(induction.ratio.log~strain+colony+exp, d.sum))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## strain       6 117.50  19.583  68.271 <2e-16 ***
    ## colony       4   0.51   0.129   0.448 0.7732    
    ## exp          3   2.33   0.778   2.713 0.0557 .  
    ## Residuals   46  13.19   0.287                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
   d.test <- tibble()
# t-tests against control
   for (s in unique(d.sum$strain)){
      if (s =="pDR110") next
      d.test <-d.sum %>% 
         filter(strain == s |strain =="pDR110") %>% 
         t.test(induction.ratio.log ~ strain, data = .) %>% 
         broom::tidy() %>% 
         mutate(strain = s) %>% 
         bind_rows(d.test, .)
   }
   
# adjust p-value for multiple testing
   d.test <- d.test %>% 
      mutate(adj.p = p.adjust(p.value, method = "BH"),
             p.lab = stars.pval(adj.p)) %>% 
      relocate(strain, p.value, adj.p, p.lab)
   
   knitr::kable(d.test, format = "html")
```

<table>
<thead>
<tr>
<th style="text-align:left;">
strain
</th>
<th style="text-align:right;">
p.value
</th>
<th style="text-align:right;">
adj.p
</th>
<th style="text-align:left;">
p.lab
</th>
<th style="text-align:right;">
estimate
</th>
<th style="text-align:right;">
estimate1
</th>
<th style="text-align:right;">
estimate2
</th>
<th style="text-align:right;">
statistic
</th>
<th style="text-align:right;">
parameter
</th>
<th style="text-align:right;">
conf.low
</th>
<th style="text-align:right;">
conf.high
</th>
<th style="text-align:left;">
method
</th>
<th style="text-align:left;">
alternative
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
sigF
</td>
<td style="text-align:right;">
0.0000227
</td>
<td style="text-align:right;">
0.0000680
</td>
<td style="text-align:left;">
\*\*\*
</td>
<td style="text-align:right;">
2.0110228
</td>
<td style="text-align:right;">
-0.0422284
</td>
<td style="text-align:right;">
-2.0532513
</td>
<td style="text-align:right;">
6.9363013
</td>
<td style="text-align:right;">
11.182110
</td>
<td style="text-align:right;">
1.3741629
</td>
<td style="text-align:right;">
2.6478828
</td>
<td style="text-align:left;">
Welch Two Sample t-test
</td>
<td style="text-align:left;">
two.sided
</td>
</tr>
<tr>
<td style="text-align:left;">
sigG
</td>
<td style="text-align:right;">
0.0160111
</td>
<td style="text-align:right;">
0.0240167
</td>
<td style="text-align:left;">

-   </td>
    <td style="text-align:right;">
    0.9974987
    </td>
    <td style="text-align:right;">
    -0.0422284
    </td>
    <td style="text-align:right;">
    -1.0397271
    </td>
    <td style="text-align:right;">
    2.9005752
    </td>
    <td style="text-align:right;">
    9.877469
    </td>
    <td style="text-align:right;">
    0.2299583
    </td>
    <td style="text-align:right;">
    1.7650390
    </td>
    <td style="text-align:left;">
    Welch Two Sample t-test
    </td>
    <td style="text-align:left;">
    two.sided
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    Goe3
    </td>
    <td style="text-align:right;">
    0.0004348
    </td>
    <td style="text-align:right;">
    0.0008697
    </td>
    <td style="text-align:left;">
    \*\*\*
    </td>
    <td style="text-align:right;">
    0.8739457
    </td>
    <td style="text-align:right;">
    -0.0422284
    </td>
    <td style="text-align:right;">
    -0.9161741
    </td>
    <td style="text-align:right;">
    4.6073686
    </td>
    <td style="text-align:right;">
    13.635705
    </td>
    <td style="text-align:right;">
    0.4660912
    </td>
    <td style="text-align:right;">
    1.2818002
    </td>
    <td style="text-align:left;">
    Welch Two Sample t-test
    </td>
    <td style="text-align:left;">
    two.sided
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    ELDg168
    </td>
    <td style="text-align:right;">
    0.0457281
    </td>
    <td style="text-align:right;">
    0.0548737
    </td>
    <td style="text-align:left;">
    .
    </td>
    <td style="text-align:right;">
    0.5380198
    </td>
    <td style="text-align:right;">
    -0.0422284
    </td>
    <td style="text-align:right;">
    -0.5802482
    </td>
    <td style="text-align:right;">
    2.1731837
    </td>
    <td style="text-align:right;">
    15.422322
    </td>
    <td style="text-align:right;">
    0.0115880
    </td>
    <td style="text-align:right;">
    1.0644515
    </td>
    <td style="text-align:left;">
    Welch Two Sample t-test
    </td>
    <td style="text-align:left;">
    two.sided
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    SP10
    </td>
    <td style="text-align:right;">
    0.8927637
    </td>
    <td style="text-align:right;">
    0.8927637
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:right;">
    0.0229639
    </td>
    <td style="text-align:right;">
    -0.0422284
    </td>
    <td style="text-align:right;">
    -0.0651923
    </td>
    <td style="text-align:right;">
    0.1379307
    </td>
    <td style="text-align:right;">
    11.108165
    </td>
    <td style="text-align:right;">
    -0.3430398
    </td>
    <td style="text-align:right;">
    0.3889676
    </td>
    <td style="text-align:left;">
    Welch Two Sample t-test
    </td>
    <td style="text-align:left;">
    two.sided
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    ELDg169
    </td>
    <td style="text-align:right;">
    0.0000000
    </td>
    <td style="text-align:right;">
    0.0000000
    </td>
    <td style="text-align:left;">
    \*\*\*
    </td>
    <td style="text-align:right;">
    4.1088530
    </td>
    <td style="text-align:right;">
    -0.0422284
    </td>
    <td style="text-align:right;">
    -4.1510814
    </td>
    <td style="text-align:right;">
    22.8371752
    </td>
    <td style="text-align:right;">
    13.586504
    </td>
    <td style="text-align:right;">
    3.7218594
    </td>
    <td style="text-align:right;">
    4.4958467
    </td>
    <td style="text-align:left;">
    Welch Two Sample t-test
    </td>
    <td style="text-align:left;">
    two.sided
    </td>
    </tr>
    </tbody>
    </table>

<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
