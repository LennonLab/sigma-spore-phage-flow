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
![](analysis_post_FCM_files/figure-gfm/plot%20overview-1.png)<!-- -->
Remove blanks

Change in response to induction, log10(IPTG/noIPTG):

![](analysis_post_FCM_files/figure-gfm/plot%20response-1.png)<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony', 'pop'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

\#stats on log ratio of % spores

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

``` r
summary(aov(induction.ratio~strain+colony+exp, d.sum))
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## strain       6  7.421  1.2369  19.119 5.05e-11 ***
    ## colony       4  0.077  0.0193   0.299   0.8774    
    ## exp          3  0.662  0.2206   3.410   0.0251 *  
    ## Residuals   46  2.976  0.0647                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
   d.test <- tibble()
# t-tests against control
   for (s in unique(d.sum$strain)){
      if (s =="pDR110") next
      d.test <-d.sum %>% 
         filter(strain == s |strain =="pDR110") %>% 
         t.test(induction.ratio ~ strain, data = .) %>% 
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
0.0001818
</td>
<td style="text-align:right;">
0.0005453
</td>
<td style="text-align:left;">
\*\*\*
</td>
<td style="text-align:right;">
0.8652696
</td>
<td style="text-align:right;">
1.02375
</td>
<td style="text-align:right;">
0.1584802
</td>
<td style="text-align:right;">
6.4317579
</td>
<td style="text-align:right;">
8.198880
</td>
<td style="text-align:right;">
0.5563460
</td>
<td style="text-align:right;">
1.1741932
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
0.0151901
</td>
<td style="text-align:right;">
0.0227851
</td>
<td style="text-align:left;">

-   </td>
    <td style="text-align:right;">
    0.5331770
    </td>
    <td style="text-align:right;">
    1.02375
    </td>
    <td style="text-align:right;">
    0.4905729
    </td>
    <td style="text-align:right;">
    2.7689175
    </td>
    <td style="text-align:right;">
    13.857782
    </td>
    <td style="text-align:right;">
    0.1197833
    </td>
    <td style="text-align:right;">
    0.9465706
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
    0.0018176
    </td>
    <td style="text-align:right;">
    0.0036351
    </td>
    <td style="text-align:left;">
    \*\*
    </td>
    <td style="text-align:right;">
    0.6034826
    </td>
    <td style="text-align:right;">
    1.02375
    </td>
    <td style="text-align:right;">
    0.4202673
    </td>
    <td style="text-align:right;">
    4.3873280
    </td>
    <td style="text-align:right;">
    8.864083
    </td>
    <td style="text-align:right;">
    0.2915912
    </td>
    <td style="text-align:right;">
    0.9153740
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
    0.0532393
    </td>
    <td style="text-align:right;">
    0.0638871
    </td>
    <td style="text-align:left;">
    .
    </td>
    <td style="text-align:right;">
    0.3626490
    </td>
    <td style="text-align:right;">
    1.02375
    </td>
    <td style="text-align:right;">
    0.6611008
    </td>
    <td style="text-align:right;">
    2.0967041
    </td>
    <td style="text-align:right;">
    15.123347
    </td>
    <td style="text-align:right;">
    -0.0057479
    </td>
    <td style="text-align:right;">
    0.7310459
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
    0.6525419
    </td>
    <td style="text-align:right;">
    0.6525419
    </td>
    <td style="text-align:left;">
    </td>
    <td style="text-align:right;">
    0.0670483
    </td>
    <td style="text-align:right;">
    1.02375
    </td>
    <td style="text-align:right;">
    0.9567015
    </td>
    <td style="text-align:right;">
    0.4635731
    </td>
    <td style="text-align:right;">
    10.358384
    </td>
    <td style="text-align:right;">
    -0.2537108
    </td>
    <td style="text-align:right;">
    0.3878074
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
    0.0001067
    </td>
    <td style="text-align:right;">
    0.0005453
    </td>
    <td style="text-align:left;">
    \*\*\*
    </td>
    <td style="text-align:right;">
    1.0071715
    </td>
    <td style="text-align:right;">
    1.02375
    </td>
    <td style="text-align:right;">
    0.0165783
    </td>
    <td style="text-align:right;">
    7.8020680
    </td>
    <td style="text-align:right;">
    7.002654
    </td>
    <td style="text-align:right;">
    0.7019449
    </td>
    <td style="text-align:right;">
    1.3123982
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

![](analysis_post_FCM_files/figure-gfm/plot%20log-ratio-1.png)<!-- -->

# ms plot

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

    ## Warning: Ignoring unknown parameters: width

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

![](analysis_post_FCM_files/figure-gfm/plot%20log%20col-1.png)<!-- -->

    ## `summarise()` has grouped output by 'strain', 'treat', 'colony'. You can override using the `.groups` argument.

![](analysis_post_FCM_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
