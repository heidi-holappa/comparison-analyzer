# General information

The purpose of this application is to provide further information based on the information produced by `gffcompare`. `gffcompare` provides comparison information on the accuracy of one or more `GFF`-files. Our interest has been in the further comparison of `GTF`-files.  

The vision is for this application is 

- to produce overview statistical data on the results provided by `gffcompare` 
- to provide further analysis on each transcript


## Statistical overview

TBA

## Further analysis - Offsets

### Defining offset

`gffcompare` classifies comparison results with class codes providing detailed information on the comparison results. One interesting aspect in the case of misalignments is the nature of the offset in the misalignment. Here we provide a definition for offset and a detail how offset is computed. 

Let $S$ be  a transcript from an analyzed read and assume that $S$ has one or more exons as children. Let the set of these exons be $E = \{e_1, \ldots , e_n\},\,n\in\mathbb{N},\,\mid E\mid \ge 1$. Let $R$ be a reference transcript to which $S$ has been compared to. Let exons children in $R$ be $X = \{x_1, \ldots, x_m\}$. 

Each exon is given as a tuple with a start index and end index $(i, j)\, i,j\in\mathbb{N},\,i < j$. Let $e_i\in E$ and $x_j\in X$ be an arbitrary exons and let $e_i = (i_e, j_e)$ and $x_j = (i_x, j_x)$.  

We assume that for a selected transcript the exon children do not overlap. That is, for arbitary two exons $e_i = (a, b) $ and $e_k = (c, d)$ it always holds that either $b < c$ or $d < a$

**Definition 1.** The distance $i_e - i_x$ is the offset at the start index and the distance $j_e - j_x$ is the offset at the end index. Notice that if the offset is 'to the left', the output value is negative and similarily if the offset is 'to the right', the output value is positive.   

The total offset is the sum off the absolute value of these two offsets

$$ \mid i_e - i_x \mid + \mid j_e - j_x \mid $$

**Definition 2.** For any $e_i$ we call the optimal match in $X$ to be the exon $x_j$ that satisfies the following conditions: 

1. The offset for $e_i$ compared to $x_j$ is the smallest possible
2. No other exon in $E$ has a smaller offset compared to $x_j$


### Computing offset

The following pseudocode computes the offset following the rules given in definition two

```
function compute_offset(list E, list X)
  x_start_index = 0
  for e_index in list E:
    result  = (inf, inf)
    for x_index in X from x_start_index to len(X):
      total_offset = compute total offset for E[e_index] and X[x_index] 
      total_offset_next = compute total offset for E[e_index + 1] and X[x_index]
      if e_index < len(E) - 1 and total_offset_next < total_offset:
          x_start_index = x_index
          continue
      if total_offset < total offset of result:
        result = (start_offset, end_offset) for E[e_index] and X[x_index]
        x_start_index = x_index + 1
      else:
        continue
```

The algorithm is quite simple and perhaps self-explanatory. Assume that we have a list of exons $E = [e_1, \ldots, e_n]$ and a list of reference exons $X=[x_1, \ldots , x_m]$. The algorithm starts from $e_1$ and iterates through the exons. Let us assume that the algorithm is now at an arbitary index $e_p$ and the next index would be $e_q$. The following steps happen:

1. set result to $(\inf, \inf)$
2. iterate through reference exons starting from the current x_start_index to the end of the reference exons. 
3. if the total offset between $e_p$ and arbitary $x_i$ is smaller than the total offset of the result
    - If $e_p$ is not the last exon, check whether the total offset between $e_q$ and $x_i$ is smaller that the total offset between $e_p$ and $x_i$. If it is, update the x_start_index to be the current index and continue to iterate over $e_q$
    - If the new total offset is smaller than the offset in result, update the result and set the x_start_index to be one index after the current index
    - if not, continue to iterate over $e_q$