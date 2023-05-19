# General information

The purpose of this application is to provide further information based on the information produced by `gffcompare`. `gffcompare` provides comparison information on the accuracy of one or more `GFF`-files. Our interest has been in the further comparison of `GTF`-files.  

The vision is for this application is 

- to produce overview statistical data on the results provided by `gffcompare` 
- to provide further analysis on each transcript


## Statistical overview

TBA

## Further analysis - Offsets

### Defining offset

`gffcompare` output contains class codes for transcripts which provide information on the quality of the alignment. One interesting aspect in the case of misalignments is the amount of the offset in the misalignment. Here we provide a definition for offset and detail how offset is calculated in this application. 

Let $S$ be  a transcript from an analyzed read and assume that $S$ has one or more exons as children. Let the set of these exons be $E = \{e_1, \ldots , e_n\},\,n\in\mathbb{N},\,\mid E\mid \ge 1$. Let $R$ be a reference transcript to which $S$ has been compared to. Let exons children in $R$ be $X = \{x_1, \ldots, x_m\}$. 

Each exon is given as a tuple with a start index and end index $(a, b)\, a,b\in\mathbb{N},\,a < b$. Let $e_i\in E$ and $x_j\in X$ be an arbitrary exons and let $e_i = (a_e, b_e)$ and $x_j = (a_x, b_x)$.  

We assume that for a selected transcript the exon children do not overlap. That is, for arbitary two exons $e_i = (a_{e_i}, b_{e_i}) $ and $e_k = (a_{e_k}, b_{e_k})$ it always holds that either $b_{e_i} < a_{e_k}$ or $b_{e_k} < a_{e_i}$

**Definition 1.** For arbitary $e_i$ and $x_j$ The distance $a_{e_i} - a_{x_j}$ is the offset at the start index and the distance $b_{e_i} - b_{x_j}$ is the offset at the end index. Notice that if the offset is 'to the left', the output value is negative and similarily if the offset is 'to the right', the output value is positive.   

The total offset $t_{ij}$ for $e_i$ and $x_j$ is the sum off the absolute value of these two offsets

$$ t_{ij} =  \mid a_{e_i} - a_{x_j} \mid + \mid b_{e_i} - b_{x_j} \mid $$

**Definition 2.** For any $e_i$ we call the optimal match in $X$ to be the exon $x_j$ that satisfies the following conditions: 

1. The offset for $e_i$ compared to $x_j$ is the smallest possible. That is for $e_i$ and $X=[x_1, \ldots x_m]$ $\min \{t_{i1}, t_{i2}, \ldots, t_{im} \} = t_{ij}$. 
2. No other exon in $E$ has a smaller offset compared to $x_j$. That is $\not\exists e_k$ for which $\min { t_{i1}, t_{i2}, \ldots, t_{im} } = t_{kj}$ and $t_{kj} < t_{ij}$


### Computing offset

The basic idea in computing the offset is to have the lists of tuples $E$ and $X$ in ascending order. The list $E$ is then iterated over to find the optimal matches from list $X$ for the exons. At each index, if possible, the items in the next index are also considered for an optimal match. Let $e_i\in E$ and $x_j\in X$ be arbitary items. Now four possible scenarios can happen

1. $e_i$ and $x_j$ are an optimal match
2. $e_{i}$ and $x_{j+1}$ are an optimal match
3. $e_{i+1}$ and $x_{j}$ are an optimal match
4. $e_{i+1}$ and $x_{j+1}$ are an optimal match

Each of these cases needs to be considered and in case of $e_{i+1}$ being an item in the optimal match, $e_{i+2}$ needs to be considered on the next iteration. 

The following pseudocode computes the offsets following the rules given in definition two:

```
function compute_offset(list E, list X)
  list_of_results = []
  x_start_index = 0
  for e_index in list E:
    result  = (inf, inf)
    for x_index in X from x_start_index to len(X):
      total_offset = compute total offset for E[e_index] and X[x_index] 
      if e_index < len(E) - 1:
        total_offset_next = compute total offset for E[e_index + 1] and X[x_index]
        if total_offset_next < total_offset:
          if x_index < len(X) - 1:
            total_offset_next_ref = compute total offset for E[e_index + 1] and X[x_index + 1]
            if total_offset_next < total_offset_next_ref:
              x_start_index = x_index
              break inner loop
      if total_offset < total offset of result:
        if total offset of result != (inf, inf):
          add (-inf, -inf) to list_of_results
        result = (start_offset, end_offset) for E[e_index] and X[x_index]
        x_start_index = x_index + 1
      else:
        break inner for loop
  add result to list_of_results
```

Assume that we have a list of exons $E = [e_1, \ldots, e_n]$ and a list of reference exons $X=[x_1, \ldots , x_m]$. The algorithm starts from $e_1$ and iterates through the exons. Let us assume that the algorithm is now at an arbitary index $e_i$. Additionally assume that we are inspecting arbitary index $x_j\in X$. The following steps happen:

1. set result to $(\inf, \inf)$
2. iterate through reference exons starting from the current x_start_index to the end of the reference exons. 
3. if the total offset between $e_i$ and arbitary $x_j$ is smaller than the total offset of the value stored in result
    - If $e_p$ is not the last exon, check whether the total offset between $e_{p+1}$ and $x_i$ is smaller that the total offset between $e_i$ and $x_j$. If it is, check further if the offset between $e_{i+1}$ and $x_{j+1}$ is smaller that the offset between $e_{i+1}$ and $x_j$. If it is, update the x_start_index to be the current index and break the inner loop
    - If the new total offset is smaller than the offset in result, check whether the values in result are less than $(\inf, \inf)$. In that case, append $(-\inf, -\inf)$ to the list of results. Update result to the new offset values. Set the x_start_index = x_index + 1. Incase the new offset is not smaller than the offset stored in result, break the inner loop. 
4. at the end append result to the list of results. 


### Offset output

1. In case of a match the offset is expressed from the point of view of the analyzed transcript in the form of a tuple of two integers $(a, b),\,a,b\in\mathbb{Z}$. A negative integer indicates that the exon $e_i$ from the analyzed data has a smaller value than the matching reference exon $x_j$ and similarily a positive value indicates that the exon $e_i$ has a higher value
2. A tuple $(\inf, \inf)$ indicates that no optimal match for an exon in analyzed data was found (i.e. there's possibly an exon in the analyzed data that is not present in the reference data)
3. A tuple $(-\inf, -\inf)$ that no optimal match for an exon in the reference data was found (i.e. there's possibly an exon is missing from the analyzed data)
